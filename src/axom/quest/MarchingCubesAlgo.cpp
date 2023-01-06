// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/quest/MarchingCubesAlgo.hpp"
#include "axom/quest/detail/marching_cubes_lookup.hpp"
#include "axom/mint/mesh/UnstructuredMesh.hpp"
#include "conduit_blueprint.hpp"
#include "axom/fmt.hpp"

#ifndef __WHERE
#define __STRINGIZE(x) __STRINGIZE2(x)
#define __STRINGIZE2(x) #x
//!@brief String literal for code location
#define __WHERE __FILE__ ":" __STRINGIZE(__LINE__) "(" + std::string(__func__) + ") "
#endif

namespace axom
{
namespace quest
{

MarchingCubesAlgo::MarchingCubesAlgo(const conduit::Node &bpMesh,
                                     const std::string &coordsetName,
                                     const std::string &fcnField,
                                     const std::string &maskField)
  : _sd()
  , _ndim(0)
  , _coordsetPath("coordsets/" + coordsetName)
  , _fcnPath("fields/" + fcnField)
  , _maskPath(maskField.empty() ? std::string() : "fields/" + maskField)
  , _surfaceMesh(nullptr)
  , _cellIdField()
  , _domainIdField()
{
  _sd.reserve(conduit::blueprint::mesh::number_of_domains(bpMesh));
  for(auto &dom : bpMesh.children())
  {
    _sd.emplace_back(new MarchingCubesAlgo1(dom, coordsetName, fcnField, maskField));
    if(_ndim == 0)
    {
      _ndim = _sd.back()->dimension();
    }
    else
    {
      SLIC_ASSERT(_ndim == _sd.back()->dimension());
    }
  }
}


/*!
  @brief Set the output surface mesh object.
*/
void MarchingCubesAlgo::set_output_mesh(axom::mint::Mesh *surfaceMesh)
{
  _surfaceMesh = surfaceMesh;

  for(auto &s : _sd)
  {
    s->set_output_mesh(_surfaceMesh);
  }
}


void MarchingCubesAlgo::compute_iso_surface(double contourVal)
{
  SLIC_ASSERT_MSG(_surfaceMesh,
                  "You must call set_output_mesh before compute_iso_surface.");

  if(!_domainIdField.empty() && !_surfaceMesh->hasField(_domainIdField, axom::mint::CELL_CENTERED))
  {
    _surfaceMesh->createField<int>(_domainIdField, axom::mint::CELL_CENTERED);
  }

  for(int dId=0; dId<_sd.size(); ++dId)
  {
    std::shared_ptr<MarchingCubesAlgo1>& single = _sd[dId];

    auto nPrev = _surfaceMesh->getNumberOfCells();
    single->compute_iso_surface(contourVal);
    auto nNew = _surfaceMesh->getNumberOfCells();

    if(nNew > nPrev && !_domainIdField.empty())
    {
      auto *domainIdPtr =
        _surfaceMesh->getFieldPtr<int>(_domainIdField, axom::mint::CELL_CENTERED);
      for(int n=nPrev; n<nNew; ++n)
      {
        domainIdPtr[n] = dId;
      }
    }
  }
}

// From app
#define NDSET3D(v,v1,v2,v3,v4,v5,v6,v7,v8)  \
   v4 = v ;   \
   v1 = v4 + 1 ;  \
   v2 = v1 + jp ; \
   v3 = v4 + jp ; \
   v5 = v1 + kp ; \
   v6 = v2 + kp ; \
   v7 = v3 + kp ; \
   v8 = v4 + kp ;

MarchingCubesAlgo1::MarchingCubesAlgo1(const conduit::Node &dom,
                                       const std::string &coordsetName,
                                       const std::string &fcnField,
                                       const std::string &maskField)
  : _dom(nullptr)
  , _ndim(0)
  , _logicalSize()
  , _logicalOrigin()
  , _coordsetPath("coordsets/" + coordsetName)
  , _fcnPath("fields/" + fcnField)
  , _maskPath(maskField.empty() ? std::string() : "fields/" + maskField)
  , _surfaceMesh(nullptr)
  , _cellIdField()
  , _contourVal(0.0)
{
  set_domain(dom);
  return;
}

void MarchingCubesAlgo1::set_domain(const conduit::Node &dom)
{
  SLIC_ASSERT_MSG(!conduit::blueprint::mesh::is_multi_domain(dom),
                  "MarchingCubesAlgo1 is single-domain only.  Try MarchingCubesAlgo.");

  SLIC_ASSERT(dom.has_path(_coordsetPath));
  SLIC_ASSERT(dom["topologies/mesh/type"].as_string() == "structured");

  SLIC_ASSERT(dom.has_path(_fcnPath));
  SLIC_ASSERT(dom[_fcnPath + "/association"].as_string() == "vertex");
  SLIC_ASSERT(dom.has_path(_fcnPath + "/values"));

  if(!_maskPath.empty())
  {
    SLIC_ASSERT(dom.has_path(_maskPath + "/values"));
  }

  _dom = &dom;

  const conduit::Node& dimsNode = _dom->fetch_existing("topologies/mesh/elements/dims");

  _ndim = dimsNode.number_of_children();

  _logicalSize.resize(_ndim);
  for(int d=0; d<_ndim; ++d)
  {
    _logicalSize[d] = dimsNode[d].as_int();
  }

  _logicalOrigin.resize(_ndim, 0);
  if (_dom->has_path("topologies/mesh/elements/origin"))
  {
    const conduit::Node& origins = _dom->fetch_existing("topologies/mesh/elements/origin");
    for(int d=0; d<_ndim; ++d)
    {
      _logicalOrigin[d] = origins[d].as_int();
    }
  }

  SLIC_ASSERT(_ndim >= 1 && _ndim <= 3);

  const conduit::Node& coordsValues = dom[_coordsetPath + "/values"];
  bool isInterleaved = conduit::blueprint::mcarray::is_interleaved(coordsValues);
  SLIC_ASSERT_MSG(!isInterleaved,
                  "MarchingCubesAlgo currently requires contiguous coordinates layout.");
}


/*!
  @brief Set the output surface mesh object.
*/
void MarchingCubesAlgo1::set_output_mesh(axom::mint::Mesh *surfaceMesh)
{
  _surfaceMesh = surfaceMesh;

  if(_ndim == 2)
  {
    using SegmentMesh = axom::mint::UnstructuredMesh< axom::mint::SINGLE_SHAPE >;
    SegmentMesh* mesh = dynamic_cast< SegmentMesh* >( _surfaceMesh );
    SLIC_ASSERT_MSG(mesh, "Surface mesh for 2D problem must be a SegmentMesh");
  }
  else
  {
    using TriangleMesh = axom::mint::UnstructuredMesh< axom::mint::SINGLE_SHAPE >;
    TriangleMesh* mesh = dynamic_cast< TriangleMesh* >( _surfaceMesh );
    SLIC_ASSERT_MSG(mesh, "Surface mesh for 3D problem must be a TriangleMesh");
  }
}


void MarchingCubesAlgo1::compute_iso_surface(double contourVal)
{
  SLIC_ASSERT_MSG(_surfaceMesh,
                  "You must call set_output_mesh before compute_iso_surface.");

  /*
    Notes: how to handle non-interleaved coordinates.
    We should be able to handle any blueprint-compatible domains.

    Need an array interface that supports strides.

    conduit DataArray uses DataType, which has strides.  But I'm not convinced
    that it actually uses those strides to compute offsets.  Need to do a practice
    problem or try it here.

    axom::ArrayView doesn't have strides.  It can support multicomponent arrays
    so I can represent an Nx3 or 3xN layout.  However, I would have to switch
    point index and dimension index depending on which layout.  Won't work.
  */

  const conduit::Node& coordValues = _dom->fetch_existing(_coordsetPath + "/values");
  const double* xPtr = coordValues["x"].as_double_ptr();
  const double* yPtr = _ndim >= 2 ? coordValues["y"].as_double_ptr() : nullptr;
  const double* zPtr = _ndim >= 3 ? coordValues["z"].as_double_ptr() : nullptr;

  auto& fcnValues = _dom->fetch_existing(_fcnPath + "/values");
  const double *fcnPtr = fcnValues.as_double_ptr();

  const int* maskPtr = nullptr;
  if(!_maskPath.empty())
  {
    auto& maskValues = _dom->fetch_existing(_maskPath + "/values");
    maskPtr = maskValues.as_int_ptr();
  }

  if(!_cellIdField.empty() && !_surfaceMesh->hasField(_cellIdField, axom::mint::CELL_CENTERED))
  {
    _surfaceMesh->createField<int>(_cellIdField, axom::mint::CELL_CENTERED);
  }

  _contourVal = contourVal;

  if ( _ndim==2 ) {
    /*
      For now, assume zero offsets and zero ghost width.
      Eventually, we'll have to support index offsets to
      handle data with ghosts.
    */
#if 1
    axom::StackArray<axom::IndexType, 2> cShape{_logicalSize[0], _logicalSize[1]};
    axom::StackArray<axom::IndexType, 2> nShape{1+_logicalSize[0], 1+_logicalSize[1]};
    axom::ArrayView<const double, 2> fcnView(fcnPtr, nShape);
    axom::ArrayView<const double, 2> xView(xPtr, nShape);
    axom::ArrayView<const double, 2> yView(yPtr, nShape);
    axom::ArrayView<const int, 2> maskView(maskPtr, cShape);
#else
    RAJA::OffsetLayout<2> fieldLayout = RAJA::make_offset_layout<2>(
      {_logicalOrigin[0], 1 + _logicalSize[0] + _logicalOrigin[0]},
      {_logicalOrigin[1], 1 + _logicalSize[1] + _logicalOrigin[1]} );
    RAJA::View<double, 2> fcnView(fcnPtr, fieldLayout);
    RAJA::View<double, 2> xView(xPtr, _logicalSize[0], _logicalSize[1]);
    RAJA::View<double, 2> yView(yPtr, _logicalSize[0], _logicalSize[1]);
#endif

#if 0
    // Write as RAJA::kernel<...>
    using KERNEL_EXEC_POL_SEQ =
      RAJA::KernelPolicy<
        RAJA::statement::For<1, RAJA::seq_exec,
                             RAJA::statement::For<0, RAJA::seq_exec,
                                                  RAJA::statement::Lambda<0>
                                                  >
                             >
      >;
    RAJA::kernel<KERNEL_EXEC_POL_SEQ>( RAJA::make_tuple(col_Range, row_Range),
                                       [=](int col, int row) {
                                         // computation for zone (row, col)
                                       } );
#elif 1
    // Write as regular nested loops.
    for(int j=0; j<_logicalSize[1]; ++j)
    {
      for(int i=0; i<_logicalSize[0]; ++i)
      {
        const bool skipZone = maskPtr && bool(maskView(i, j));
        if ( !skipZone ) {

           double vfs[4];
           double xx[4];
           double yy[4];

           vfs[0] = fcnView(i+1, j);
           vfs[1] = fcnView(i+1, j+1);
           vfs[2] = fcnView(i, j+1);
           vfs[3] = fcnView(i, j);

           xx[0] = xView(i+1, j);
           xx[1] = xView(i+1, j+1);
           xx[2] = xView(i, j+1);
           xx[3] = xView(i, j);

           yy[0] = yView(i+1, j);
           yy[1] = yView(i+1, j+1);
           yy[2] = yView(i, j+1);
           yy[3] = yView(i, j);

           auto nPrev = _surfaceMesh->getNumberOfCells();
           this->contourCell2D( xx, yy, vfs );
           auto nNew = _surfaceMesh->getNumberOfCells();

           if(nNew > nPrev && !_cellIdField.empty())
           {
             int zoneIdx = i + j * _logicalSize[0]; // TODO: Fix for ghost layer size.
             auto *cellIdPtr =
               _surfaceMesh->getFieldPtr<int>(_cellIdField, axom::mint::CELL_CENTERED);
             for(int n=nPrev; n<nNew; ++n)
             {
               cellIdPtr[n] = zoneIdx;
             }
           }
        } // END if
      }
    }
#endif

  } else {

#if 1
    axom::StackArray<axom::IndexType, 3> cShape{_logicalSize[0], _logicalSize[1], _logicalSize[2]};
    axom::StackArray<axom::IndexType, 3> nShape{1+_logicalSize[0], 1+_logicalSize[1], 1+_logicalSize[2]};
    axom::ArrayView<const double, 3> fcnView(fcnPtr, nShape);
    axom::ArrayView<const double, 3> xView(xPtr, nShape);
    axom::ArrayView<const double, 3> yView(yPtr, nShape);
    axom::ArrayView<const double, 3> zView(zPtr, nShape);
    axom::ArrayView<const int, 3> maskView(maskPtr, cShape);
    // Write as regular nested loops.
    for(int k=0; k<_logicalSize[2]; ++k)
    {
      for(int j=0; j<_logicalSize[1]; ++j)
      {
        for(int i=0; i<_logicalSize[0]; ++i)
        {
          const bool skipZone = maskPtr && bool(maskView(i, j, k));
          if ( !skipZone ) {

            double vfs[8];
            double xx[8];
            double yy[8];
            double zz[8];

            vfs[0] = fcnView(i+1, j  , k);
            vfs[1] = fcnView(i+1, j+1, k);
            vfs[2] = fcnView(i  , j+1, k);
            vfs[3] = fcnView(i  , j  , k);
            vfs[4] = fcnView(i+1, j  , k+1);
            vfs[5] = fcnView(i+1, j+1, k+1);
            vfs[6] = fcnView(i  , j+1, k+1);
            vfs[7] = fcnView(i  , j  , k+1);

            xx[0] = xView(i+1, j  , k);
            xx[1] = xView(i+1, j+1, k);
            xx[2] = xView(i  , j+1, k);
            xx[3] = xView(i  , j  , k);
            xx[4] = xView(i+1, j  , k+1);
            xx[5] = xView(i+1, j+1, k+1);
            xx[6] = xView(i  , j+1, k+1);
            xx[7] = xView(i  , j  , k+1);

            yy[0] = yView(i+1, j  , k);
            yy[1] = yView(i+1, j+1, k);
            yy[2] = yView(i  , j+1, k);
            yy[3] = yView(i  , j  , k);
            yy[4] = yView(i+1, j  , k+1);
            yy[5] = yView(i+1, j+1, k+1);
            yy[6] = yView(i  , j+1, k+1);
            yy[7] = yView(i  , j  , k+1);

            zz[0] = zView(i+1, j  , k);
            zz[1] = zView(i+1, j+1, k);
            zz[2] = zView(i  , j+1, k);
            zz[3] = zView(i  , j  , k);
            zz[4] = zView(i+1, j  , k+1);
            zz[5] = zView(i+1, j+1, k+1);
            zz[6] = zView(i  , j+1, k+1);
            zz[7] = zView(i  , j  , k+1);

            auto nPrev = _surfaceMesh->getNumberOfCells();
            this->contourCell3D( xx, yy, zz, vfs );
            auto nNew = _surfaceMesh->getNumberOfCells();


           if(nNew > nPrev && !_cellIdField.empty())
           {
             int zoneIdx = i + j * _logicalSize[0] + k * _logicalSize[0] * _logicalSize[1]; // TODO: Fix for ghost layer size.
             auto *cellIdPtr =
               _surfaceMesh->getFieldPtr<int>(_cellIdField, axom::mint::CELL_CENTERED);
             for(int n=nPrev; n<nNew; ++n)
             {
               cellIdPtr[n] = zoneIdx;
             }
           }
          } // END if
        }
      }
    }
#else
    double* zPtr = _ndim > 3 ? coordValues["z"].as_double_ptr() : nullptr;
    double* vf1, *vf2, *vf3, *vf4, *vf5, *vf6, *vf7, *vf8;
    double* x1,  *x2,  *x3,  *x4, *x5, *x6, *x7, *x8;
    double* y1,  *y2,  *y3,  *y4, *y5, *y6, *y7, *y8;
    double* z1,  *z2,  *z3,  *z4, *z5, *z6, *z7, *z8;
    int jp = _dom->fetch_existing(_coordsetPath + "/dims/i").as_int();
    int kp = _dom->fetch_existing(_coordsetPath + "/dims/j").as_int() * jp;

    NDSET3D( vf, vf1, vf2, vf3, vf4, vf5, vf6, vf7, vf8 );
    NDSET3D( x, x1, x2, x3, x4, x5, x6, x7, x8 );
    NDSET3D( y, y1, y2, y3, y4, y5, y6, y7, y8 );
    NDSET3D( z, z1, z2, z3, z4, z5, z6, z7, z8 );

    for_all_zones< policy::seq >( domain, [&](int zoneIdx) {

        const bool skipZone = finest && (finest[ zoneIdx ] != level);
        if ( !skipZone ) {

            double vfs[ 8 ];
            double xx[ 8 ];
            double yy[ 8 ];
            double zz[ 8 ];

            vfs[ 0 ] = vf1[ zoneIdx ];
            vfs[ 1 ] = vf2[ zoneIdx ];
            vfs[ 2 ] = vf3[ zoneIdx ];
            vfs[ 3 ] = vf4[ zoneIdx ];
            vfs[ 4 ] = vf5[ zoneIdx ];
            vfs[ 5 ] = vf6[ zoneIdx ];
            vfs[ 6 ] = vf7[ zoneIdx ];
            vfs[ 7 ] = vf8[ zoneIdx ];

            xx[ 0 ] = x1[ zoneIdx ];
            xx[ 1 ] = x2[ zoneIdx ];
            xx[ 2 ] = x3[ zoneIdx ];
            xx[ 3 ] = x4[ zoneIdx ];
            xx[ 4 ] = x5[ zoneIdx ];
            xx[ 5 ] = x6[ zoneIdx ];
            xx[ 6 ] = x7[ zoneIdx ];
            xx[ 7 ] = x8[ zoneIdx ];

            yy[ 0 ] = y1[ zoneIdx ];
            yy[ 1 ] = y2[ zoneIdx ];
            yy[ 2 ] = y3[ zoneIdx ];
            yy[ 3 ] = y4[ zoneIdx ];
            yy[ 4 ] = y5[ zoneIdx ];
            yy[ 5 ] = y6[ zoneIdx ];
            yy[ 6 ] = y7[ zoneIdx ];
            yy[ 7 ] = y8[ zoneIdx ];

            zz[ 0 ] = z1[ zoneIdx ];
            zz[ 1 ] = z2[ zoneIdx ];
            zz[ 2 ] = z3[ zoneIdx ];
            zz[ 3 ] = z4[ zoneIdx ];
            zz[ 4 ] = z5[ zoneIdx ];
            zz[ 5 ] = z6[ zoneIdx ];
            zz[ 6 ] = z7[ zoneIdx ];
            zz[ 7 ] = z8[ zoneIdx ];

            this->contourCell3D( xx, yy, zz, vfs );
        } // END if

    } );

#endif
  } // END else

}

//------------------------------------------------------------------------------
void MarchingCubesAlgo1::linear_interp( int edgeIdx,
                                       const double* xx,
                                       const double* yy,
                                       const double* zz,
                                       const double* f,
                                       double* xyz )
{
  SLIC_ASSERT( xx != NULL );
  SLIC_ASSERT( yy != NULL );
  SLIC_ASSERT( f  != NULL );

  // STEP 0: get the edge node indices
  // 2 nodes define the edge.  n1 and n2 are the indices of
  // the nodes w.r.t. the square or cubic zone.  There is a
  // agreed-on ordering of these indices in the arrays xx, yy,
  // zz, f, xyz.
  int n1 = edgeIdx;
  int n2 = (edgeIdx==3)? 0 : edgeIdx+1;

  if ( _ndim==3 ) {
    SLIC_ASSERT( zz != NULL );

    const int hex_edge_table[ ] = {
             0,1, 1,2, 2,3, 3,0, // base
             4,5, 5,6, 6,7, 7,4, // top
             0,4, 1,5, 2,6, 3,7  // vertical
     };

     n1 = hex_edge_table[ edgeIdx*2   ];
     n2 = hex_edge_table[ edgeIdx*2+1 ];

  } // END if 3-D

  // STEP 1: get the fields and coordinates from the two points
  const double f1 = f[ n1 ];
  const double f2 = f[ n2 ];

  double p1[ 3 ];
  p1[ 0 ] = xx[ n1 ];
  p1[ 1 ] = yy[ n1 ];
  p1[ 2 ] = 0.0;

  double p2[ 3 ];
  p2[ 0 ] = xx[ n2 ];
  p2[ 1 ] = yy[ n2 ];
  p2[ 2 ] = 0.0;

  if ( _ndim==3 ) {
    // set the z--coordinate if in 3-D
    p1[ 2 ] = zz[ n1 ];
    p2[ 2 ] = zz[ n2 ];
  }

  // STEP 2: check whether the interpolated point is at one of the two corners.
  if ( axom::utilities::isNearlyEqual( _contourVal, f1 ) ||
       axom::utilities::isNearlyEqual( f1,f2) ) {

      memcpy( xyz, p1, _ndim*sizeof(double) );
      return;
  }

  if ( axom::utilities::isNearlyEqual(_contourVal, f2 ) ) {

     memcpy( xyz, p2, _ndim*sizeof(double) );
     return;
  }

  // STEP 3: point is in between the edge points, interpolate its position
  constexpr double ptiny = 1.0e-80;
  const double df = f2 - f1 + ptiny; //add ptiny to avoid division by zero
  const double w  = (_contourVal-f1) / df;
  for ( int i=0; i < _ndim; ++i ) {
     xyz[ i ] = p1[ i ] + w*( p2[ i ] - p1[ i ] );
  }

}

//------------------------------------------------------------------------------
int MarchingCubesAlgo1::computeIndex( const double* f )
{
  const int numNodes = ( _ndim==3 ) ? 8 : 4;

  int index = 0;
  for ( int i=0; i < numNodes; ++i ) {

      if ( f[ i ] >= _contourVal ) {
         const int mask = (1 << i);
         index |= mask;
      }
  }
  return ( index );
}

//------------------------------------------------------------------------------
void MarchingCubesAlgo1::contourCell2D( double xx[4], double yy[4], double cellValues[4] )
{
  SLIC_ASSERT( xx != NULL );
  SLIC_ASSERT( yy != NULL );
  SLIC_ASSERT( cellValues != NULL );

  // compute index
  int index = MarchingCubesAlgo1::computeIndex( cellValues );
  SLIC_ASSERT( (index >= 0) && (index < 16) );

  // short-circuit
  if ( detail::num_segments[ index ] == 0 ) {
    return;
  } // END if


  // Generate line segments
  using SegmentMesh = axom::mint::UnstructuredMesh< axom::mint::SINGLE_SHAPE >;
  SegmentMesh* mesh = static_cast< SegmentMesh* >( _surfaceMesh );
  SLIC_ASSERT( mesh != NULL );
  SLIC_ASSERT( mesh->getCellType() == axom::mint::SEGMENT );
  SLIC_ASSERT( mesh->getDimension() == 2 );

  int idx = mesh->getNumberOfNodes();
  int cell[2];
  double p[2];

  const int nsegs = detail::num_segments[ index ];

  for ( int i=0; i < nsegs; ++i ) {

     const int e1 = detail::cases2D[ index ][ i*2   ];
     const int e2 = detail::cases2D[ index ][ i*2+1 ];

     MarchingCubesAlgo1::linear_interp( e1, xx, yy, NULL, cellValues, p );
     mesh->appendNode( p[0], p[1] );
     cell[ 0 ] = idx;
     ++idx;

     MarchingCubesAlgo1::linear_interp( e2, xx, yy, NULL, cellValues, p );
     mesh->appendNode( p[0], p[1] );
     cell[ 1 ] = idx;
     ++idx;

     mesh->appendCell( cell );

  } // END for all segments

}

//------------------------------------------------------------------------------
void MarchingCubesAlgo1::contourCell3D( double xx[8], double yy[8], double zz[8], double f[8] )
{
  SLIC_ASSERT( xx != NULL );
  SLIC_ASSERT( yy != NULL );
  SLIC_ASSERT( zz != NULL );
  SLIC_ASSERT( f != NULL );

  // compute index
  int index = MarchingCubesAlgo1::computeIndex( f );
  SLIC_ASSERT( (index >= 0) && (index < 256) );

  // short-circuit
  if ( detail::num_triangles[ index ] == 0 ) {
    return;
  } // END if

  // Generate triangles
  using TriangleMesh = axom::mint::UnstructuredMesh< axom::mint::SINGLE_SHAPE >;
  TriangleMesh* mesh = static_cast< TriangleMesh* >( _surfaceMesh );
  SLIC_ASSERT( mesh != NULL );
  SLIC_ASSERT( mesh->getCellType() == axom::mint::TRIANGLE );
  SLIC_ASSERT( mesh->getDimension() == 3 );

  const int numTriangles = detail::num_triangles[ index ];

  int idx = mesh->getNumberOfNodes();
  int cell[3];
  double p[3];

  for ( int i=0; i < numTriangles; ++i )
  {
    const int e1 = detail::cases3D[ index ][ i*3   ];
    const int e2 = detail::cases3D[ index ][ i*3+1 ];
    const int e3 = detail::cases3D[ index ][ i*3+2 ];

    MarchingCubesAlgo1::linear_interp( e1, xx, yy, zz, f, p );
    mesh->appendNode( p[0], p[1], p[2] );
    cell[ 0 ] = idx;
    ++idx;

    MarchingCubesAlgo1::linear_interp( e2, xx, yy, zz, f, p );
    mesh->appendNode( p[0], p[1], p[2] );
    cell[ 1 ] = idx;
    ++idx;

    MarchingCubesAlgo1::linear_interp( e3, xx, yy, zz, f, p );
    mesh->appendNode( p[0], p[1], p[2] );
    cell[ 2 ] = idx;
    ++idx;

    mesh->appendCell( cell );

  } // END for all triangles

}



#if 0
// Code from app

struct Attributes {
  std::vector< int > Zone;
  std::vector< int > DomainID;
  std::vector< int > Level;
};

//------------------------------------------------------------------------------
IsoContour::IsoContour( const std::string& field, double isoval ) :
        m_field( field ),
        m_isoval( isoval ),
        m_saveAttributes( true ),
        m_contour_attributes( new Attributes() )
{

}

//------------------------------------------------------------------------------
IsoContour::~IsoContour()
{
  delete m_contour_attributes;
  m_contour_attributes = NULL;
}

//------------------------------------------------------------------------------
const int* IsoContour::getZoneAttribute() const
{
  if ( m_saveAttributes ) {
     return ( &(m_contour_attributes->Zone)[0] );
  }
  return NULL;
}

//------------------------------------------------------------------------------
const int* IsoContour::getLevelAttribute() const
{
  if ( m_saveAttributes ) {
    return ( &(m_contour_attributes->Level)[0] );
  }
  return NULL;
}

//------------------------------------------------------------------------------
const int* IsoContour::getDomainIDAttribute() const
{
  if ( m_saveAttributes ) {
    return ( &(m_contour_attributes->DomainID)[0] );
  }
  return NULL;
}
#endif


}  // end namespace quest
}  // end namespace axom

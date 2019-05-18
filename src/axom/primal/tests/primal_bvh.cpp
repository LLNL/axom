/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#include "gtest/gtest.h"

#include "axom/core/Types.hpp"
#include "axom/core/memory_management.hpp"

// primal includes
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/spatial_acceleration/BVH.hpp"
#include "axom/primal/spatial_acceleration/UniformGrid.hpp"

using namespace axom;
//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------
namespace
{

//------------------------------------------------------------------------------
void generate_aabbs_and_centroids2d( const double* origin,
                                     const double* h,
                                     IndexType Ni,
                                     IndexType Nj,
                                     IndexType& ncells,
                                     double*& aabbs,
                                     IndexType*& cellIds,
                                     double*& xc,
                                     double*& yc )
{
  // pre-condition sanity checks
  EXPECT_TRUE( aabbs == nullptr );
  EXPECT_TRUE( cellIds == nullptr );
  EXPECT_TRUE( xc == nullptr );
  EXPECT_TRUE( yc == nullptr );

  // calculate some constants
  constexpr int ndims    = 2;
  constexpr int stride   = 2 * ndims;
  const IndexType NCi    = Ni - 1;
  const IndexType NCj    = Nj - 1;
  const IndexType jp     = NCi;
  ncells =  NCi * NCj;

  // allocate output arrays
  aabbs   = axom::allocate< double >( ncells * stride );
  cellIds = axom::allocate< IndexType >( ncells );
  xc      = axom::allocate< double >( ncells );
  yc      = axom::allocate< double >( ncells );

  // compute aabbs & cenrtoids for each cell
  for ( IndexType j=0 ; j < NCj ; ++j )
  {
    const IndexType jp_offset = j * jp;
    for ( IndexType i=0 ; i < NCi ; ++i )
    {
      // compute linear cell index
      const IndexType cellIdx = i + jp_offset;
      cellIds[ cellIdx ]      = cellIdx;

      // compute min/max bounds of aabb
      const IndexType offset = cellIdx * stride;
      const double xmin      = origin[ 0 ] + i     * h[ 0 ] ;
      const double xmax      = origin[ 0 ] + (i+1) * h[ 0 ] ;
      const double ymin      = origin[ 1 ] + j     * h[ 1 ] ;
      const double ymax      = origin[ 1 ] + (j+1) * h[ 1 ] ;

      aabbs[ offset     ]    =  xmin ;
      aabbs[ offset + 1 ]    =  xmax ;
      aabbs[ offset + 2 ]    =  ymin ;
      aabbs[ offset + 3 ]    =  ymax ;

      // compute centroid
      xc[ cellIdx ] = 0.5 * ( xmin + xmax ) ;
      yc[ cellIdx ] = 0.5 * ( ymin + ymax ) ;

    } // END for all i
  } // END for all j

  // post-condition sanity checks
  EXPECT_TRUE( aabbs != nullptr );
  EXPECT_TRUE( cellIds != nullptr );
  EXPECT_TRUE( xc != nullptr );
  EXPECT_TRUE( yc != nullptr );
}

//------------------------------------------------------------------------------
void generate_aabbs_and_centroids3d( const double* origin,
                                     const double* h,
                                     IndexType Ni,
                                     IndexType Nj,
                                     IndexType Nk,
                                     IndexType& ncells,
                                     double*& aabbs,
                                     IndexType*& cellIds,
                                     double*& xc,
                                     double*& yc,
                                     double*& zc )
{
  // pre-condition sanity checks
  EXPECT_TRUE( aabbs == nullptr );
  EXPECT_TRUE( cellIds == nullptr );
  EXPECT_TRUE( xc == nullptr );
  EXPECT_TRUE( yc == nullptr );
  EXPECT_TRUE( zc == nullptr );

  // calculate some constants
  constexpr int ndims  = 3;
  constexpr int stride = 2 * ndims;
  const IndexType NCi  = Ni - 1;
  const IndexType NCj  = Nj - 1;
  const IndexType NCk  = Nk - 1;
  const IndexType jp   = NCi;
  const IndexType kp   = NCi * NCj;
  ncells = NCi * NCj * NCk;

  // allocate output arrays
  aabbs   = axom::allocate< double >( ncells * stride );
  cellIds = axom::allocate< IndexType >( ncells );
  xc      = axom::allocate< double >( ncells );
  yc      = axom::allocate< double >( ncells );
  zc      = axom::allocate< double >( ncells );

  // compute aabbs & centroids for each cell
  for ( IndexType k=0 ; k < NCk ; ++k )
  {
    const IndexType kp_offset = k * kp ;
    for ( IndexType j=0 ; j < NCj ; ++j )
    {
      const IndexType jp_offset = j * jp ;
      for ( IndexType i=0 ; i < NCi ; ++i )
      {
        // compute linear cell index
        const IndexType cellIdx = i + jp_offset + kp_offset ;
        cellIds[ cellIdx ]      = cellIdx ;

        // compute min/max bounds of aabb
        const IndexType offset  = cellIdx * stride ;
        const double xmin       = origin[ 0 ] + i     * h[ 0 ] ;
        const double xmax       = origin[ 0 ] + (i+1) * h[ 0 ] ;
        const double ymin       = origin[ 1 ] + j     * h[ 1 ] ;
        const double ymax       = origin[ 1 ] + (j+1) * h[ 1 ] ;
        const double zmin       = origin[ 2 ] + k     * h[ 2 ] ;
        const double zmax       = origin[ 2 ] + (k+1) * h[ 2 ] ;

        aabbs[ offset     ] = xmin ;
        aabbs[ offset + 1 ] = xmax ;
        aabbs[ offset + 2 ] = ymin ;
        aabbs[ offset + 3 ] = ymax ;
        aabbs[ offset + 4 ] = zmin ;
        aabbs[ offset + 5 ] = zmax ;

        // compute centroid
        xc[ cellIdx ] = 0.5 * ( xmin + xmax );
        yc[ cellIdx ] = 0.5 * ( ymin + ymax );
        zc[ cellIdx ] = 0.5 * ( zmin + zmax );

      } // END for all i
    } // END for all j
  } // END for all k

  // post-condition sanity checks
  EXPECT_TRUE( aabbs != nullptr );
  EXPECT_TRUE( cellIds != nullptr );
  EXPECT_TRUE( xc != nullptr );
  EXPECT_TRUE( yc != nullptr );
  EXPECT_TRUE( zc != nullptr );
}

//------------------------------------------------------------------------------
template < typename FloatType >
void check_build_bvh2d( )
{
  constexpr int NUM_BOXES = 2;
  constexpr int NDIMS     = 2;
  FloatType boxes[ ]      = { 0., 0., 1., 1.,
                              1., 1., 2., 2. };

  primal::BVH< NDIMS, FloatType > bvh( boxes, NUM_BOXES );
  bvh.build( );

  FloatType lo[ NDIMS ];
  FloatType hi[ NDIMS ];
  bvh.getBounds( lo, hi );

  for ( int idim=0; idim < NDIMS; ++idim )
  {
    EXPECT_DOUBLE_EQ( lo[ idim ], 0.0 );
    EXPECT_DOUBLE_EQ( hi[ idim ], 2.0 );
  }

}

//------------------------------------------------------------------------------
template < typename FloatType >
void check_build_bvh3d( )
{
  constexpr int NUM_BOXES = 2;
  constexpr int NDIMS     = 3;
  FloatType boxes[ ]      = { 0., 0., 0., 1., 1., 1.,
                              1., 1., 1., 2., 2., 2. };

  primal::BVH< NDIMS, FloatType > bvh( boxes, NUM_BOXES );
  bvh.build( );

  FloatType lo[ NDIMS ];
  FloatType hi[ NDIMS ];
  bvh.getBounds( lo, hi );

  for ( int idim=0; idim < NDIMS; ++idim )
  {
    EXPECT_DOUBLE_EQ( lo[ idim ], 0.0 );
    EXPECT_DOUBLE_EQ( hi[ idim ], 2.0 );
  }

}

} /* end unnamed namespace */

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
TEST( primal_bvh, contruct2D)
{
  check_build_bvh2d< float >( );
  check_build_bvh2d< double >( );
}

//------------------------------------------------------------------------------
TEST( primal_bvh, contruct3D)
{
  check_build_bvh3d< float >( );
  check_build_bvh3d< double >( );
}


//------------------------------------------------------------------------------
//TEST( primal_bvh, check_find_3d )
//{
//  using PointType         = primal::Point< double, 3 >;
//  constexpr int NDIMS     = 3;
//  constexpr IndexType N   = 4;
//
//  double lo[ NDIMS ] = { 0.0, 0.0, 0.0 };
//  double hi[ NDIMS ] = { 3.0, 3.0, 3.0 };
//  int res[ NDIMS ]   = { N-1, N-1, N-1 };
//
//  double mesh_origin[ NDIMS ]  = { 0.0, 0.0, 0.0 };
//  double mesh_spacing[ NDIMS ] = { 1.0, 1.0, 1.0 };
//
//  double* aabbs      = nullptr;
//  IndexType* cellIds = nullptr;
//  double* xc         = nullptr;
//  double* yc         = nullptr;
//  double* zc         = nullptr;
//  IndexType ncells   = 0;
//
//  // generate bounding boxes
//  generate_aabbs_and_centroids3d(
//      mesh_origin, mesh_spacing, N, N, N, ncells, aabbs, cellIds, xc, yc, zc );
//  EXPECT_EQ( ncells, (N-1) * (N-1) * (N-1) );
//
//  // construct the BVH
//  primal::BVH bvh( NDIMS, aabbs, ncells );
//  bvh.build( );
//
//  // traverse the BVH to find the candidates for all the centroids
//  IndexType* offsets    = axom::allocate< IndexType >( ncells );
//  IndexType* candidates = nullptr;
//  bvh.find( offsets, candidates, ncells, xc, yc, zc );
//
//  EXPECT_TRUE( candidates != nullptr );
//
//  primal::UniformGrid< IndexType, NDIMS > ug( lo, hi, res );
//
//  for ( IndexType i=0; i < ncells; ++i )
//  {
//    PointType q            = PointType::make_point( xc[ i ],yc[ i ],zc[ i ] );
//    const int donorCellIdx = ug.getBinIndex( q );
//    EXPECT_EQ( donorCellIdx, cellIds[ i ] );
//
//  } // END for all cell centroids
//
//  axom::deallocate( offsets );
//// TODO: add line below
////  axom::free( candidates );
//  axom::deallocate( aabbs );
//  axom::deallocate( cellIds );
//  axom::deallocate( xc );
//  axom::deallocate( yc );
//  axom::deallocate( zc );
//}

//------------------------------------------------------------------------------
//TEST( primal_bvh, check_find_2d )
//{
//  using PointType         = primal::Point< double, 2 >;
//  constexpr int NDIMS     = 2;
//  constexpr IndexType N   = 4;
//
//  double lo[ NDIMS ] = { 0.0, 0.0 };
//  double hi[ NDIMS ] = { 3.0, 3.0 };
//  int res[ NDIMS ]   = { N-1, N-1 };
//
//  double mesh_origin[ NDIMS ]  = { 0.0, 0.0 };
//  double mesh_spacing[ NDIMS ] = { 1.0, 1.0 };
//
//  double* aabbs      = nullptr;
//  IndexType* cellIds = nullptr;
//  double* xc         = nullptr;
//  double* yc         = nullptr;
//  IndexType ncells   = 0;
//
//  // generate bounding boxes
//  generate_aabbs_and_centroids2d(
//      mesh_origin, mesh_spacing, N, N, ncells, aabbs, cellIds, xc, yc );
//  EXPECT_EQ( ncells, (N-1) * (N-1) );
//
//  // construct the BVH
//  primal::BVH bvh( NDIMS, aabbs, ncells );
//  bvh.build( );
//
//  // traverse the BVH to find the candidates for all the centroids
//  IndexType* offsets    = axom::allocate< IndexType >( ncells );
//  IndexType* candidates = nullptr;
//  bvh.find( offsets, candidates, ncells, xc, yc );
//// TODO: add check below
////  EXPECT_TRUE( candidates != nullptr );
//
//  primal::UniformGrid< IndexType, NDIMS > ug( lo, hi, res );
//
//  for ( IndexType i=0 ; i < ncells ; ++i )
//  {
//    PointType q            = PointType::make_point( xc[ i ], yc[ i ] );
//    const int donorCellIdx = ug.getBinIndex( q );
//    EXPECT_EQ( donorCellIdx, cellIds[ i ] );
//
//    // TODO: check candidates found by BVH
//
//  } // END for all cell centroids
//
//  axom::deallocate( offsets );
//// TODO: add line below
////  axom::free( candidates );
//  axom::deallocate( aabbs );
//  axom::deallocate( cellIds );
//  axom::deallocate( xc );
//  axom::deallocate( yc );
//}

//------------------------------------------------------------------------------
#include "axom/slic/core/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}

// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// axom/core includes
#include "axom/config.hpp"

#include "axom/core/Types.hpp"
#include "axom/core/execution/execution_space.hpp"
#include "axom/core/memory_management.hpp"
#include "axom/core/numerics/Matrix.hpp"

// axom/primal includes
#include "axom/primal/geometry/Point.hpp"

// axom/spin includes
#include "axom/spin/BVH.hpp"
#include "axom/spin/UniformGrid.hpp"

#include "axom/spin/internal/linear_bvh/QueryAccessor.hpp"
#include "axom/spin/internal/linear_bvh/TraversalPredicates.hpp"

// axom/mint includes
#include "axom/mint/mesh/Mesh.hpp"
#include "axom/mint/mesh/UniformMesh.hpp"
#include "axom/mint/execution/interface.hpp"
#include "axom/mint/utils/vtk_utils.hpp"

// gtest includes
#include "gtest/gtest.h"

using namespace axom;
namespace xargs  = mint::xargs;

//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------
namespace
{

//------------------------------------------------------------------------------

/*!
 * \brief Give a 2D mesh object, this method generates an array of axis-aligned
 *  bounding boxes corresponding to each constituent cell of the mesh and
 *  a corresponding centroid.
 *
 * \param [in]  mesh pointer to the mesh object.
 * \param [out] aabbs flat array of bounding boxes [xmin,ymin,xmax,ymax....]
 * \param [out] xc buffer to store the x-component of the cell centroids
 * \param [out] yc buffer to store the y-component of the cell centroids
 *
 * \note The intent of this method is to generate synthetic input test data
 *  to test the functionality of the BVH.
 *
 * \warning This method allocates aabbs internally. The caller is responsible
 *  for properly deallocating aabbs.
 *
 * \pre aabbs == nullptr
 * \pre xc != nullptr
 * \pre yc != nullptr
 *
 * \post aabbs != nullptr
 */
template < typename FloatType >
void generate_aabbs_and_centroids2d( const mint::Mesh* mesh,
                                     FloatType*& aabbs,
                                     FloatType* xc,
                                     FloatType* yc )
{
  // sanity checks
  EXPECT_TRUE( aabbs == nullptr );
  EXPECT_TRUE( xc != nullptr );
  EXPECT_TRUE( yc != nullptr );

  // calculate some constants
  constexpr double ONE_OVER_4 = 1.f / 4.f;
  constexpr int ndims         = 2;
  constexpr int stride        = 2 * ndims;

  EXPECT_EQ( ndims, mesh->getDimension() );

  // allocate output arrays
  const IndexType ncells = mesh->getNumberOfCells();
  aabbs = axom::allocate< FloatType >( ncells * stride );

  using exec_policy = axom::SEQ_EXEC;
  mint::for_all_cells< exec_policy, xargs::coords >(
      mesh, AXOM_LAMBDA( IndexType cellIdx,
                         numerics::Matrix< double >& coords,
                         const IndexType* AXOM_NOT_USED(nodeIds) )
  {

    spin::internal::linear_bvh::Range< double > xrange;
    spin::internal::linear_bvh::Range< double > yrange;

    double xsum = 0.0;
    double ysum = 0.0;

    for ( IndexType inode=0; inode < 4; ++inode )
    {
      const double* node = coords.getColumn( inode );
      xsum += node[ mint::X_COORDINATE ];
      ysum += node[ mint::Y_COORDINATE ];

      xrange.include( node[ mint::X_COORDINATE ] );
      yrange.include( node[ mint::Y_COORDINATE ] );
    } // END for all cells nodes

    xc[ cellIdx ] = xsum * ONE_OVER_4;
    yc[ cellIdx ] = ysum * ONE_OVER_4;
    const IndexType offset  = cellIdx * stride ;
    aabbs[ offset     ] = xrange.min();
    aabbs[ offset + 1 ] = yrange.min();
    aabbs[ offset + 2 ] = xrange.max();
    aabbs[ offset + 3 ] = yrange.max();

  } );

  // post-condition sanity checks
  EXPECT_TRUE( aabbs != nullptr );

}

//------------------------------------------------------------------------------

/*!
 * \brief Give a 3D mesh object, this method generates an array of axis-aligned
 *  bounding boxes (AABBS) corresponding to each constituent cell of the mesh
 *  and a corresponding centroid.
 *
 * \param [in]  mesh pointer to the mesh object.
 * \param [out] aabbs flat array of AABBS [xmin,ymin,zmin,xmax,ymax,zmax....]
 * \param [out] xc buffer to store the x-component of the cell centroids
 * \param [out] yc buffer to store the y-component of the cell centroids
 * \param [out] zc buffer to store the z-component of the cell centroids
 *
 * \note The intent of this method is to generate synthetic input test data
 *  to test the functionality of the BVH.
 *
 * \warning This method allocates aabbs internally. The caller is responsible
 *  for properly deallocating aabbs.
 *
 * \pre aabbs == nullptr
 * \pre xc != nullptr
 * \pre yc != nullptr
 * \pre zc != nullptr
 *
 * \post aabbs != nullptr
 */
template < typename FloatType >
void generate_aabbs_and_centroids3d( const mint::Mesh* mesh,
                                     FloatType*& aabbs,
                                     FloatType* xc,
                                     FloatType* yc,
                                     FloatType* zc )
{
  // sanity checks
  EXPECT_TRUE( aabbs == nullptr );
  EXPECT_TRUE( xc != nullptr );
  EXPECT_TRUE( yc != nullptr );
  EXPECT_TRUE( zc != nullptr );

  // calculate some constants
  constexpr double ONE_OVER_8 = 1.f / 8.f;
  constexpr int ndims         = 3;
  constexpr int stride        = 2 * ndims;

  EXPECT_EQ( ndims, mesh->getDimension() );

  // allocate output arrays
  const IndexType ncells = mesh->getNumberOfCells();
  aabbs = axom::allocate< FloatType >( ncells * stride );

  using exec_policy = axom::SEQ_EXEC;
  mint::for_all_cells< exec_policy, xargs::coords >(
      mesh, AXOM_LAMBDA( IndexType cellIdx,
                         numerics::Matrix< double >& coords,
                         const IndexType* AXOM_NOT_USED(nodeIds) )
  {

    spin::internal::linear_bvh::Range< double > xrange;
    spin::internal::linear_bvh::Range< double > yrange;
    spin::internal::linear_bvh::Range< double > zrange;

    double xsum = 0.0;
    double ysum = 0.0;
    double zsum = 0.0;

    for ( IndexType inode=0; inode < 8; ++inode )
    {
      const double* node = coords.getColumn( inode );
      xsum += node[ mint::X_COORDINATE ];
      ysum += node[ mint::Y_COORDINATE ];
      zsum += node[ mint::Z_COORDINATE ];

      xrange.include( node[ mint::X_COORDINATE ] );
      yrange.include( node[ mint::Y_COORDINATE ] );
      zrange.include( node[ mint::Z_COORDINATE ] );
    } // END for all cells nodes

    xc[ cellIdx ] = xsum * ONE_OVER_8;
    yc[ cellIdx ] = ysum * ONE_OVER_8;
    zc[ cellIdx ] = zsum * ONE_OVER_8;

    const IndexType offset  = cellIdx * stride ;
    aabbs[ offset     ] = xrange.min();
    aabbs[ offset + 1 ] = yrange.min();
    aabbs[ offset + 2 ] = zrange.min();
    aabbs[ offset + 3 ] = xrange.max();
    aabbs[ offset + 4 ] = yrange.max();
    aabbs[ offset + 5 ] = zrange.max();

  } );

  // post-condition sanity checks
  EXPECT_TRUE( aabbs != nullptr );
}

//------------------------------------------------------------------------------

/*!
 * \brief Tests the construction of the BVH in 2D by inserting two bounding
 *  boxes in the BVH and ensuring that the bounds of the BVH are as expected.
 */
template < typename ExecSpace, typename FloatType >
void check_build_bvh2d( )
{
  constexpr int NUM_BOXES = 2;
  constexpr int NDIMS     = 2;

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator( axom::execution_space<ExecSpace>::allocatorID());

  FloatType* boxes = axom::allocate< FloatType >( 8 );
  boxes[ 0 ] = boxes[ 1 ] = 0.;
  boxes[ 2 ] = boxes[ 3 ] = 1.;
  boxes[ 4 ] = boxes[ 5 ] = 1.;
  boxes[ 6 ] = boxes[ 7 ] = 2.;

  spin::BVH< NDIMS, ExecSpace, FloatType > bvh( boxes, NUM_BOXES );
  bvh.setScaleFactor( 1.0 ); // i.e., no scaling
  bvh.build( );

  FloatType lo[ NDIMS ];
  FloatType hi[ NDIMS ];
  bvh.getBounds( lo, hi );

  for ( int idim=0 ; idim < NDIMS ; ++idim )
  {
    EXPECT_DOUBLE_EQ( lo[ idim ], 0.0 );
    EXPECT_DOUBLE_EQ( hi[ idim ], 2.0 );
  }

  axom::deallocate( boxes );
  axom::setDefaultAllocator( current_allocator );
}

//------------------------------------------------------------------------------

/*!
 * \brief Tests the construction of the BVH in 3D by inserting two bounding
 *  boxes in the BVH and ensuring that the bounds of the BVH are as expected.
 */
template < typename ExecSpace, typename FloatType >
void check_build_bvh3d( )
{
  constexpr int NUM_BOXES = 2;
  constexpr int NDIMS     = 3;

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator( axom::execution_space<ExecSpace>::allocatorID());

  FloatType* boxes = axom::allocate< FloatType >( 12 );
  boxes[ 0 ] = boxes[  1 ] = boxes[  2 ] = 0.;
  boxes[ 3 ] = boxes[  4 ] = boxes[  5 ] = 1.;
  boxes[ 6 ] = boxes[  7 ] = boxes[  8 ] = 1.;
  boxes[ 9 ] = boxes[ 10 ] = boxes[ 11 ] = 2.;

  spin::BVH< NDIMS, ExecSpace, FloatType > bvh( boxes, NUM_BOXES );
  bvh.setScaleFactor( 1.0 ); // i.e., no scaling
  bvh.build( );

  FloatType lo[ NDIMS ];
  FloatType hi[ NDIMS ];
  bvh.getBounds( lo, hi );

  for ( int idim=0 ; idim < NDIMS ; ++idim )
  {
    EXPECT_DOUBLE_EQ( lo[ idim ], 0.0 );
    EXPECT_DOUBLE_EQ( hi[ idim ], 2.0 );
  }

  axom::deallocate( boxes );
  axom::setDefaultAllocator( current_allocator );
}

//------------------------------------------------------------------------------

/*!
 * \brief Tests the find algorithm of the BVH in 3D.
 *
 *  A uniform mesh is used for this test where the query points are generated
 *  by taking the centroids of the constituent cells of the mesh. Since the
 *  mesh is a uniform, cartesian mesh the find algorithm should return exactly
 *  one candidate for each query point corresponding to the cell on the mesh
 *  that generated the centroid. This property is checked by using the
 *  spin::UniformGrid class.
 *
 *  In addition, the test shifts the points by an offset to ensure that points
 *  outside the mesh return no candidate.
 *
 */
template < typename ExecSpace, typename FloatType >
void check_find3d( )
{
  constexpr int NDIMS   = 3;
  constexpr IndexType N = 4;

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator( axom::execution_space<ExecSpace>::allocatorID());

  using PointType = primal::Point< double, NDIMS >;

  double lo[ NDIMS ] = { 0.0, 0.0, 0.0 };
  double hi[ NDIMS ] = { 3.0, 3.0, 3.0 };
  int res[ NDIMS ]   = { N-1, N-1, N-1 };

  mint::UniformMesh mesh( lo, hi, N, N, N );
  FloatType* xc = mesh.createField< FloatType >( "xc", mint::CELL_CENTERED );
  FloatType* yc = mesh.createField< FloatType >( "yc", mint::CELL_CENTERED );
  FloatType* zc = mesh.createField< FloatType >( "zc", mint::CELL_CENTERED );
  const IndexType ncells = mesh.getNumberOfCells();

  FloatType* aabbs = nullptr;
  generate_aabbs_and_centroids3d( &mesh, aabbs, xc, yc, zc );

  // construct the BVH
  spin::BVH< NDIMS, ExecSpace, FloatType > bvh( aabbs, ncells );
  bvh.setScaleFactor( 1.0 ); // i.e., no scaling
  bvh.build( );

  FloatType min[ NDIMS ];
  FloatType max[ NDIMS ];
  bvh.getBounds( min, max );
  for ( int i=0 ; i < NDIMS ; ++i )
  {
    EXPECT_DOUBLE_EQ( min[ i ], lo[ i ] );
    EXPECT_DOUBLE_EQ( max[ i ], hi[ i ] );
  }

  // traverse the BVH to find the candidates for all the centroids
  IndexType* offsets    = axom::allocate< IndexType >( ncells );
  IndexType* counts     = axom::allocate< IndexType >( ncells );
  IndexType* candidates = nullptr;
  bvh.findPoints( offsets, counts, candidates, ncells, xc, yc, zc );

  EXPECT_TRUE( candidates != nullptr );

  spin::UniformGrid< IndexType, NDIMS > ug( lo, hi, res );

  for ( IndexType i=0 ; i < ncells ; ++i )
  {
    PointType q            = PointType::make_point( xc[ i ],yc[ i ],zc[ i ] );
    const int donorCellIdx = ug.getBinIndex( q );
    EXPECT_EQ( counts[ i ], 1 );
    EXPECT_EQ( donorCellIdx, candidates[ offsets[ i ] ] );
  } // END for all cell centroids

  axom::deallocate( candidates );

  // check points that are outside by shifting the query points
  constexpr double OFFSET = 10.0;
  for ( IndexType i=0 ; i < ncells ; ++i )
  {
    xc[ i ] += OFFSET;
    yc[ i ] += OFFSET;
    zc[ i ] += OFFSET;
  }

  bvh.findPoints( offsets, counts, candidates, ncells, xc, yc, zc );

  for ( IndexType i=0 ; i < ncells ; ++i )
  {
    EXPECT_EQ( counts[ i ], 0 );
  }

  axom::deallocate( offsets );
  axom::deallocate( candidates );
  axom::deallocate( counts );
  axom::deallocate( aabbs );

  axom::setDefaultAllocator( current_allocator );
}

//------------------------------------------------------------------------------

/*!
 * \brief Tests the find algorithm of the BVH in 2D.
 *
 *  A uniform mesh is used for this test where the query points are generated
 *  by taking the centroids of the constituent cells of the mesh. Since the
 *  mesh is a uniform, cartesian mesh the find algorithm should return exactly
 *  one candidate for each query point corresponding to the cell on the mesh
 *  that generated the centroid. This property is checked by using the
 *  spin::UniformGrid class.
 *
 *  In addition, the test shifts the points by an offset to ensure that points
 *  outside the mesh return no candidate.
 *
 */
template < typename ExecSpace, typename FloatType >
void check_find2d( )
{
  constexpr int NDIMS   = 2;
  constexpr IndexType N = 4;

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator( axom::execution_space<ExecSpace>::allocatorID());

  using PointType = primal::Point< double, NDIMS >;

  double lo[ NDIMS ] = { 0.0, 0.0 };
  double hi[ NDIMS ] = { 3.0, 3.0 };
  int res[ NDIMS ]   = { N-1, N-1 };

  mint::UniformMesh mesh( lo, hi, N, N );
  FloatType* xc = mesh.createField< FloatType >( "xc", mint::CELL_CENTERED );
  FloatType* yc = mesh.createField< FloatType >( "yc", mint::CELL_CENTERED );
  const IndexType ncells = mesh.getNumberOfCells();

  FloatType* aabbs = nullptr;
  generate_aabbs_and_centroids2d( &mesh, aabbs, xc, yc );

  // construct the BVH
  spin::BVH< NDIMS, ExecSpace, FloatType > bvh( aabbs, ncells );
  bvh.setScaleFactor( 1.0 ); // i.e., no scaling
  bvh.build( );

  FloatType min[ NDIMS ];
  FloatType max[ NDIMS ];
  bvh.getBounds( min, max );
  for ( int i=0 ; i < NDIMS ; ++i )
  {
    EXPECT_DOUBLE_EQ( min[ i ], lo[ i ] );
    EXPECT_DOUBLE_EQ( max[ i ], hi[ i ] );
  }

  // traverse the BVH to find the candidates for all the centroids
  IndexType* offsets    = axom::allocate< IndexType >( ncells );
  IndexType* counts     = axom::allocate< IndexType >( ncells );
  IndexType* candidates = nullptr;
  bvh.findPoints( offsets, counts, candidates, ncells, xc, yc );

  EXPECT_TRUE( candidates != nullptr );

  spin::UniformGrid< IndexType, NDIMS > ug( lo, hi, res );

  for ( IndexType i=0 ; i < ncells ; ++i )
  {
    PointType q            = PointType::make_point( xc[ i ],yc[ i ] );
    const int donorCellIdx = ug.getBinIndex( q );
    EXPECT_EQ( counts[ i ], 1 );
    EXPECT_EQ( donorCellIdx, candidates[ offsets[ i ] ] );
  } // END for all cell centroids

  axom::deallocate( candidates );

  // check points that are outside by shifting the query points
  constexpr double OFFSET = 10.0;
  for ( IndexType i=0 ; i < ncells ; ++i )
  {
    xc[ i ] += OFFSET;
    yc[ i ] += OFFSET;
  }

  bvh.findPoints( offsets, counts, candidates, ncells, xc, yc );

  for ( IndexType i=0 ; i < ncells ; ++i )
  {
    EXPECT_EQ( counts[ i ], 0 );
  }

  axom::deallocate( offsets );
  axom::deallocate( candidates );
  axom::deallocate( counts );
  axom::deallocate( aabbs );

  axom::setDefaultAllocator( current_allocator );
}

//------------------------------------------------------------------------------

/*!
 * \brief Checks that the BVH behaves properly when user supplies a single box.
 *
 *  The Test inserts a single bounding box that spans [0,1] x [0,1] and tests
 *  that the centroid xc=(0.5,0.5) can be found. Moreover, it shifts the point
 *  by an offset and ensures that a point that is outside is not found.
 */
template < typename ExecSpace, typename FloatType >
void check_single_box2d( )
{
  constexpr int NUM_BOXES = 1;
  constexpr int NDIMS     = 2;

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator( axom::execution_space<ExecSpace>::allocatorID());

  // single bounding box in [0,1] x [0,1]
  FloatType* boxes = axom::allocate< FloatType >( 4 );
  boxes[ 0 ] = boxes[ 1 ] = 0.;
  boxes[ 2 ] = boxes[ 3 ] = 1.;

  // construct a BVH with a single box
  spin::BVH< NDIMS, ExecSpace, FloatType > bvh( boxes, NUM_BOXES );
  bvh.setScaleFactor( 1.0 ); // i.e., no scaling
  bvh.build( );

  // check the bounds -- should match the bounds of the input bounding box
  FloatType lo[ NDIMS ];
  FloatType hi[ NDIMS ];
  bvh.getBounds( lo, hi );

  for ( int idim=0 ; idim < NDIMS ; ++idim )
  {
    EXPECT_DOUBLE_EQ( lo[ idim ], 0.0 );
    EXPECT_DOUBLE_EQ( hi[ idim ], 1.0 );
  }

  // run the find algorithm w/ the centroid of the bounding box as input.
  // Should return one and only one candidate that corresponds to the
  // single bounding box.
  FloatType* xc = axom::allocate< FloatType >( NUM_BOXES );
  FloatType* yc = axom::allocate< FloatType >( NUM_BOXES );
  xc[ 0 ] = yc[ 0 ] = 0.5;

  IndexType* offsets    = axom::allocate< IndexType >( NUM_BOXES );
  IndexType* counts     = axom::allocate< IndexType >( NUM_BOXES );
  IndexType* candidates = nullptr;
  bvh.findPoints( offsets, counts, candidates, NUM_BOXES, xc, yc );
  EXPECT_TRUE( candidates != nullptr );
  EXPECT_EQ( counts[ 0 ], 1 );
  EXPECT_EQ( 0, candidates[ offsets[ 0 ] ] );
  axom::deallocate( candidates );

  // shift centroid outside of the BVH, should return no candidates.
  xc[ 0 ] += 10.0; yc[ 0 ] += 10.0;
  bvh.findPoints( offsets, counts, candidates, NUM_BOXES, xc, yc );
  EXPECT_EQ( counts[ 0 ], 0 );

  axom::deallocate( xc );
  axom::deallocate( yc );
  axom::deallocate( boxes );
  axom::deallocate( offsets );
  axom::deallocate( counts );
  axom::deallocate( candidates );
  axom::setDefaultAllocator( current_allocator );
}

//------------------------------------------------------------------------------
/*!
 * \brief Checks that the BVH behaves properly when user supplies a single box.
 *
 *  The Test inserts a single bounding box that spans [0,1] x [0,1] x [0,1] and
 *  tests that the centroid xc=(0.5,0.5) can be found. Moreover, it shifts the
 *  point by an offset and ensures that a point that is outside is not found.
 */
template < typename ExecSpace, typename FloatType >
void check_single_box3d( )
{
  constexpr int NUM_BOXES = 1;
  constexpr int NDIMS     = 3;

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator( axom::execution_space<ExecSpace>::allocatorID());

  // single bounding box in [0,1] x [0,1] x [0,1]
  FloatType* boxes = axom::allocate< FloatType >( 6 );
  boxes[ 0 ] = boxes[  1 ] = boxes[  2 ] = 0.;
  boxes[ 3 ] = boxes[  4 ] = boxes[  5 ] = 1.;

  // construct a BVH with a single box
  spin::BVH< NDIMS, ExecSpace, FloatType > bvh( boxes, NUM_BOXES );
  bvh.setScaleFactor( 1.0 ); // i.e., no scaling
  bvh.build( );

  // check the bounds -- should match the bounds of the input bounding box
  FloatType lo[ NDIMS ];
  FloatType hi[ NDIMS ];
  bvh.getBounds( lo, hi );

  for ( int idim=0 ; idim < NDIMS ; ++idim )
  {
    EXPECT_DOUBLE_EQ( lo[ idim ], 0.0 );
    EXPECT_DOUBLE_EQ( hi[ idim ], 1.0 );
  }

  // run the find algorithm w/ the centroid of the bounding box as input.
  // Should return one and only one candidate that corresponds to the
  // single bounding box.
  FloatType* xc = axom::allocate< FloatType >( NUM_BOXES );
  FloatType* yc = axom::allocate< FloatType >( NUM_BOXES );
  FloatType* zc = axom::allocate< FloatType >( NUM_BOXES );
  xc[ 0 ] = yc[ 0 ] = zc[ 0 ] = 0.5;

  IndexType* offsets    = axom::allocate< IndexType >( NUM_BOXES );
  IndexType* counts     = axom::allocate< IndexType >( NUM_BOXES );
  IndexType* candidates = nullptr;
  bvh.findPoints( offsets, counts, candidates, NUM_BOXES, xc, yc, zc );
  EXPECT_TRUE( candidates != nullptr );
  EXPECT_EQ( counts[ 0 ], 1 );
  EXPECT_EQ( 0, candidates[ offsets[ 0 ] ] );
  axom::deallocate( candidates );

  // shift centroid outside of the BVH, should return no candidates.
  xc[ 0 ] += 10.0; yc[ 0 ] += 10.0;
  bvh.findPoints( offsets, counts, candidates, NUM_BOXES, xc, yc, zc );
  EXPECT_EQ( counts[ 0 ], 0 );

  axom::deallocate( xc );
  axom::deallocate( yc );
  axom::deallocate( zc );
  axom::deallocate( boxes );
  axom::deallocate( offsets );
  axom::deallocate( counts );
  axom::deallocate( candidates );
  axom::setDefaultAllocator( current_allocator );
}

} /* end unnamed namespace */

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
TEST( spin_bvh, traversal_predicates_rayIntersectsLeftBin )
{
  namespace bvh               = axom::spin::internal::linear_bvh;
  using TraversalPredicates2D = bvh::TraversalPredicates< 2, double >;
  using TraversalPredicates3D = bvh::TraversalPredicates< 3, double >;
  using VecType               = bvh::Vec< double, 4 >;
  using RayType               = bvh::Vec< double, 6 >;

  VecType s1, s2;
  s1[ 0 ] = 0.; // LeftBin.xmin
  s1[ 1 ] = 0.; // LeftBin.ymin
  s1[ 2 ] = 0.; // LeftBin.zmin

  s1[ 3 ] = 1.; // LeftBin.xmax
  s2[ 0 ] = 1.; // LeftBin.ymax
  s2[ 1 ] = 1.; // LeftBin.zmax

  RayType ray2d;
  ray2d[ 0 ] = -1.0;
  ray2d[ 1 ] = -1.0;
  ray2d[ 2 ] =  1.0;
  ray2d[ 3 ] =  1.0;
  EXPECT_TRUE( TraversalPredicates2D::rayIntersectsLeftBin( ray2d, s1, s2 ) );

  // flip normal
  ray2d[ 2 ] = -1.0;
  ray2d[ 3 ] = -1.0;
  EXPECT_FALSE( TraversalPredicates2D::rayIntersectsLeftBin( ray2d, s1, s2 ) );

  RayType ray3d;
  ray3d[ 0 ] = -1.0;
  ray3d[ 1 ] = -1.0;
  ray3d[ 2 ] = -1.0;
  ray3d[ 3 ] = 1.0;
  ray3d[ 4 ] = 1.0;
  ray3d[ 5 ] = 1.0;
  EXPECT_TRUE( TraversalPredicates3D::rayIntersectsLeftBin( ray3d, s1, s2) );

  // flip normal
  ray3d[ 3 ] = -1.0;
  ray3d[ 4 ] = -1.0;
  ray3d[ 5 ] = -1.0;
  EXPECT_FALSE( TraversalPredicates3D::rayIntersectsLeftBin( ray3d, s1, s2 ) );
}

//------------------------------------------------------------------------------
TEST( spin_bvh, traversal_predicates_rayIntersectsRightBin )
{
  namespace bvh               = axom::spin::internal::linear_bvh;
  using TraversalPredicates2D = bvh::TraversalPredicates< 2, double >;
  using TraversalPredicates3D = bvh::TraversalPredicates< 3, double >;
  using VecType               = bvh::Vec< double, 4 >;
  using RayType               = bvh::Vec< double, 6 >;

  VecType s2, s3;
  s2[ 2 ] = 0.; // RightBin.xmin
  s2[ 3 ] = 0.; // RightBin.ymin
  s3[ 0 ] = 0.; // RightBin.zmin

  s3[ 1 ] = 1.; // RightBin.xmax
  s3[ 2 ] = 1.; // RightBin.ymax
  s3[ 3 ] = 1.; // RightBin.zmax

  RayType ray2d;
  ray2d[ 0 ] = -1.0;
  ray2d[ 1 ] = -1.0;
  ray2d[ 2 ] =  1.0;
  ray2d[ 3 ] =  1.0;
  EXPECT_TRUE( TraversalPredicates2D::rayIntersectsRightBin( ray2d, s2, s3 ) );

  // flip normal
  ray2d[ 2 ] = -1.0;
  ray2d[ 3 ] = -1.0;
  EXPECT_FALSE( TraversalPredicates2D::rayIntersectsRightBin( ray2d, s2, s3 ) );

  RayType ray3d;
  ray3d[ 0 ] = -1.0;
  ray3d[ 1 ] = -1.0;
  ray3d[ 2 ] = -1.0;
  ray3d[ 3 ] = 1.0;
  ray3d[ 4 ] = 1.0;
  ray3d[ 5 ] = 1.0;
  EXPECT_TRUE( TraversalPredicates3D::rayIntersectsRightBin( ray3d, s2, s3) );

  // flip normal
  ray3d[ 3 ] = -1.0;
  ray3d[ 4 ] = -1.0;
  ray3d[ 5 ] = -1.0;
  EXPECT_FALSE( TraversalPredicates3D::rayIntersectsRightBin( ray3d, s2, s3 ) );
}

//------------------------------------------------------------------------------
TEST( spin_bvh, traversal_predicates_pointInLeftBin )
{
  namespace bvh               = axom::spin::internal::linear_bvh;
  using TraversalPredicates2D = bvh::TraversalPredicates< 2, double >;
  using TraversalPredicates3D = bvh::TraversalPredicates< 3, double >;
  using VecType               = bvh::Vec< double, 4 >;

  VecType in_point, out_point, s1, s2;

  in_point[ 0 ]  = in_point[ 1 ]  = in_point[ 2 ]  = 0.5;
  out_point[ 0 ] = out_point[ 1 ] = out_point[ 2 ] = 1.5;

  s1[ 0 ] = 0.; // LeftBin.xmin
  s1[ 1 ] = 0.; // LeftBin.ymin
  s1[ 2 ] = 0.; // LeftBin.zmin

  s1[ 3 ] = 1.; // LeftBin.xmax
  s2[ 0 ] = 1.; // LeftBin.ymax
  s2[ 1 ] = 1.; // LeftBin.zmax

  EXPECT_TRUE( TraversalPredicates2D::pointInLeftBin( in_point, s1, s2 ) );
  EXPECT_TRUE( TraversalPredicates3D::pointInLeftBin( in_point, s1, s2 ) );

  EXPECT_FALSE( TraversalPredicates2D::pointInLeftBin( out_point, s1, s2 ) );
  EXPECT_FALSE( TraversalPredicates3D::pointInLeftBin( out_point, s1, s2 ) );
}

//------------------------------------------------------------------------------
TEST( spin_bvh, traversal_predicates_pointInRightBin )
{
  namespace bvh               = axom::spin::internal::linear_bvh;
  using TraversalPredicates2D = bvh::TraversalPredicates< 2, double >;
  using TraversalPredicates3D = bvh::TraversalPredicates< 3, double >;
  using VecType               = bvh::Vec< double, 4 >;

  VecType in_point, out_point, s2, s3;

  in_point[ 0 ]  = in_point[ 1 ]  = in_point[ 2 ]  = 0.5;
  out_point[ 0 ] = out_point[ 1 ] = out_point[ 2 ] = 1.5;

  s2[ 2 ] = 0.; // RightBin.xmin
  s2[ 3 ] = 0.; // RightBin.ymin
  s3[ 0 ] = 0.; // RightBin.zmin

  s3[ 1 ] = 1.; // RightBin.xmax
  s3[ 2 ] = 1.; // RightBin.ymax
  s3[ 3 ] = 1.; // RightBin.zmax

  EXPECT_TRUE( TraversalPredicates2D::pointInRightBin( in_point, s2, s3 ) );
  EXPECT_TRUE( TraversalPredicates3D::pointInRightBin( in_point, s2, s3 ) );

  EXPECT_FALSE( TraversalPredicates2D::pointInRightBin( out_point, s2, s3 ) );
  EXPECT_FALSE( TraversalPredicates3D::pointInRightBin( out_point, s2, s3 ) );
}

//------------------------------------------------------------------------------
TEST( spin_bvh, contruct2D_sequential )
{
  check_build_bvh2d< axom::SEQ_EXEC, float >( );
  check_build_bvh2d< axom::SEQ_EXEC, double >( );
}

//------------------------------------------------------------------------------
TEST( spin_bvh, contruct3D_sequential )
{
  check_build_bvh3d< axom::SEQ_EXEC, float >( );
  check_build_bvh3d< axom::SEQ_EXEC, double >( );
}

//------------------------------------------------------------------------------
TEST( spin_bvh, find_3d_sequential )
{
  check_find3d< axom::SEQ_EXEC, float >( );
  check_find3d< axom::SEQ_EXEC, double >( );
}

//------------------------------------------------------------------------------
TEST( spin_bvh, find_2d_sequential )
{
  check_find2d< axom::SEQ_EXEC, float >( );
  check_find2d< axom::SEQ_EXEC, double >( );
}

//------------------------------------------------------------------------------
TEST( spin_bvh, single_box2d_sequential )
{
  check_single_box2d< axom::SEQ_EXEC, float >( );
  check_single_box2d< axom::SEQ_EXEC, double >( );
}

//------------------------------------------------------------------------------
TEST( spin_bvh, single_box3d_sequential )
{
  check_single_box3d< axom::SEQ_EXEC, float >( );
  check_single_box3d< axom::SEQ_EXEC, double >( );
}

//------------------------------------------------------------------------------
#ifdef AXOM_USE_OPENMP

TEST( spin_bvh, contruct2D_omp )
{
  check_build_bvh2d< axom::OMP_EXEC, float >( );
  check_build_bvh2d< axom::OMP_EXEC, double >( );
}

//------------------------------------------------------------------------------
TEST( spin_bvh, contruct3D_omp )
{
  check_build_bvh3d< axom::OMP_EXEC, float >( );
  check_build_bvh3d< axom::OMP_EXEC, double >( );
}

//------------------------------------------------------------------------------
TEST( spin_bvh, find_3d_omp )
{
  check_find3d< axom::OMP_EXEC, float >( );
  check_find3d< axom::OMP_EXEC, double >( );
}

//------------------------------------------------------------------------------
TEST( spin_bvh, find_2d_omp )
{
  check_find2d< axom::OMP_EXEC, float >( );
  check_find2d< axom::OMP_EXEC, double >( );
}

//------------------------------------------------------------------------------
TEST( spin_bvh, single_box2d_omp )
{
  check_single_box2d< axom::OMP_EXEC, float >( );
  check_single_box2d< axom::OMP_EXEC, double >( );
}

//------------------------------------------------------------------------------
TEST( spin_bvh, single_box3d_omp )
{
  check_single_box3d< axom::OMP_EXEC, float >( );
  check_single_box3d< axom::OMP_EXEC, double >( );
}

#endif

//------------------------------------------------------------------------------
#ifdef AXOM_USE_CUDA

AXOM_CUDA_TEST( spin_bvh, contruct2D_cuda )
{
  constexpr int BLOCK_SIZE = 256;
  using exec  = axom::CUDA_EXEC< BLOCK_SIZE >;

  check_build_bvh2d< exec, float >( );
  check_build_bvh2d< exec, double >( );
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST( spin_bvh, contruct3D_cuda )
{
  constexpr int BLOCK_SIZE = 256;
  using exec  = axom::CUDA_EXEC< BLOCK_SIZE >;

  check_build_bvh3d< exec, float >( );
  check_build_bvh3d< exec, double >( );
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST( spin_bvh, find_3d_cuda )
{
  constexpr int BLOCK_SIZE = 256;
  using exec  = axom::CUDA_EXEC< BLOCK_SIZE >;

  check_find3d< exec, float >( );
  check_find3d< exec, double >( );
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST( spin_bvh, find_2d_cuda )
{
  constexpr int BLOCK_SIZE = 256;
  using exec  = axom::CUDA_EXEC< BLOCK_SIZE >;

  check_find2d< exec, float >( );
  check_find2d< exec, double >( );
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST( spin_bvh, single_box2d_cuda )
{
  constexpr int BLOCK_SIZE = 256;
  using exec  = axom::CUDA_EXEC< BLOCK_SIZE >;

  check_single_box2d< exec, float >( );
  check_single_box2d< exec, double >( );
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST( spin_bvh, single_box3d_cuda )
{
  constexpr int BLOCK_SIZE = 256;
  using exec  = axom::CUDA_EXEC< BLOCK_SIZE >;

  check_single_box3d< exec, float >( );
  check_single_box3d< exec, double >( );
}


#endif

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

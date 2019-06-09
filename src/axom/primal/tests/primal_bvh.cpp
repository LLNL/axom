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

// axom/core includes
#include "axom/config.hpp"
#include "axom/core/Types.hpp"
#include "axom/core/numerics/Matrix.hpp"
#include "axom/core/memory_management.hpp"

// axom/primal includes
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/spatial_acceleration/ExecutionSpace.hpp"
#include "axom/primal/spatial_acceleration/linear_bvh/range.hpp"
#include "axom/primal/spatial_acceleration/BVH.hpp"
#include "axom/primal/spatial_acceleration/UniformGrid.hpp"

// axom/mint includes
#include "axom/mint/mesh/Mesh.hpp"
#include "axom/mint/mesh/UniformMesh.hpp"
#include "axom/mint/execution/interface.hpp"
#include "axom/mint/utils/vtk_utils.hpp"

// gtest includes
#include "gtest/gtest.h"

using namespace axom;
namespace xargs  = mint::xargs;
namespace policy = mint::policy;

//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------
namespace
{

//------------------------------------------------------------------------------
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

  using exec_policy = policy::serial;
  mint::for_all_cells< exec_policy, xargs::coords >(
      mesh, AXOM_LAMBDA( IndexType cellIdx,
                         numerics::Matrix< double >& coords,
                         const IndexType* AXOM_NOT_USED(nodeIds) )
  {

    primal::bvh::Range< double > xrange;
    primal::bvh::Range< double > yrange;

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

  using exec_policy = policy::serial;
  mint::for_all_cells< exec_policy, xargs::coords >(
      mesh, AXOM_LAMBDA( IndexType cellIdx,
                         numerics::Matrix< double >& coords,
                         const IndexType* AXOM_NOT_USED(nodeIds) )
  {

    primal::bvh::Range< double > xrange;
    primal::bvh::Range< double > yrange;
    primal::bvh::Range< double > zrange;

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
template < typename ExecSpace, typename FloatType >
void check_build_bvh2d( )
{
  constexpr int NUM_BOXES = 2;
  constexpr int NDIMS     = 2;

  umpire::Allocator current_allocator = axom::getDefaultAllocator();
  axom::setDefaultAllocator( primal::execution_space<ExecSpace>::allocatorID());

  FloatType* boxes = axom::allocate< FloatType >( 8 );
  boxes[ 0 ] = boxes[ 1 ] = 0.;
  boxes[ 2 ] = boxes[ 3 ] = 1.;
  boxes[ 4 ] = boxes[ 5 ] = 1.;
  boxes[ 6 ] = boxes[ 7 ] = 2.;

  primal::BVH< NDIMS, ExecSpace, FloatType > bvh( boxes, NUM_BOXES );
  bvh.build( );

  FloatType lo[ NDIMS ];
  FloatType hi[ NDIMS ];
  bvh.getBounds( lo, hi );

  for ( int idim=0; idim < NDIMS; ++idim )
  {
    EXPECT_DOUBLE_EQ( lo[ idim ], 0.0 );
    EXPECT_DOUBLE_EQ( hi[ idim ], 2.0 );
  }

  axom::deallocate( boxes );
  axom::setDefaultAllocator( current_allocator );
}

//------------------------------------------------------------------------------
template < typename ExecSpace, typename FloatType >
void check_build_bvh3d( )
{
  constexpr int NUM_BOXES = 2;
  constexpr int NDIMS     = 3;

  umpire::Allocator current_allocator = axom::getDefaultAllocator();
  axom::setDefaultAllocator( primal::execution_space<ExecSpace>::allocatorID());

  FloatType* boxes = axom::allocate< FloatType >( 12 );
  boxes[ 0 ] = boxes[  1 ] = boxes[  2 ] = 0.;
  boxes[ 3 ] = boxes[  4 ] = boxes[  5 ] = 1.;
  boxes[ 6 ] = boxes[  7 ] = boxes[  8 ] = 1.;
  boxes[ 9 ] = boxes[ 10 ] = boxes[ 11 ] = 2.;

  primal::BVH< NDIMS, ExecSpace, FloatType > bvh( boxes, NUM_BOXES );
  bvh.build( );

  FloatType lo[ NDIMS ];
  FloatType hi[ NDIMS ];
  bvh.getBounds( lo, hi );

  for ( int idim=0; idim < NDIMS; ++idim )
  {
    EXPECT_DOUBLE_EQ( lo[ idim ], 0.0 );
    EXPECT_DOUBLE_EQ( hi[ idim ], 2.0 );
  }

  axom::deallocate( boxes );
  axom::setDefaultAllocator( current_allocator );
}

//------------------------------------------------------------------------------
template < typename ExecSpace, typename FloatType >
void check_find3d( )
{
  constexpr int NDIMS   = 3;
  constexpr IndexType N = 4;

  umpire::Allocator current_allocator = axom::getDefaultAllocator();
  axom::setDefaultAllocator( primal::execution_space<ExecSpace>::allocatorID());

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
  primal::BVH< NDIMS, ExecSpace, FloatType > bvh( aabbs, ncells );
  bvh.build( );

  FloatType min[ NDIMS ];
  FloatType max[ NDIMS ];
  bvh.getBounds( min, max );
  for ( int i=0; i < NDIMS; ++i )
  {
    EXPECT_DOUBLE_EQ( min[ i ], lo[ i ] );
    EXPECT_DOUBLE_EQ( max[ i ], hi[ i ] );
  }

  // traverse the BVH to find the candidates for all the centroids
  IndexType* offsets    = axom::allocate< IndexType >( ncells );
  IndexType* counts     = axom::allocate< IndexType >( ncells );
  IndexType* candidates = nullptr;
  bvh.find( offsets, counts, candidates, ncells, xc, yc, zc );

  EXPECT_TRUE( candidates != nullptr );

  primal::UniformGrid< IndexType, NDIMS > ug( lo, hi, res );

  for ( IndexType i=0; i < ncells; ++i )
  {
    PointType q            = PointType::make_point( xc[ i ],yc[ i ],zc[ i ] );
    const int donorCellIdx = ug.getBinIndex( q );
    EXPECT_EQ( counts[ i ], 1 );
    EXPECT_EQ( donorCellIdx, candidates[ offsets[ i ] ] );
  } // END for all cell centroids

  axom::deallocate( offsets );
  axom::deallocate( candidates );
  axom::deallocate( counts );
  axom::deallocate( aabbs );

  axom::setDefaultAllocator( current_allocator );
}

//------------------------------------------------------------------------------
template < typename ExecSpace, typename FloatType >
void check_find2d( )
{
  constexpr int NDIMS   = 2;
  constexpr IndexType N = 4;

  umpire::Allocator current_allocator = axom::getDefaultAllocator();
  axom::setDefaultAllocator( primal::execution_space<ExecSpace>::allocatorID());

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
  primal::BVH< NDIMS, ExecSpace, FloatType > bvh( aabbs, ncells );
  bvh.build( );

  FloatType min[ NDIMS ];
  FloatType max[ NDIMS ];
  bvh.getBounds( min, max );
  for ( int i=0; i < NDIMS; ++i )
  {
    EXPECT_DOUBLE_EQ( min[ i ], lo[ i ] );
    EXPECT_DOUBLE_EQ( max[ i ], hi[ i ] );
  }

  // traverse the BVH to find the candidates for all the centroids
  IndexType* offsets    = axom::allocate< IndexType >( ncells );
  IndexType* counts     = axom::allocate< IndexType >( ncells );
  IndexType* candidates = nullptr;
  bvh.find( offsets, counts, candidates, ncells, xc, yc );

  EXPECT_TRUE( candidates != nullptr );

  primal::UniformGrid< IndexType, NDIMS > ug( lo, hi, res );

  for ( IndexType i=0; i < ncells; ++i )
  {
    PointType q            = PointType::make_point( xc[ i ],yc[ i ] );
    const int donorCellIdx = ug.getBinIndex( q );
    EXPECT_EQ( counts[ i ], 1 );
    EXPECT_EQ( donorCellIdx, candidates[ offsets[ i ] ] );
  } // END for all cell centroids

  axom::deallocate( offsets );
  axom::deallocate( candidates );
  axom::deallocate( counts );
  axom::deallocate( aabbs );

  axom::setDefaultAllocator( current_allocator );
}

} /* end unnamed namespace */

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
TEST( primal_bvh, contruct2D_sequential )
{
  check_build_bvh2d< primal::SEQ_EXEC, float >( );
  check_build_bvh2d< primal::SEQ_EXEC, double >( );
}

//------------------------------------------------------------------------------
TEST( primal_bvh, contruct3D_sequential )
{
  check_build_bvh3d< primal::SEQ_EXEC, float >( );
  check_build_bvh3d< primal::SEQ_EXEC, double >( );
}

//------------------------------------------------------------------------------
TEST( primal_bvh, find_3d_sequential )
{
  check_find3d< primal::SEQ_EXEC, float >( );
  check_find3d< primal::SEQ_EXEC, double >( );
}

//------------------------------------------------------------------------------
TEST( primal_bvh, find_2d_sequential )
{
  check_find2d< primal::SEQ_EXEC, float >( );
  check_find2d< primal::SEQ_EXEC, double >( );
}

//------------------------------------------------------------------------------
#ifdef AXOM_USE_OPENMP

TEST( primal_bvh, contruct2D_omp )
{
  check_build_bvh2d< primal::OMP_EXEC, float >( );
  check_build_bvh2d< primal::OMP_EXEC, double >( );
}

//------------------------------------------------------------------------------
TEST( primal_bvh, contruct3D_omp )
{
  check_build_bvh3d< primal::OMP_EXEC, float >( );
  check_build_bvh3d< primal::OMP_EXEC, double >( );
}

//------------------------------------------------------------------------------
TEST( primal_bvh, find_3d_omp )
{
  check_find3d< primal::OMP_EXEC, float >( );
  check_find3d< primal::OMP_EXEC, double >( );
}

//------------------------------------------------------------------------------
TEST( primal_bvh, find_2d_omp )
{
  check_find2d< primal::OMP_EXEC, float >( );
  check_find2d< primal::OMP_EXEC, double >( );
}

#endif

//------------------------------------------------------------------------------
#ifdef AXOM_USE_CUDA

AXOM_CUDA_TEST( primal_bvh, contruct2D_cuda )
{
  constexpr int BLOCK_SIZE = 256;
  using exec  = primal::CUDA_EXEC< BLOCK_SIZE >;

  check_build_bvh2d< exec, float >( );
  check_build_bvh2d< exec, double >( );
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST( primal_bvh, contruct3D_cuda )
{
  constexpr int BLOCK_SIZE = 256;
  using exec  = primal::CUDA_EXEC< BLOCK_SIZE >;

  check_build_bvh3d< exec, float >( );
  check_build_bvh3d< exec, double >( );
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST( primal_bvh, find_3d_cuda )
{
  constexpr int BLOCK_SIZE = 256;
  using exec  = primal::CUDA_EXEC< BLOCK_SIZE >;

  check_find3d< exec, float >( );
  check_find3d< exec, double >( );
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST( primal_bvh, find_2d_cuda )
{
  constexpr int BLOCK_SIZE = 256;
  using exec  = primal::CUDA_EXEC< BLOCK_SIZE >;

  check_find2d< exec, float >( );
  check_find2d< exec, double >( );
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

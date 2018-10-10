/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
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

// Axom includes
#include "axom/config.hpp"                   // for compile-time definitions
#include "axom/core/utilities/Utilities.hpp" // for alloc() /free()

// Mint includes
#include "axom/mint/config.hpp"              // mint compile-time definitions
#include "axom/mint/execution/policy.hpp"    // mint execution policies/traits
#include "axom/mint/execution/interface.hpp" // for_all()

// Slic includes
#include "axom/slic.hpp" // for SLIC macros

#include "mint_test_utilities.hpp"

// gtest includes
#include "gtest/gtest.h" // for gtest

namespace axom
{
namespace mint
{

//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------
namespace
{

template < typename ExecPolicy, int MeshType >
void check_for_all_nodes_idx( int dimension )
{
  SLIC_INFO( "dimension=" << dimension << ", policy="
            << policy_traits< ExecPolicy >::name() << ", mesh_type="
            << internal::mesh_type_name< MeshType >::name() );

  constexpr int MAGIC_VAL = 42;

  const IndexType Ni = 20;
  const IndexType Nj = (dimension >= 2) ? 20 : -1;
  const IndexType Nk = (dimension == 3) ? 20 : -1;

  const double lo[] = { -10, -10, -10 };
  const double hi[] = {  10,  10,  10 };
  UniformMesh uniform_mesh( lo, hi, Ni, Nj, Nk );

  Mesh* test_mesh = internal::create_mesh< MeshType >( uniform_mesh );
  EXPECT_TRUE( test_mesh != nullptr );

  int* field = test_mesh->createField< int >( "n1", NODE_CENTERED );

  for_all_nodes< ExecPolicy >(
    test_mesh, AXOM_LAMBDA(IndexType nodeIdx)
    {
      field[ nodeIdx ] = MAGIC_VAL;
    } );

  const IndexType numNodes = test_mesh->getNumberOfNodes();
  for ( IndexType inode=0 ; inode < numNodes ; ++inode )
  {
    EXPECT_EQ( field[ inode ], MAGIC_VAL );
  }

  delete test_mesh;
  test_mesh = nullptr;
}

//------------------------------------------------------------------------------
template < typename ExecPolicy, int MeshType >
void check_for_all_nodes_ij( )
{
  SLIC_INFO( "policy=" << policy_traits< ExecPolicy >::name() << ", mesh_type="
            << internal::mesh_type_name< MeshType >::name() );
  
  constexpr IndexType N = 20;
  const double lo[] = { -10, -10 };
  const double hi[] = {  10,  10 };
  UniformMesh uniform_mesh( lo, hi, N, N );

  // STEP 0: create the test mesh
  Mesh* test_mesh = internal::create_mesh< MeshType >( uniform_mesh );
  EXPECT_TRUE( test_mesh != nullptr );

  const IndexType numNodes = test_mesh->getNumberOfNodes();

  IndexType* icoords = utilities::alloc< IndexType >( numNodes );
  IndexType* jcoords = utilities::alloc< IndexType >( numNodes );

  for_all_nodes< ExecPolicy, xargs::ij >( test_mesh,
                                                      AXOM_LAMBDA( IndexType
                                                                   nodeIdx,
                                                                   IndexType i,
                                                                   IndexType j )
    {
      icoords[ nodeIdx ] = i;
      jcoords[ nodeIdx ] = j;
    } );

  IndexType inode = 0;
  for ( IndexType j=0 ; j < N ; ++j )
  {
    for ( IndexType i=0 ; i < N ; ++i )
    {

      EXPECT_EQ( icoords[ inode ], i );
      EXPECT_EQ( jcoords[ inode ], j );
      ++inode;

    } // END for all i
  } // END for all j

  delete test_mesh;
  test_mesh = nullptr;

  utilities::free( icoords );
  utilities::free( jcoords );
}

//------------------------------------------------------------------------------
template < typename ExecPolicy, int MeshType >
void check_for_all_nodes_ijk( )
{
  SLIC_INFO( "policy=" << policy_traits< ExecPolicy >::name() << ", mesh_type="
            << internal::mesh_type_name< MeshType >::name() );

  constexpr IndexType N = 20;
  const double lo[] = { -10, -10, -10 };
  const double hi[] = {  10,  10,  10 };
  UniformMesh uniform_mesh( lo, hi, N, N, N );

  // STEP 0: create the test mesh
  Mesh* test_mesh = internal::create_mesh< MeshType >( uniform_mesh );
  EXPECT_TRUE( test_mesh != nullptr );

  const IndexType numNodes = test_mesh->getNumberOfNodes();

  IndexType* icoords = utilities::alloc< IndexType >( numNodes );
  IndexType* jcoords = utilities::alloc< IndexType >( numNodes );
  IndexType* kcoords = utilities::alloc< IndexType >( numNodes );

  for_all_nodes< ExecPolicy, xargs::ijk >( test_mesh,
                                                       AXOM_LAMBDA( IndexType
                                                                    nodeIdx,
                                                                    IndexType i,
                                                                    IndexType j,
                                                                    IndexType k )
    {
      icoords[ nodeIdx ] = i;
      jcoords[ nodeIdx ] = j;
      kcoords[ nodeIdx ] = k;
    } );

  IndexType inode = 0;
  for ( IndexType k=0 ; k < N ; ++k )
  {
    for ( IndexType j=0 ; j < N ; ++j )
    {
      for ( IndexType i=0 ; i < N ; ++i )
      {

        EXPECT_EQ( icoords[ inode ], i );
        EXPECT_EQ( jcoords[ inode ], j );
        EXPECT_EQ( kcoords[ inode ], k );
        ++inode;

      } // END for all i
    } // END for all j
  } // END for all k

  delete test_mesh;
  test_mesh = nullptr;

  utilities::free( icoords );
  utilities::free( jcoords );
  utilities::free( kcoords );
}

//------------------------------------------------------------------------------
template < typename ExecPolicy, int MeshType >
void check_for_all_nodes_xyz( )
{
  SLIC_INFO( "policy=" << policy_traits< ExecPolicy >::name() << ", mesh_type="
            << internal::mesh_type_name< MeshType >::name() );

  constexpr IndexType N = 20;
  const double lo[] = { -10, -10, -10 };
  const double hi[] = {  10,  10,  10 };
  UniformMesh uniform_mesh( lo, hi, N, N, N );

  // STEP 0: create the test mesh
  Mesh* test_mesh = internal::create_mesh< MeshType >( uniform_mesh );
  EXPECT_TRUE( test_mesh != nullptr );

  // STEP 1: generate test coordinate arrays
  const IndexType numNodes = test_mesh->getNumberOfNodes();
  double* x = axom::utilities::alloc< double >( numNodes );
  double* y = axom::utilities::alloc< double >( numNodes );
  double* z = axom::utilities::alloc< double >( numNodes );
  for_all_nodes< ExecPolicy, xargs::xyz >( test_mesh,
    AXOM_LAMBDA( IndexType idx, double xx, double yy, double zz )
    {
      x[ idx ] = xx;
      y[ idx ] = yy;
      z[ idx ] = zz;
    }
  );

  // STEP 2:check coordinate arrays
  for ( int inode=0 ; inode < numNodes ; ++inode )
  {
    double node[ 3 ];
    test_mesh->getNode( inode, node );
    EXPECT_DOUBLE_EQ( x[ inode ], node[ X_COORDINATE ] );
    EXPECT_DOUBLE_EQ( y[ inode ], node[ Y_COORDINATE ] );
    EXPECT_DOUBLE_EQ( z[ inode ], node[ Z_COORDINATE ] );
  } // END for all nodes

  // STEP 3: clean up
  delete test_mesh;
  test_mesh = nullptr;

  axom::utilities::free( x );
  axom::utilities::free( y );
  axom::utilities::free( z );
}

//------------------------------------------------------------------------------
template < typename ExecPolicy, int MeshType >
void check_for_all_nodes_xy( )
{
  SLIC_INFO( "policy=" << policy_traits< ExecPolicy >::name() << ", mesh_type="
            << internal::mesh_type_name< MeshType >::name() );

  constexpr IndexType N = 20;
  const double lo[] = { -10, -10 };
  const double hi[] = {  10,  10 };
  UniformMesh uniform_mesh( lo, hi, N, N );

  // STEP 0: create the test mesh
  Mesh* test_mesh = internal::create_mesh< MeshType >( uniform_mesh );
  EXPECT_TRUE( test_mesh != nullptr );

  // STEP 1: generate test coordinate arrays
  const IndexType numNodes = test_mesh->getNumberOfNodes();
  double* x = axom::utilities::alloc< double >( numNodes );
  double* y = axom::utilities::alloc< double >( numNodes );
  for_all_nodes< ExecPolicy, xargs::xy >( test_mesh,
    AXOM_LAMBDA( IndexType idx, double xx, double yy )
    {
      x[ idx ] = xx;
      y[ idx ] = yy;
    }
  );

  // STEP 2:check coordinate arrays
  for ( int inode=0 ; inode < numNodes ; ++inode )
  {
    double node[ 2 ];
    test_mesh->getNode( inode, node );
    EXPECT_DOUBLE_EQ( x[ inode ], node[ X_COORDINATE ] );
    EXPECT_DOUBLE_EQ( y[ inode ], node[ Y_COORDINATE ] );
  } // END for all nodes

  // STEP 3: clean up
  delete test_mesh;
  test_mesh = nullptr;

  axom::utilities::free( x );
  axom::utilities::free( y );
}

//------------------------------------------------------------------------------
template < typename ExecPolicy, int MeshType >
void check_for_all_nodes_x( )
{
  SLIC_INFO( "policy=" << policy_traits< ExecPolicy >::name() << ", mesh_type="
            << internal::mesh_type_name< MeshType >::name() );

  constexpr IndexType Ni = 20;
  const double lo[] = { -10 };
  const double hi[] = {  10 };
  UniformMesh uniform_mesh( lo, hi, Ni );

  // STEP 0: create the test mesh
  Mesh* test_mesh = internal::create_mesh< MeshType >( uniform_mesh );
  EXPECT_TRUE( test_mesh != nullptr );

  if ( MeshType != PARTICLE_MESH )
  {
    EXPECT_EQ( test_mesh->getNumberOfCells(), uniform_mesh.getNumberOfCells() );
  }

  // STEP 1: generate test coordinate arrays
  const IndexType numNodes = uniform_mesh.getNumberOfNodes();
  double* x = axom::utilities::alloc< double >( numNodes );
  for_all_nodes< ExecPolicy, xargs::x >(
    test_mesh, AXOM_LAMBDA( IndexType idx, double xx )
    {
      x[ idx ] = xx;
    } );

  // STEP 2:check coordinate arrays
  for ( int inode=0 ; inode < numNodes ; ++inode )
  {
    double node[ 1 ];
    uniform_mesh.getNode( inode, node );
    EXPECT_DOUBLE_EQ( x[ inode ], node[ X_COORDINATE ] );
  } // END for all nodes

  // STEP 3: clean up
  delete test_mesh;
  test_mesh = nullptr;
  axom::utilities::free( x );
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------

TEST( mint_execution_node_traversals, for_all_nodes_xyz )
{
  using seq_exec = policy::serial;
  check_for_all_nodes_xyz< seq_exec, STRUCTURED_UNIFORM_MESH >();
  check_for_all_nodes_xyz< seq_exec, STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_nodes_xyz< seq_exec, STRUCTURED_RECTILINEAR_MESH >();
  check_for_all_nodes_xyz< seq_exec, PARTICLE_MESH >();
  check_for_all_nodes_xyz< seq_exec, UNSTRUCTURED_SINGLE_SHAPE_MESH >();
  check_for_all_nodes_xyz< seq_exec, UNSTRUCTURED_MIXED_SHAPE_MESH >();

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && \
  defined(RAJA_ENABLE_OPENMP)

  using openmp_exec = policy::parallel_cpu;
  check_for_all_nodes_xyz< openmp_exec, STRUCTURED_UNIFORM_MESH >();
  check_for_all_nodes_xyz< openmp_exec, STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_nodes_xyz< openmp_exec, STRUCTURED_RECTILINEAR_MESH >();
  check_for_all_nodes_xyz< openmp_exec, PARTICLE_MESH >();
  check_for_all_nodes_xyz< openmp_exec, UNSTRUCTURED_SINGLE_SHAPE_MESH >();
  check_for_all_nodes_xyz< openmp_exec, UNSTRUCTURED_MIXED_SHAPE_MESH >();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && \
  defined(RAJA_ENABLE_CUDA)

  using cuda_exec = policy::parallel_gpu;
  check_for_all_nodes_xyz< cuda_exec, STRUCTURED_UNIFORM_MESH >();
  check_for_all_nodes_xyz< cuda_exec, STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_nodes_xyz< cuda_exec, STRUCTURED_RECTILINEAR_MESH >();
  check_for_all_nodes_xyz< cuda_exec, PARTICLE_MESH >();
  check_for_all_nodes_xyz< cuda_exec, UNSTRUCTURED_SINGLE_SHAPE_MESH >();
  check_for_all_nodes_xyz< cuda_exec, UNSTRUCTURED_MIXED_SHAPE_MESH >();

#endif
}

//------------------------------------------------------------------------------
TEST( mint_execution_node_traversals, for_all_nodes_xy )
{
  using seq_exec = policy::serial;
  check_for_all_nodes_xy< seq_exec, STRUCTURED_UNIFORM_MESH >();
  check_for_all_nodes_xy< seq_exec, STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_nodes_xy< seq_exec, STRUCTURED_RECTILINEAR_MESH >();
  check_for_all_nodes_xy< seq_exec, PARTICLE_MESH >();
  check_for_all_nodes_xy< seq_exec, UNSTRUCTURED_SINGLE_SHAPE_MESH >();
  check_for_all_nodes_xy< seq_exec, UNSTRUCTURED_MIXED_SHAPE_MESH >();

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && \
  defined(RAJA_ENABLE_OPENMP)

  using openmp_exec = policy::parallel_cpu;
  check_for_all_nodes_xy< openmp_exec, STRUCTURED_UNIFORM_MESH >();
  check_for_all_nodes_xy< openmp_exec, STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_nodes_xy< openmp_exec, STRUCTURED_RECTILINEAR_MESH >();
  check_for_all_nodes_xy< openmp_exec, PARTICLE_MESH >();
  check_for_all_nodes_xy< openmp_exec, UNSTRUCTURED_SINGLE_SHAPE_MESH >();
  check_for_all_nodes_xy< openmp_exec, UNSTRUCTURED_MIXED_SHAPE_MESH >();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && \
  defined(RAJA_ENABLE_CUDA)

  using cuda_exec = policy::parallel_gpu;
  check_for_all_nodes_xy< cuda_exec, STRUCTURED_UNIFORM_MESH >();
  check_for_all_nodes_xy< cuda_exec, STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_nodes_xy< cuda_exec, STRUCTURED_RECTILINEAR_MESH >();
  check_for_all_nodes_xy< cuda_exec, PARTICLE_MESH >();
  check_for_all_nodes_xy< cuda_exec, UNSTRUCTURED_SINGLE_SHAPE_MESH >();
  check_for_all_nodes_xy< cuda_exec, UNSTRUCTURED_MIXED_SHAPE_MESH >();

#endif
}

//------------------------------------------------------------------------------
TEST( mint_execution_node_traversals, for_all_nodes_x )
{
  using seq_exec = policy::serial;
  check_for_all_nodes_x< seq_exec, STRUCTURED_UNIFORM_MESH >();
  check_for_all_nodes_x< seq_exec, STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_nodes_x< seq_exec, STRUCTURED_RECTILINEAR_MESH >();
  check_for_all_nodes_x< seq_exec, PARTICLE_MESH >();
  check_for_all_nodes_x< seq_exec, UNSTRUCTURED_SINGLE_SHAPE_MESH >();
  check_for_all_nodes_x< seq_exec, UNSTRUCTURED_MIXED_SHAPE_MESH >();

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && \
  defined(RAJA_ENABLE_OPENMP)

  using openmp_exec = policy::parallel_cpu;
  check_for_all_nodes_x< openmp_exec, STRUCTURED_UNIFORM_MESH >();
  check_for_all_nodes_x< openmp_exec, STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_nodes_x< openmp_exec, STRUCTURED_RECTILINEAR_MESH >();
  check_for_all_nodes_x< openmp_exec, PARTICLE_MESH >();
  check_for_all_nodes_x< openmp_exec, UNSTRUCTURED_SINGLE_SHAPE_MESH >();
  check_for_all_nodes_x< openmp_exec, UNSTRUCTURED_MIXED_SHAPE_MESH >();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && \
  defined(RAJA_ENABLE_CUDA)

  using cuda_exec = policy::parallel_gpu;
  check_for_all_nodes_x< cuda_exec, STRUCTURED_UNIFORM_MESH >();
  check_for_all_nodes_x< cuda_exec, STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_nodes_x< cuda_exec, STRUCTURED_RECTILINEAR_MESH >();
  check_for_all_nodes_x< cuda_exec, PARTICLE_MESH >();
  check_for_all_nodes_x< cuda_exec, UNSTRUCTURED_SINGLE_SHAPE_MESH >();
  check_for_all_nodes_x< cuda_exec, UNSTRUCTURED_MIXED_SHAPE_MESH >();

#endif

}

//------------------------------------------------------------------------------
TEST( mint_execution_node_traversals, for_all_nodes_ijk )
{
  using seq_exec = policy::serial;
  check_for_all_nodes_ijk< seq_exec, STRUCTURED_UNIFORM_MESH >();
  check_for_all_nodes_ijk< seq_exec, STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_nodes_ijk< seq_exec, STRUCTURED_RECTILINEAR_MESH >();

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && \
  defined(RAJA_ENABLE_OPENMP)

  using openmp_exec = policy::parallel_cpu;
  check_for_all_nodes_ijk< openmp_exec, STRUCTURED_UNIFORM_MESH >();
  check_for_all_nodes_ijk< openmp_exec, STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_nodes_ijk< openmp_exec, STRUCTURED_RECTILINEAR_MESH >();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && \
  defined(RAJA_ENABLE_CUDA)

  using cuda_exec = policy::parallel_gpu;
  check_for_all_nodes_ijk< cuda_exec, STRUCTURED_UNIFORM_MESH >();
  check_for_all_nodes_ijk< cuda_exec, STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_nodes_ijk< cuda_exec, STRUCTURED_RECTILINEAR_MESH >();

#endif
}

//------------------------------------------------------------------------------
TEST( mint_execution_node_traversals, for_all_nodes_ij )
{
  using seq_exec = policy::serial;
  check_for_all_nodes_ij< seq_exec, STRUCTURED_UNIFORM_MESH >();
  check_for_all_nodes_ij< seq_exec, STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_nodes_ij< seq_exec, STRUCTURED_RECTILINEAR_MESH >();

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && \
  defined(RAJA_ENABLE_OPENMP)

  using openmp_exec = policy::parallel_cpu;
  check_for_all_nodes_ij< openmp_exec, STRUCTURED_UNIFORM_MESH >();
  check_for_all_nodes_ij< openmp_exec, STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_nodes_ij< openmp_exec, STRUCTURED_RECTILINEAR_MESH >();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && \
  defined(RAJA_ENABLE_CUDA)

  using cuda_exec = policy::parallel_gpu;
  check_for_all_nodes_ij< cuda_exec, STRUCTURED_UNIFORM_MESH >();
  check_for_all_nodes_ij< cuda_exec, STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_nodes_ij< cuda_exec, STRUCTURED_RECTILINEAR_MESH >();

#endif
}

//------------------------------------------------------------------------------
TEST( mint_execution_node_traversals, for_all_nodes_index )
{
  constexpr int NDIMS = 3;
  for ( int i=1 ; i <= NDIMS ; ++i )
  {

    using seq_exec = policy::serial;
    check_for_all_nodes_idx< seq_exec, STRUCTURED_UNIFORM_MESH >(i);
    check_for_all_nodes_idx< seq_exec, STRUCTURED_CURVILINEAR_MESH >(i);
    check_for_all_nodes_idx< seq_exec, STRUCTURED_RECTILINEAR_MESH >(i);
    check_for_all_nodes_idx< seq_exec, PARTICLE_MESH >(i);
    check_for_all_nodes_idx< seq_exec, UNSTRUCTURED_SINGLE_SHAPE_MESH >(i);
    check_for_all_nodes_idx< seq_exec, UNSTRUCTURED_MIXED_SHAPE_MESH >(i);

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && \
    defined(RAJA_ENABLE_OPENMP)

    using omp_exec = policy::parallel_cpu;
    check_for_all_nodes_idx< omp_exec, STRUCTURED_UNIFORM_MESH >(i);
    check_for_all_nodes_idx< omp_exec, STRUCTURED_CURVILINEAR_MESH >(i);
    check_for_all_nodes_idx< omp_exec, STRUCTURED_RECTILINEAR_MESH >(i);
    check_for_all_nodes_idx< omp_exec, PARTICLE_MESH >(i);
    check_for_all_nodes_idx< omp_exec, UNSTRUCTURED_SINGLE_SHAPE_MESH >(i);
    check_for_all_nodes_idx< omp_exec, UNSTRUCTURED_MIXED_SHAPE_MESH >(i);

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && \
    defined(RAJA_ENABLE_CUDA)

    using cuda_exec = policy::parallel_gpu;
    check_for_all_nodes_idx< cuda_exec, STRUCTURED_UNIFORM_MESH >(i);
    check_for_all_nodes_idx< cuda_exec, STRUCTURED_CURVILINEAR_MESH >(i);
    check_for_all_nodes_idx< cuda_exec, STRUCTURED_RECTILINEAR_MESH >(i);
    check_for_all_nodes_idx< cuda_exec, PARTICLE_MESH >(i);
    check_for_all_nodes_idx< cuda_exec, UNSTRUCTURED_SINGLE_SHAPE_MESH >(i);
    check_for_all_nodes_idx< cuda_exec, UNSTRUCTURED_MIXED_SHAPE_MESH >(i);

#endif

  } // END for all dimensions
}

} /* namespace mint */
} /* namespace axom */

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

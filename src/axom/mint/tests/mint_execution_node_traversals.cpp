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
#include "axom/mint/execution/interface.hpp" // mint::for_all()

// Slic includes
#include "axom/slic.hpp" // for SLIC macros

#include "mint_test_utilities.hpp"

// gtest includes
#include "gtest/gtest.h" // for gtest

// namespace aliases
namespace mint      = axom::mint;
namespace policy    = mint::policy;
namespace utilities = axom::utilities;

//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------
namespace
{

template < typename ExecPolicy, int MeshType >
void check_for_all_nodes_idx( int dimension )
{
  SLIC_INFO( "dimension=" << dimension );
  SLIC_INFO( "checking [for_all_nodes_index] policy=" <<
             mint::policy_traits< ExecPolicy >::name() << " " <<
             "mesh_type=" << mesh_type_name< MeshType >::name() );

  using IndexType = mint::IndexType;

  constexpr int MAGIC_VAL = 42;

  const IndexType Ni = 20;
  const IndexType Nj = (dimension >= 2) ? 20 : -1;
  const IndexType Nk = (dimension == 3) ? 20 : -1;

  const double lo[] = { -10, -10, -10 };
  const double hi[] = {  10,  10,  10 };
  mint::UniformMesh uniform_mesh( lo, hi, Ni, Nj, Nk );

  mint::Mesh* test_mesh = create_mesh< MeshType >( uniform_mesh );
  EXPECT_TRUE( test_mesh != nullptr );

  int* field = test_mesh->createField< int >( "n1", mint::NODE_CENTERED );

  mint::for_all_nodes< ExecPolicy >(
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
  SLIC_INFO( "checking [for_all_nodes_ij] policy=" <<
             mint::policy_traits< ExecPolicy >::name() << " " <<
             "mesh_type=" << mesh_type_name< MeshType >::name() );
  using IndexType = mint::IndexType;

  constexpr IndexType N = 20;
  const double lo[] = { -10, -10 };
  const double hi[] = {  10,  10 };
  mint::UniformMesh uniform_mesh( lo, hi, N, N );

  // STEP 0: create the test mesh
  mint::Mesh* test_mesh = create_mesh< MeshType >( uniform_mesh );
  EXPECT_TRUE( test_mesh != nullptr );

  const IndexType numNodes = test_mesh->getNumberOfNodes();

  IndexType* icoords = utilities::alloc< mint::IndexType >( numNodes );
  IndexType* jcoords = utilities::alloc< mint::IndexType >( numNodes );

  mint::for_all_nodes< ExecPolicy, mint::xargs::ij >( test_mesh,
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
  SLIC_INFO( "checking [for_all_nodes_ijk] policy=" <<
             mint::policy_traits< ExecPolicy >::name() << " " <<
             "mesh_type=" << mesh_type_name< MeshType >::name() );
  using IndexType = mint::IndexType;

  constexpr IndexType N = 20;
  const double lo[] = { -10, -10, -10 };
  const double hi[] = {  10,  10,  10 };
  mint::UniformMesh uniform_mesh( lo, hi, N, N, N );

  // STEP 0: create the test mesh
  mint::Mesh* test_mesh = create_mesh< MeshType >( uniform_mesh );
  EXPECT_TRUE( test_mesh != nullptr );

  const IndexType numNodes = test_mesh->getNumberOfNodes();

  IndexType* icoords = utilities::alloc< mint::IndexType >( numNodes );
  IndexType* jcoords = utilities::alloc< mint::IndexType >( numNodes );
  IndexType* kcoords = utilities::alloc< mint::IndexType >( numNodes );

  mint::for_all_nodes< ExecPolicy, mint::xargs::ijk >( test_mesh,
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
  SLIC_INFO( "checking [for_all_nodes_xyz] policy=" <<
             mint::policy_traits< ExecPolicy >::name() << " " <<
             "mesh_type=" << mesh_type_name< MeshType >::name() );

  constexpr mint::IndexType N = 20;
  const double lo[] = { -10, -10, -10 };
  const double hi[] = {  10,  10,  10 };
  mint::UniformMesh uniform_mesh( lo, hi, N, N, N );

  // STEP 0: create the test mesh
  mint::Mesh* test_mesh = create_mesh< MeshType >( uniform_mesh );
  EXPECT_TRUE( test_mesh != nullptr );

  if ( MeshType != mint::PARTICLE_MESH )
  {
    EXPECT_EQ( test_mesh->getNumberOfCells(), uniform_mesh.getNumberOfCells() );
  }

  // STEP 1: generate test coordinate arrays
  const mint::IndexType numNodes = uniform_mesh.getNumberOfNodes();
  double* x = axom::utilities::alloc< double >( numNodes );
  double* y = axom::utilities::alloc< double >( numNodes );
  double* z = axom::utilities::alloc< double >( numNodes );
  mint::for_all_nodes< ExecPolicy, mint::xargs::xyz >(test_mesh,
                                                      AXOM_LAMBDA( mint::
                                                                   IndexType idx,
                                                                   double xx,
                                                                   double yy,
                                                                   double zz )
    {
      x[ idx ] = xx;
      y[ idx ] = yy;
      z[ idx ] = zz;
    } );

  // STEP 2:check coordinate arrays
  for ( int inode=0 ; inode < numNodes ; ++inode )
  {
    double node[ 3 ];
    uniform_mesh.getNode( inode, node );
    EXPECT_DOUBLE_EQ( x[ inode ], node[ mint::X_COORDINATE ] );
    EXPECT_DOUBLE_EQ( y[ inode ], node[ mint::Y_COORDINATE ] );
    EXPECT_DOUBLE_EQ( z[ inode ], node[ mint::Z_COORDINATE ] );
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
  SLIC_INFO( "checking [for_all_nodes_xy] policy=" <<
             mint::policy_traits< ExecPolicy >::name() << " " <<
             "mesh_type=" << mesh_type_name< MeshType >::name() );

  constexpr mint::IndexType N = 20;
  const double lo[] = { -10, -10 };
  const double hi[] = {  10,  10 };
  mint::UniformMesh uniform_mesh( lo, hi, N, N );

  // STEP 0: create the test mesh
  mint::Mesh* test_mesh = create_mesh< MeshType >( uniform_mesh );
  EXPECT_TRUE( test_mesh != nullptr );

  if ( MeshType != mint::PARTICLE_MESH )
  {
    EXPECT_EQ( test_mesh->getNumberOfCells(), uniform_mesh.getNumberOfCells() );
  }

  // STEP 1: generate test coordinate arrays
  const mint::IndexType numNodes = uniform_mesh.getNumberOfNodes();
  double* x = axom::utilities::alloc< double >( numNodes );
  double* y = axom::utilities::alloc< double >( numNodes );
  mint::for_all_nodes< ExecPolicy, mint::xargs::xy >(
    test_mesh, AXOM_LAMBDA( mint::IndexType idx, double xx, double yy )
    {
      x[ idx ] = xx;
      y[ idx ] = yy;
    } );

  // STEP 2:check coordinate arrays
  for ( int inode=0 ; inode < numNodes ; ++inode )
  {
    double node[ 2 ];
    uniform_mesh.getNode( inode, node );
    EXPECT_DOUBLE_EQ( x[ inode ], node[ mint::X_COORDINATE ] );
    EXPECT_DOUBLE_EQ( y[ inode ], node[ mint::Y_COORDINATE ] );
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
  SLIC_INFO( "checking [for_all_nodes_x] policy=" <<
             mint::policy_traits< ExecPolicy >::name() << " " <<
             "mesh_type=" << mesh_type_name< MeshType >::name() );

  constexpr mint::IndexType Ni = 20;
  const double lo[] = { -10 };
  const double hi[] = {  10 };
  mint::UniformMesh uniform_mesh( lo, hi, Ni );

  // STEP 0: create the test mesh
  mint::Mesh* test_mesh = create_mesh< MeshType >( uniform_mesh );
  EXPECT_TRUE( test_mesh != nullptr );

  if ( MeshType != mint::PARTICLE_MESH )
  {
    EXPECT_EQ( test_mesh->getNumberOfCells(), uniform_mesh.getNumberOfCells() );
  }

  // STEP 1: generate test coordinate arrays
  const mint::IndexType numNodes = uniform_mesh.getNumberOfNodes();
  double* x = axom::utilities::alloc< double >( numNodes );
  mint::for_all_nodes< ExecPolicy, mint::xargs::x >(
    test_mesh, AXOM_LAMBDA( mint::IndexType idx, double xx )
    {
      x[ idx ] = xx;
    } );

  // STEP 2:check coordinate arrays
  for ( int inode=0 ; inode < numNodes ; ++inode )
  {
    double node[ 1 ];
    uniform_mesh.getNode( inode, node );
    EXPECT_DOUBLE_EQ( x[ inode ], node[ mint::X_COORDINATE ] );
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
  check_for_all_nodes_xyz< seq_exec, mint::STRUCTURED_UNIFORM_MESH >();
  check_for_all_nodes_xyz< seq_exec, mint::STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_nodes_xyz< seq_exec, mint::STRUCTURED_RECTILINEAR_MESH >();
  check_for_all_nodes_xyz< seq_exec, mint::PARTICLE_MESH >();
  check_for_all_nodes_xyz< seq_exec, mint::UNSTRUCTURED_MESH >();

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && \
  defined(RAJA_ENABLE_OPENMP)

  using openmp_exec = policy::parallel_cpu;
  check_for_all_nodes_xyz< openmp_exec, mint::STRUCTURED_UNIFORM_MESH >();
  check_for_all_nodes_xyz< openmp_exec, mint::STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_nodes_xyz< openmp_exec, mint::STRUCTURED_RECTILINEAR_MESH >();
  check_for_all_nodes_xyz< openmp_exec, mint::PARTICLE_MESH >();
  check_for_all_nodes_xyz< openmp_exec, mint::UNSTRUCTURED_MESH >();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && \
  defined(RAJA_ENABLE_CUDA)

  using cuda_exec = policy::parallel_gpu;
  check_for_all_nodes_xyz< cuda_exec, mint::STRUCTURED_UNIFORM_MESH >();
  check_for_all_nodes_xyz< cuda_exec, mint::STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_nodes_xyz< cuda_exec, mint::STRUCTURED_RECTILINEAR_MESH >();
  check_for_all_nodes_xyz< cuda_exec, mint::PARTICLE_MESH >();
  check_for_all_nodes_xyz< cuda_exec, mint::UNSTRUCTURED_MESH >();

#endif
}

//------------------------------------------------------------------------------
TEST( mint_execution_node_traversals, for_all_nodes_xy )
{
  using seq_exec = policy::serial;
  check_for_all_nodes_xy< seq_exec, mint::STRUCTURED_UNIFORM_MESH >();
  check_for_all_nodes_xy< seq_exec, mint::STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_nodes_xy< seq_exec, mint::STRUCTURED_RECTILINEAR_MESH >();
  check_for_all_nodes_xy< seq_exec, mint::PARTICLE_MESH >();
  check_for_all_nodes_xy< seq_exec, mint::UNSTRUCTURED_MESH >();

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && \
  defined(RAJA_ENABLE_OPENMP)

  using openmp_exec = policy::parallel_cpu;
  check_for_all_nodes_xy< openmp_exec, mint::STRUCTURED_UNIFORM_MESH >();
  check_for_all_nodes_xy< openmp_exec, mint::STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_nodes_xy< openmp_exec, mint::STRUCTURED_RECTILINEAR_MESH >();
  check_for_all_nodes_xy< openmp_exec, mint::PARTICLE_MESH >();
  check_for_all_nodes_xy< openmp_exec, mint::UNSTRUCTURED_MESH >();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && \
  defined(RAJA_ENABLE_CUDA)

  using cuda_exec = policy::parallel_gpu;
  check_for_all_nodes_xy< cuda_exec, mint::STRUCTURED_UNIFORM_MESH >();
  check_for_all_nodes_xy< cuda_exec, mint::STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_nodes_xy< cuda_exec, mint::STRUCTURED_RECTILINEAR_MESH >();
  check_for_all_nodes_xy< cuda_exec, mint::PARTICLE_MESH >();
  check_for_all_nodes_xy< cuda_exec, mint::UNSTRUCTURED_MESH >();

#endif
}

//------------------------------------------------------------------------------
TEST( mint_execution_node_traversals, for_all_nodes_x )
{
  using seq_exec = policy::serial;
  check_for_all_nodes_x< seq_exec, mint::STRUCTURED_UNIFORM_MESH >();
  check_for_all_nodes_x< seq_exec, mint::STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_nodes_x< seq_exec, mint::STRUCTURED_RECTILINEAR_MESH >();
  check_for_all_nodes_x< seq_exec, mint::PARTICLE_MESH >();
  check_for_all_nodes_x< seq_exec, mint::UNSTRUCTURED_MESH >();

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && \
  defined(RAJA_ENABLE_OPENMP)

  using openmp_exec = policy::parallel_cpu;
  check_for_all_nodes_x< openmp_exec, mint::STRUCTURED_UNIFORM_MESH >();
  check_for_all_nodes_x< openmp_exec, mint::STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_nodes_x< openmp_exec, mint::STRUCTURED_RECTILINEAR_MESH >();
  check_for_all_nodes_x< openmp_exec, mint::PARTICLE_MESH >();
  check_for_all_nodes_x< openmp_exec, mint::UNSTRUCTURED_MESH >();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && \
  defined(RAJA_ENABLE_CUDA)

  using cuda_exec = policy::parallel_gpu;
  check_for_all_nodes_x< cuda_exec, mint::STRUCTURED_UNIFORM_MESH >();
  check_for_all_nodes_x< cuda_exec, mint::STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_nodes_x< cuda_exec, mint::STRUCTURED_RECTILINEAR_MESH >();
  check_for_all_nodes_x< cuda_exec, mint::PARTICLE_MESH >();
  check_for_all_nodes_x< cuda_exec, mint::UNSTRUCTURED_MESH >();

#endif

}

//------------------------------------------------------------------------------
TEST( mint_execution_node_traversals, for_all_nodes_ijk )
{
  using seq_exec = policy::serial;
  check_for_all_nodes_ijk< seq_exec, mint::STRUCTURED_UNIFORM_MESH >();
  check_for_all_nodes_ijk< seq_exec, mint::STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_nodes_ijk< seq_exec, mint::STRUCTURED_RECTILINEAR_MESH >();

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && \
  defined(RAJA_ENABLE_OPENMP)

  using openmp_exec = policy::parallel_cpu;
  check_for_all_nodes_ijk< openmp_exec, mint::STRUCTURED_UNIFORM_MESH >();
  check_for_all_nodes_ijk< openmp_exec, mint::STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_nodes_ijk< openmp_exec, mint::STRUCTURED_RECTILINEAR_MESH >();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && \
  defined(RAJA_ENABLE_CUDA)

  using cuda_exec = policy::parallel_gpu;
  check_for_all_nodes_ijk< cuda_exec, mint::STRUCTURED_UNIFORM_MESH >();
  check_for_all_nodes_ijk< cuda_exec, mint::STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_nodes_ijk< cuda_exec, mint::STRUCTURED_RECTILINEAR_MESH >();

#endif
}

//------------------------------------------------------------------------------
TEST( mint_execution_node_traversals, for_all_nodes_ij )
{
  using seq_exec = policy::serial;
  check_for_all_nodes_ij< seq_exec, mint::STRUCTURED_UNIFORM_MESH >();
  check_for_all_nodes_ij< seq_exec, mint::STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_nodes_ij< seq_exec, mint::STRUCTURED_RECTILINEAR_MESH >();

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && \
  defined(RAJA_ENABLE_OPENMP)

  using openmp_exec = policy::parallel_cpu;
  check_for_all_nodes_ij< openmp_exec, mint::STRUCTURED_UNIFORM_MESH >();
  check_for_all_nodes_ij< openmp_exec, mint::STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_nodes_ij< openmp_exec, mint::STRUCTURED_RECTILINEAR_MESH >();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && \
  defined(RAJA_ENABLE_CUDA)

  using cuda_exec = policy::parallel_gpu;
  check_for_all_nodes_ij< cuda_exec, mint::STRUCTURED_UNIFORM_MESH >();
  check_for_all_nodes_ij< cuda_exec, mint::STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_nodes_ij< cuda_exec, mint::STRUCTURED_RECTILINEAR_MESH >();

#endif
}

//------------------------------------------------------------------------------
TEST( mint_execution_node_traversals, for_all_nodes_index )
{
  constexpr int NDIMS = 3;
  for ( int i=1 ; i <= NDIMS ; ++i )
  {

    using seq_exec = policy::serial;
    check_for_all_nodes_idx< seq_exec, mint::STRUCTURED_UNIFORM_MESH >(i);
    check_for_all_nodes_idx< seq_exec, mint::STRUCTURED_CURVILINEAR_MESH >(i);
    check_for_all_nodes_idx< seq_exec, mint::STRUCTURED_RECTILINEAR_MESH >(i);
    check_for_all_nodes_idx< seq_exec, mint::PARTICLE_MESH >(i);
    check_for_all_nodes_idx< seq_exec, mint::UNSTRUCTURED_MESH >(i);

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && \
    defined(RAJA_ENABLE_OPENMP)

    using omp_exec = policy::parallel_cpu;
    check_for_all_nodes_idx< omp_exec, mint::STRUCTURED_UNIFORM_MESH >(i);
    check_for_all_nodes_idx< omp_exec, mint::STRUCTURED_CURVILINEAR_MESH >(i);
    check_for_all_nodes_idx< omp_exec, mint::STRUCTURED_RECTILINEAR_MESH >(i);
    check_for_all_nodes_idx< omp_exec, mint::PARTICLE_MESH >(i);
    check_for_all_nodes_idx< omp_exec, mint::UNSTRUCTURED_MESH >(i);

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && \
    defined(RAJA_ENABLE_CUDA)

    using cuda_exec = policy::parallel_gpu;
    check_for_all_nodes_idx< cuda_exec, mint::STRUCTURED_UNIFORM_MESH >(i);
    check_for_all_nodes_idx< cuda_exec, mint::STRUCTURED_CURVILINEAR_MESH >(i);
    check_for_all_nodes_idx< cuda_exec, mint::STRUCTURED_RECTILINEAR_MESH >(i);
    check_for_all_nodes_idx< cuda_exec, mint::PARTICLE_MESH >(i);
    check_for_all_nodes_idx< cuda_exec, mint::UNSTRUCTURED_MESH >(i);

#endif

  } // END for all dimensions

}

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

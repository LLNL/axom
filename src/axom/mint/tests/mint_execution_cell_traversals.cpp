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
void check_for_all_cells_idx( int dimension )
{
  SLIC_INFO( "dimension=" << dimension );
  SLIC_INFO( "checking [for_all_cells_index] policy=" <<
             mint::policy_traits< ExecPolicy >::name() << " " <<
             "mesh_type=" << mesh_type_name< MeshType >::name() );

  using IndexType = mint::IndexType;

  constexpr int MAGIC_VAL = 42;

  const IndexType Ni = 20;
  const IndexType Nj = (dimension >= 2) ? Ni : -1;
  const IndexType Nk = (dimension == 3) ? Ni : -1;

  const double lo[] = { -10, -10, -10 };
  const double hi[] = {  10,  10,  10 };
  mint::UniformMesh uniform_mesh( lo, hi, Ni, Nj, Nk );

  mint::Mesh* test_mesh = nullptr;
  create_mesh< MeshType >( &uniform_mesh, test_mesh );
  EXPECT_TRUE( test_mesh != nullptr );
  EXPECT_EQ( test_mesh->getMeshType(), MeshType );
  EXPECT_EQ( test_mesh->getDimension(), uniform_mesh.getDimension() );
  EXPECT_EQ( test_mesh->getNumberOfNodes(), uniform_mesh.getNumberOfNodes() );
  EXPECT_EQ( test_mesh->getNumberOfCells(), uniform_mesh.getNumberOfCells() );

  int* field = test_mesh->createField< int >( "c1", mint::CELL_CENTERED );

  mint::for_all_cells< ExecPolicy >(
    test_mesh, AXOM_LAMBDA( IndexType cellIdx )
    {
      field[ cellIdx ] = MAGIC_VAL;
    } );

  const IndexType numCells = test_mesh->getNumberOfCells();
  for ( IndexType icell=0 ; icell < numCells ; ++icell )
  {
    EXPECT_EQ( field[ icell ], MAGIC_VAL );
  }

  delete test_mesh;
  test_mesh = nullptr;
}

//------------------------------------------------------------------------------
template < typename ExecPolicy, int MeshType >
void check_for_all_cells_ij( )
{
  SLIC_INFO( "checking [for_all_cells_ij] policy=" <<
             mint::policy_traits< ExecPolicy >::name() << " " <<
             "mesh_type=" << mesh_type_name< MeshType >::name() );
  using IndexType = mint::IndexType;

  constexpr IndexType N = 20;
  const double lo[] = { -10, -10 };
  const double hi[] = {  10,  10 };
  mint::UniformMesh uniform_mesh( lo, hi, N, N );

  // STEP 0: create the test mesh
  mint::Mesh* test_mesh = nullptr;
  create_mesh< MeshType >( &uniform_mesh, test_mesh );
  EXPECT_TRUE( test_mesh != nullptr );
  EXPECT_EQ( test_mesh->getMeshType(), MeshType );
  EXPECT_EQ( test_mesh->getDimension(), uniform_mesh.getDimension() );
  EXPECT_EQ( test_mesh->getNumberOfNodes(), uniform_mesh.getNumberOfNodes() );
  EXPECT_EQ( test_mesh->getNumberOfCells(), uniform_mesh.getNumberOfCells() );

  const IndexType numCells = test_mesh->getNumberOfCells();

  IndexType* icoords = utilities::alloc< mint::IndexType >( numCells );
  IndexType* jcoords = utilities::alloc< mint::IndexType >( numCells );

  mint::for_all_cells< ExecPolicy, mint::xargs::ij >( test_mesh,
                                                      AXOM_LAMBDA( IndexType
                                                                   cellIdx,
                                                                   IndexType i,
                                                                   IndexType j )
    {
      icoords[ cellIdx ] = i;
      jcoords[ cellIdx ] = j;
    } );

  IndexType icell = 0;
  for ( IndexType j=0 ; j < (N-1) ; ++j )
  {
    for ( IndexType i=0 ; i < (N-1) ; ++i )
    {
      EXPECT_EQ( icoords[ icell ], i );
      EXPECT_EQ( jcoords[ icell ], j );
      ++icell;

    } // END for all i
  } // END for all j

  delete test_mesh;
  test_mesh = nullptr;

  utilities::free( icoords );
  utilities::free( jcoords );
}

//------------------------------------------------------------------------------
template < typename ExecPolicy, int MeshType >
void check_for_all_cells_ijk( )
{
  SLIC_INFO( "checking [for_all_cells_ijk] policy=" <<
             mint::policy_traits< ExecPolicy >::name() << " " <<
             "mesh_type=" << mesh_type_name< MeshType >::name() );
  using IndexType = mint::IndexType;

  constexpr IndexType N = 20;
  const double lo[] = { -10, -10, -10 };
  const double hi[] = {  10,  10,  10 };
  mint::UniformMesh uniform_mesh( lo, hi, N, N, N );

  // STEP 0: create the test mesh
  mint::Mesh* test_mesh = nullptr;
  create_mesh< MeshType >( &uniform_mesh, test_mesh );
  EXPECT_TRUE( test_mesh != nullptr );
  EXPECT_EQ( test_mesh->getMeshType(), MeshType );
  EXPECT_EQ( test_mesh->getDimension(), uniform_mesh.getDimension() );
  EXPECT_EQ( test_mesh->getNumberOfNodes(), uniform_mesh.getNumberOfNodes() );
  EXPECT_EQ( test_mesh->getNumberOfCells(), uniform_mesh.getNumberOfCells() );

  const IndexType numCells = test_mesh->getNumberOfCells();

  IndexType* icoords = utilities::alloc< mint::IndexType >( numCells );
  IndexType* jcoords = utilities::alloc< mint::IndexType >( numCells );
  IndexType* kcoords = utilities::alloc< mint::IndexType >( numCells );

  mint::for_all_cells< ExecPolicy, mint::xargs::ijk >( test_mesh,
                                                       AXOM_LAMBDA( IndexType
                                                                    cellIdx,
                                                                    IndexType i,
                                                                    IndexType j,
                                                                    IndexType k )
    {
      icoords[ cellIdx ] = i;
      jcoords[ cellIdx ] = j;
      kcoords[ cellIdx ] = k;
    } );

  IndexType icell = 0;
  for ( IndexType k=0 ; k < (N-1) ; ++k )
  {
    for ( IndexType j=0 ; j < (N-1) ; ++j )
    {
      for ( IndexType i=0 ; i < (N-1) ; ++i )
      {
        EXPECT_EQ( icoords[ icell ], i );
        EXPECT_EQ( jcoords[ icell ], j );
        EXPECT_EQ( kcoords[ icell ], k );
        ++icell;

      } // END for all i
    } // END for all j
  } // END for all k

  delete test_mesh;
  test_mesh = nullptr;

  utilities::free( icoords );
  utilities::free( jcoords );
}

//------------------------------------------------------------------------------
template < typename ExecPolicy, int MeshType >
void check_for_all_cell_nodes( int dimension )
{
  SLIC_INFO( "dimension=" << dimension );
  SLIC_INFO( "checking [for_all_cells_index] policy=" <<
             mint::policy_traits< ExecPolicy >::name() << " " <<
             "mesh_type=" << mesh_type_name< MeshType >::name() );

  using IndexType = mint::IndexType;

  const IndexType Ni = 20;
  const IndexType Nj = (dimension >= 2) ? Ni : -1;
  const IndexType Nk = (dimension == 3) ? Ni : -1;

  const double lo[] = { -10, -10, -10 };
  const double hi[] = {  10,  10,  10 };
  mint::UniformMesh uniform_mesh( lo, hi, Ni, Nj, Nk );

  mint::Mesh* test_mesh = nullptr;
  create_mesh< MeshType >( &uniform_mesh, test_mesh );
  EXPECT_TRUE( test_mesh != nullptr );
  EXPECT_EQ( test_mesh->getMeshType(), MeshType );
  EXPECT_EQ( test_mesh->getDimension(), uniform_mesh.getDimension() );
  EXPECT_EQ( test_mesh->getNumberOfNodes(), uniform_mesh.getNumberOfNodes() );
  EXPECT_EQ( test_mesh->getNumberOfCells(), uniform_mesh.getNumberOfCells() );

  IndexType numCells        = test_mesh->getNumberOfCells();
  IndexType numNodesPerCell = test_mesh->getNumberOfCellNodes();
  IndexType* conn = utilities::alloc< IndexType >( numNodesPerCell*numCells );

  mint::for_all_cells< ExecPolicy, mint::xargs::nodeids >( test_mesh,
                                                           AXOM_LAMBDA(
                                                             IndexType cellIdx,
                                                             const IndexType*
                                                             ccon, IndexType N)
    {

      EXPECT_EQ( N, numNodesPerCell );

      for ( int i=0 ; i < N ; ++i )
      {
        conn[ cellIdx * N + i ] = ccon[ i ];
      } // END for all cell nodes

    } );

  for ( IndexType icell=0 ; icell < numCells ; ++icell )
  {
    const IndexType N = uniform_mesh.getNumberOfCellNodes( );
    EXPECT_EQ( N, numNodesPerCell );

    IndexType cellNodes[ 8 ];
    uniform_mesh.getCellNodeIDs( icell, cellNodes );

    for ( int i=0 ; i < N ; ++i )
    {
      EXPECT_EQ( conn[ icell * N + i ], cellNodes[ i ] );
    }

  }  // END for all cells

  /* clean up */
  delete test_mesh;
  test_mesh = nullptr;

  utilities::free( conn );
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------

TEST( mint_execution_cell_traversals, for_all_cells_nodeids )
{
  constexpr int NDIMS = 3;
  for ( int i=1 ; i <= NDIMS ; ++i )
  {

    using seq_exec = policy::serial;
    check_for_all_cell_nodes< seq_exec, mint::STRUCTURED_UNIFORM_MESH >(i);
    check_for_all_cell_nodes< seq_exec, mint::STRUCTURED_CURVILINEAR_MESH >(i);
    check_for_all_cell_nodes< seq_exec, mint::STRUCTURED_RECTILINEAR_MESH >(i);
    check_for_all_cell_nodes< seq_exec, mint::UNSTRUCTURED_MESH >(i);

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && \
    defined(RAJA_ENABLE_OPENMP)

    using omp_exec = policy::parallel_cpu;
    check_for_all_cell_nodes< omp_exec, mint::STRUCTURED_UNIFORM_MESH >(i);
    check_for_all_cell_nodes< omp_exec, mint::STRUCTURED_CURVILINEAR_MESH >(i);
    check_for_all_cell_nodes< omp_exec, mint::STRUCTURED_RECTILINEAR_MESH >(i);
    check_for_all_cell_nodes< omp_exec, mint::UNSTRUCTURED_MESH >(i);

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && \
    defined(RAJA_ENABLE_CUDA)

    using cuda_exec = policy::parallel_gpu;
    check_for_all_cell_nodes< cuda_exec, mint::STRUCTURED_UNIFORM_MESH >(i);
    check_for_all_cell_nodes< cuda_exec, mint::STRUCTURED_CURVILINEAR_MESH >(i);
    check_for_all_cell_nodes< cuda_exec, mint::STRUCTURED_RECTILINEAR_MESH >(i);
    check_for_all_cell_nodes< cuda_exec, mint::UNSTRUCTURED_MESH >(i);

#endif

  } // END for all dimensions
}

//------------------------------------------------------------------------------
TEST( mint_execution_cell_traversals, for_all_cells_ij )
{
  using seq_exec = policy::serial;
  check_for_all_cells_ij< seq_exec, mint::STRUCTURED_UNIFORM_MESH >();
  check_for_all_cells_ij< seq_exec, mint::STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_cells_ij< seq_exec, mint::STRUCTURED_RECTILINEAR_MESH >();

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && \
  defined(RAJA_ENABLE_OPENMP)

  using omp_exec = policy::parallel_cpu;
  check_for_all_cells_ij< omp_exec, mint::STRUCTURED_UNIFORM_MESH >();
  check_for_all_cells_ij< omp_exec, mint::STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_cells_ij< omp_exec, mint::STRUCTURED_RECTILINEAR_MESH >();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && \
  defined(RAJA_ENABLE_CUDA)

  using cuda_exec = policy::parallel_gpu;
  check_for_all_cells_ij< cuda_exec, mint::STRUCTURED_UNIFORM_MESH >();
  check_for_all_cells_ij< cuda_exec, mint::STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_cells_ij< cuda_exec, mint::STRUCTURED_RECTILINEAR_MESH >();

#endif

}

//------------------------------------------------------------------------------
TEST( mint_execution_cell_traversals, for_all_cells_ijk )
{
  using seq_exec = policy::serial;
  check_for_all_cells_ijk< seq_exec, mint::STRUCTURED_UNIFORM_MESH >();
  check_for_all_cells_ijk< seq_exec, mint::STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_cells_ijk< seq_exec, mint::STRUCTURED_RECTILINEAR_MESH >();

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && \
  defined(RAJA_ENABLE_OPENMP)

  using omp_exec = policy::parallel_cpu;
  check_for_all_cells_ijk< omp_exec, mint::STRUCTURED_UNIFORM_MESH >();
  check_for_all_cells_ijk< omp_exec, mint::STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_cells_ijk< omp_exec, mint::STRUCTURED_RECTILINEAR_MESH >();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && \
  defined(RAJA_ENABLE_CUDA)

  using cuda_exec = policy::parallel_gpu;
  check_for_all_cells_ijk< cuda_exec, mint::STRUCTURED_UNIFORM_MESH >();
  check_for_all_cells_ijk< cuda_exec, mint::STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_cells_ijk< cuda_exec, mint::STRUCTURED_RECTILINEAR_MESH >();

#endif

}

//------------------------------------------------------------------------------
TEST( mint_execution_cell_traversals, for_all_cells_index )
{
  constexpr int NDIMS = 3;
  for ( int i=1 ; i <= NDIMS ; ++i )
  {

    using seq_exec = policy::serial;
    check_for_all_cells_idx< seq_exec, mint::STRUCTURED_UNIFORM_MESH >(i);
    check_for_all_cells_idx< seq_exec, mint::STRUCTURED_CURVILINEAR_MESH >(i);
    check_for_all_cells_idx< seq_exec, mint::STRUCTURED_RECTILINEAR_MESH >(i);
    check_for_all_cells_idx< seq_exec, mint::UNSTRUCTURED_MESH >(i);

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && \
    defined(RAJA_ENABLE_OPENMP)

    using omp_exec = policy::parallel_cpu;
    check_for_all_cells_idx< omp_exec, mint::STRUCTURED_UNIFORM_MESH >(i);
    check_for_all_cells_idx< omp_exec, mint::STRUCTURED_CURVILINEAR_MESH >(i);
    check_for_all_cells_idx< omp_exec, mint::STRUCTURED_RECTILINEAR_MESH >(i);
    check_for_all_cells_idx< omp_exec, mint::UNSTRUCTURED_MESH >(i);

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && \
    defined(RAJA_ENABLE_CUDA)

    using cuda_exec = policy::parallel_gpu;
    check_for_all_cells_idx< cuda_exec, mint::STRUCTURED_UNIFORM_MESH >(i);
    check_for_all_cells_idx< cuda_exec, mint::STRUCTURED_CURVILINEAR_MESH >(i);
    check_for_all_cells_idx< cuda_exec, mint::STRUCTURED_RECTILINEAR_MESH >(i);
    check_for_all_cells_idx< cuda_exec, mint::UNSTRUCTURED_MESH >(i);

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

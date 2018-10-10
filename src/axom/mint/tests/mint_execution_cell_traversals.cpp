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
void check_for_all_cells_idx( int dimension )
{
  SLIC_INFO( "dimension=" << dimension << ", policy="
            << policy_traits< ExecPolicy >::name() << ", mesh_type="
            << internal::mesh_type_name< MeshType >::name() );

  const IndexType Ni = 20;
  const IndexType Nj = (dimension >= 2) ? Ni : -1;
  const IndexType Nk = (dimension == 3) ? Ni : -1;

  const double lo[] = { -10, -10, -10 };
  const double hi[] = {  10,  10,  10 };
  UniformMesh uniform_mesh( lo, hi, Ni, Nj, Nk );

  Mesh* test_mesh = internal::create_mesh< MeshType >( uniform_mesh );
  EXPECT_TRUE( test_mesh != nullptr );

  IndexType* field = test_mesh->createField< IndexType >( "c1", CELL_CENTERED );

  for_all_cells< ExecPolicy >( test_mesh, 
    AXOM_LAMBDA( IndexType cellID )
    {
      field[ cellID ] = cellID;
    }
  );

  const IndexType numCells = test_mesh->getNumberOfCells();
  for ( IndexType cellID = 0 ; cellID < numCells ; ++cellID )
  {
    EXPECT_EQ( field[ cellID ], cellID );
  }

  delete test_mesh;
  test_mesh = nullptr;
}

//------------------------------------------------------------------------------
template < typename ExecPolicy, int MeshType >
void check_for_all_cells_ij( )
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

  IndexType* icoords = test_mesh->createField< IndexType >( "i", CELL_CENTERED );
  IndexType* jcoords = test_mesh->createField< IndexType >( "j", CELL_CENTERED );

  for_all_cells< ExecPolicy, xargs::ij >( test_mesh,
    AXOM_LAMBDA( IndexType cellIdx, IndexType i, IndexType j )
    {
      icoords[ cellIdx ] = i;
      jcoords[ cellIdx ] = j;
    }
  );

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
}

//------------------------------------------------------------------------------
template < typename ExecPolicy, int MeshType >
void check_for_all_cells_ijk( )
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

  IndexType* icoords = test_mesh->createField< IndexType >( "i", CELL_CENTERED );
  IndexType* jcoords = test_mesh->createField< IndexType >( "j", CELL_CENTERED );
  IndexType* kcoords = test_mesh->createField< IndexType >( "k", CELL_CENTERED );

  for_all_cells< ExecPolicy, xargs::ijk >( test_mesh,
    AXOM_LAMBDA( IndexType cellIdx, IndexType i, IndexType j, IndexType k )
    {
      icoords[ cellIdx ] = i;
      jcoords[ cellIdx ] = j;
      kcoords[ cellIdx ] = k;
    }
  );

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
}

//------------------------------------------------------------------------------
template < typename ExecPolicy, int MeshType >
void check_for_all_cell_nodes( int dimension )
{
  SLIC_INFO( "dimension=" << dimension << ", policy="
            << policy_traits< ExecPolicy >::name() << ", mesh_type="
            << internal::mesh_type_name< MeshType >::name() );

  const IndexType Ni = 20;
  const IndexType Nj = (dimension >= 2) ? Ni : -1;
  const IndexType Nk = (dimension == 3) ? Ni : -1;

  const double lo[] = { -10, -10, -10 };
  const double hi[] = {  10,  10,  10 };
  UniformMesh uniform_mesh( lo, hi, Ni, Nj, Nk );

  Mesh* test_mesh = internal::create_mesh< MeshType >( uniform_mesh );
  EXPECT_TRUE( test_mesh != nullptr );

  const IndexType numCells = test_mesh->getNumberOfCells();
  IndexType* conn = test_mesh->createField< IndexType >( "conn", CELL_CENTERED, MAX_CELL_NODES );

  for_all_cells< ExecPolicy, xargs::nodeids >( test_mesh,
    AXOM_LAMBDA( IndexType cellID, const IndexType* nodes, IndexType N)
    {
      for ( int i = 0 ; i < N ; ++i )
      {
        conn[ cellID * MAX_CELL_NODES + i ] = nodes[ i ];
      } // END for all cell nodes
    } 
  );

  IndexType cellNodes[ MAX_CELL_NODES ];
  for ( IndexType cellID = 0 ; cellID < numCells ; ++cellID )
  {
    const IndexType N = test_mesh->getCellNodeIDs( cellID, cellNodes );
    for ( int i = 0 ; i < N ; ++i )
    {
      EXPECT_EQ( conn[ cellID * MAX_CELL_NODES + i ], cellNodes[ i ] );
    }
  }  // END for all cells

  /* clean up */
  delete test_mesh;
  test_mesh = nullptr;
}

//------------------------------------------------------------------------------
template < typename ExecPolicy, int MeshType >
void check_for_all_cell_faces( int dimension )
{
  SLIC_INFO( "dimension=" << dimension << ", policy="
            << policy_traits< ExecPolicy >::name() << ", mesh_type="
            << internal::mesh_type_name< MeshType >::name() );

  const IndexType Ni = 20;
  const IndexType Nj = (dimension >= 2) ? Ni : -1;
  const IndexType Nk = (dimension == 3) ? Ni : -1;

  const double lo[] = { -10, -10, -10 };
  const double hi[] = {  10,  10,  10 };
  UniformMesh uniform_mesh( lo, hi, Ni, Nj, Nk );

  Mesh* test_mesh = internal::create_mesh< MeshType >( uniform_mesh );
  EXPECT_TRUE( test_mesh != nullptr );

  const IndexType numCells        = test_mesh->getNumberOfCells();
  const IndexType numFacesPerCell = test_mesh->getNumberOfCellFaces();
  IndexType* cellFaces = test_mesh->createField< IndexType >( "cellFaces", 
                                                              CELL_CENTERED,
                                                              numFacesPerCell );

  for_all_cells< ExecPolicy, xargs::faceids >( test_mesh,
    AXOM_LAMBDA( IndexType cellID, const IndexType* faces, IndexType N)
    {
      EXPECT_EQ( N, numFacesPerCell );

      for ( int i = 0 ; i < N ; ++i )
      {
        cellFaces[ cellID * numFacesPerCell + i ] = faces[ i ];
      }
    }
  );

  IndexType faces[ 8 ];
  for ( IndexType cellID = 0 ; cellID < numCells ; ++cellID )
  {
    const IndexType N = test_mesh->getCellFaceIDs( cellID, faces );
    EXPECT_EQ( N, numFacesPerCell );

    for ( IndexType i=0 ; i < N ; ++i )
    {
      EXPECT_EQ( cellFaces[ cellID * numFacesPerCell + i ], faces[ i ] );
    }
  } 

  /* clean up */
  delete test_mesh;
  test_mesh = nullptr;
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
    check_for_all_cell_nodes< seq_exec, STRUCTURED_UNIFORM_MESH >(i);
    check_for_all_cell_nodes< seq_exec, STRUCTURED_CURVILINEAR_MESH >(i);
    check_for_all_cell_nodes< seq_exec, STRUCTURED_RECTILINEAR_MESH >(i);
    check_for_all_cell_nodes< seq_exec, UNSTRUCTURED_SINGLE_SHAPE_MESH >(i);
    check_for_all_cell_nodes< seq_exec, UNSTRUCTURED_MIXED_SHAPE_MESH >(i);

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && \
    defined(RAJA_ENABLE_OPENMP)

    using omp_exec = policy::parallel_cpu;
    check_for_all_cell_nodes< omp_exec, STRUCTURED_UNIFORM_MESH >(i);
    check_for_all_cell_nodes< omp_exec, STRUCTURED_CURVILINEAR_MESH >(i);
    check_for_all_cell_nodes< omp_exec, STRUCTURED_RECTILINEAR_MESH >(i);
    check_for_all_cell_nodes< omp_exec, UNSTRUCTURED_SINGLE_SHAPE_MESH >(i);
    check_for_all_cell_nodes< omp_exec, UNSTRUCTURED_MIXED_SHAPE_MESH >(i);

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && \
    defined(RAJA_ENABLE_CUDA)

    using cuda_exec = policy::parallel_gpu;
    check_for_all_cell_nodes< cuda_exec, STRUCTURED_UNIFORM_MESH >(i);
    check_for_all_cell_nodes< cuda_exec, STRUCTURED_CURVILINEAR_MESH >(i);
    check_for_all_cell_nodes< cuda_exec, STRUCTURED_RECTILINEAR_MESH >(i);
    check_for_all_cell_nodes< cuda_exec, UNSTRUCTURED_SINGLE_SHAPE_MESH >(i);
    check_for_all_cell_nodes< cuda_exec, UNSTRUCTURED_MIXED_SHAPE_MESH >(i);
#endif

  } // END for all dimensions
}

TEST( mint_execution_cell_traversals, for_all_cells_faceids )
{
  constexpr int NDIMS = 3;
  for ( int i=2 ; i <= NDIMS ; ++i )
  {

    using seq_exec = policy::serial;
    check_for_all_cell_faces< seq_exec, STRUCTURED_UNIFORM_MESH >(i);
    check_for_all_cell_faces< seq_exec, STRUCTURED_CURVILINEAR_MESH >(i);
    check_for_all_cell_faces< seq_exec, STRUCTURED_RECTILINEAR_MESH >(i);

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && \
    defined(RAJA_ENABLE_OPENMP)

    using omp_exec = policy::parallel_cpu;
    check_for_all_cell_faces< omp_exec, STRUCTURED_UNIFORM_MESH >(i);
    check_for_all_cell_faces< omp_exec, STRUCTURED_CURVILINEAR_MESH >(i);
    check_for_all_cell_faces< omp_exec, STRUCTURED_RECTILINEAR_MESH >(i);

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && \
    defined(RAJA_ENABLE_CUDA)

    using cuda_exec = policy::parallel_gpu;
    check_for_all_cell_faces< cuda_exec, STRUCTURED_UNIFORM_MESH >(i);
    check_for_all_cell_faces< cuda_exec, STRUCTURED_CURVILINEAR_MESH >(i);
    check_for_all_cell_faces< cuda_exec, STRUCTURED_RECTILINEAR_MESH >(i);

#endif

  } // END for all dimensions
}

//------------------------------------------------------------------------------
TEST( mint_execution_cell_traversals, for_all_cells_ij )
{
  using seq_exec = policy::serial;
  check_for_all_cells_ij< seq_exec, STRUCTURED_UNIFORM_MESH >();
  check_for_all_cells_ij< seq_exec, STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_cells_ij< seq_exec, STRUCTURED_RECTILINEAR_MESH >();

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && \
  defined(RAJA_ENABLE_OPENMP)

  using omp_exec = policy::parallel_cpu;
  check_for_all_cells_ij< omp_exec, STRUCTURED_UNIFORM_MESH >();
  check_for_all_cells_ij< omp_exec, STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_cells_ij< omp_exec, STRUCTURED_RECTILINEAR_MESH >();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && \
  defined(RAJA_ENABLE_CUDA)

  using cuda_exec = policy::parallel_gpu;
  check_for_all_cells_ij< cuda_exec, STRUCTURED_UNIFORM_MESH >();
  check_for_all_cells_ij< cuda_exec, STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_cells_ij< cuda_exec, STRUCTURED_RECTILINEAR_MESH >();

#endif

}

//------------------------------------------------------------------------------
TEST( mint_execution_cell_traversals, for_all_cells_ijk )
{
  using seq_exec = policy::serial;
  check_for_all_cells_ijk< seq_exec, STRUCTURED_UNIFORM_MESH >();
  check_for_all_cells_ijk< seq_exec, STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_cells_ijk< seq_exec, STRUCTURED_RECTILINEAR_MESH >();

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && \
  defined(RAJA_ENABLE_OPENMP)

  using omp_exec = policy::parallel_cpu;
  check_for_all_cells_ijk< omp_exec, STRUCTURED_UNIFORM_MESH >();
  check_for_all_cells_ijk< omp_exec, STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_cells_ijk< omp_exec, STRUCTURED_RECTILINEAR_MESH >();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && \
  defined(RAJA_ENABLE_CUDA)

  using cuda_exec = policy::parallel_gpu;
  check_for_all_cells_ijk< cuda_exec, STRUCTURED_UNIFORM_MESH >();
  check_for_all_cells_ijk< cuda_exec, STRUCTURED_CURVILINEAR_MESH >();
  check_for_all_cells_ijk< cuda_exec, STRUCTURED_RECTILINEAR_MESH >();

#endif

}

//------------------------------------------------------------------------------
TEST( mint_execution_cell_traversals, for_all_cells_index )
{
  constexpr int NDIMS = 3;
  for ( int i=1 ; i <= NDIMS ; ++i )
  {

    using seq_exec = policy::serial;
    check_for_all_cells_idx< seq_exec, STRUCTURED_UNIFORM_MESH >(i);
    check_for_all_cells_idx< seq_exec, STRUCTURED_CURVILINEAR_MESH >(i);
    check_for_all_cells_idx< seq_exec, STRUCTURED_RECTILINEAR_MESH >(i);
    check_for_all_cells_idx< seq_exec, UNSTRUCTURED_SINGLE_SHAPE_MESH >(i);
    check_for_all_cells_idx< seq_exec, UNSTRUCTURED_MIXED_SHAPE_MESH >(i);

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && \
    defined(RAJA_ENABLE_OPENMP)

    using omp_exec = policy::parallel_cpu;
    check_for_all_cells_idx< omp_exec, STRUCTURED_UNIFORM_MESH >(i);
    check_for_all_cells_idx< omp_exec, STRUCTURED_CURVILINEAR_MESH >(i);
    check_for_all_cells_idx< omp_exec, STRUCTURED_RECTILINEAR_MESH >(i);
    check_for_all_cells_idx< omp_exec, UNSTRUCTURED_SINGLE_SHAPE_MESH >(i);
    check_for_all_cells_idx< omp_exec, UNSTRUCTURED_MIXED_SHAPE_MESH >(i);

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && \
    defined(RAJA_ENABLE_CUDA)

    using cuda_exec = policy::parallel_gpu;
    check_for_all_cells_idx< cuda_exec, STRUCTURED_UNIFORM_MESH >(i);
    check_for_all_cells_idx< cuda_exec, STRUCTURED_CURVILINEAR_MESH >(i);
    check_for_all_cells_idx< cuda_exec, STRUCTURED_RECTILINEAR_MESH >(i);
    check_for_all_cells_idx< cuda_exec, UNSTRUCTURED_SINGLE_SHAPE_MESH >(i);
    check_for_all_cells_idx< cuda_exec, UNSTRUCTURED_MIXED_SHAPE_MESH >(i);

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

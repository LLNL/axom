// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/config.hpp"
#include "axom/core/execution/execution_space.hpp"
#include "axom/core/utilities/Utilities.hpp"
#include "axom/core/numerics/Matrix.hpp"

// Mint includes
#include "axom/mint/config.hpp"
#include "axom/mint/execution/interface.hpp"

// Slic includes
#include "axom/slic.hpp"

#include "mint_test_utilities.hpp"

// gtest includes
#include "gtest/gtest.h"

namespace axom
{
namespace mint
{
//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------
namespace
{
template <typename ExecPolicy, int MeshType, int Topology = SINGLE_SHAPE>
void check_for_all_cells_idx(int dimension)
{
  constexpr char* mesh_name = internal::mesh_type<MeshType, Topology>::name();
  SLIC_INFO("dimension=" << dimension
                         << ", policy=" << execution_space<ExecPolicy>::name()
                         << ", mesh_type=" << mesh_name);

  const IndexType Ni = 20;
  const IndexType Nj = (dimension >= 2) ? Ni : -1;
  const IndexType Nk = (dimension == 3) ? Ni : -1;

  const double lo[] = {-10, -10, -10};
  const double hi[] = {10, 10, 10};
  UniformMesh uniform_mesh(lo, hi, Ni, Nj, Nk);

  using MESH = typename internal::mesh_type<MeshType, Topology>::MeshType;
  MESH* test_mesh =
    dynamic_cast<MESH*>(internal::create_mesh<MeshType, Topology>(uniform_mesh));
  EXPECT_TRUE(test_mesh != nullptr);

  IndexType* field =
    test_mesh->template createField<IndexType>("c1", CELL_CENTERED);

  for_all_cells<ExecPolicy>(
    test_mesh,
    AXOM_LAMBDA(IndexType cellID) { field[cellID] = cellID; });

  const IndexType numCells = test_mesh->getNumberOfCells();
  for(IndexType cellID = 0; cellID < numCells; ++cellID)
  {
    EXPECT_EQ(field[cellID], cellID);
  }

  delete test_mesh;
  test_mesh = nullptr;
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, int MeshType>
void check_for_all_cells_ij()
{
  SLIC_INFO("policy=" << execution_space<ExecPolicy>::name() << ", mesh_type="
                      << internal::mesh_type<MeshType>::name());

  constexpr IndexType N = 20;
  const double lo[] = {-10, -10};
  const double hi[] = {10, 10};
  UniformMesh uniform_mesh(lo, hi, N, N);

  // STEP 0: create the test mesh
  using MESH = typename internal::mesh_type<MeshType>::MeshType;
  MESH* test_mesh =
    dynamic_cast<MESH*>(internal::create_mesh<MeshType>(uniform_mesh));
  EXPECT_TRUE(test_mesh != nullptr);

  IndexType* icoords =
    test_mesh->template createField<IndexType>("i", CELL_CENTERED);
  IndexType* jcoords =
    test_mesh->template createField<IndexType>("j", CELL_CENTERED);

  for_all_cells<ExecPolicy, xargs::ij>(
    test_mesh,
    AXOM_LAMBDA(IndexType cellIdx, IndexType i, IndexType j) {
      icoords[cellIdx] = i;
      jcoords[cellIdx] = j;
    });

  IndexType icell = 0;
  for(IndexType j = 0; j < (N - 1); ++j)
  {
    for(IndexType i = 0; i < (N - 1); ++i)
    {
      EXPECT_EQ(icoords[icell], i);
      EXPECT_EQ(jcoords[icell], j);
      ++icell;
    }  // END for all i
  }    // END for all j

  delete test_mesh;
  test_mesh = nullptr;
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, int MeshType>
void check_for_all_cells_ijk()
{
  SLIC_INFO("policy=" << execution_space<ExecPolicy>::name() << ", mesh_type="
                      << internal::mesh_type<MeshType>::name());

  constexpr IndexType N = 20;
  const double lo[] = {-10, -10, -10};
  const double hi[] = {10, 10, 10};
  UniformMesh uniform_mesh(lo, hi, N, N, N);

  // STEP 0: create the test mesh
  using MESH = typename internal::mesh_type<MeshType>::MeshType;
  MESH* test_mesh =
    dynamic_cast<MESH*>(internal::create_mesh<MeshType>(uniform_mesh));
  EXPECT_TRUE(test_mesh != nullptr);

  IndexType* icoords =
    test_mesh->template createField<IndexType>("i", CELL_CENTERED);
  IndexType* jcoords =
    test_mesh->template createField<IndexType>("j", CELL_CENTERED);
  IndexType* kcoords =
    test_mesh->template createField<IndexType>("k", CELL_CENTERED);

  for_all_cells<ExecPolicy, xargs::ijk>(
    test_mesh,
    AXOM_LAMBDA(IndexType cellIdx, IndexType i, IndexType j, IndexType k) {
      icoords[cellIdx] = i;
      jcoords[cellIdx] = j;
      kcoords[cellIdx] = k;
    });

  IndexType icell = 0;
  for(IndexType k = 0; k < (N - 1); ++k)
  {
    for(IndexType j = 0; j < (N - 1); ++j)
    {
      for(IndexType i = 0; i < (N - 1); ++i)
      {
        EXPECT_EQ(icoords[icell], i);
        EXPECT_EQ(jcoords[icell], j);
        EXPECT_EQ(kcoords[icell], k);
        ++icell;
      }  // END for all i
    }    // END for all j
  }      // END for all k

  delete test_mesh;
  test_mesh = nullptr;
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, int MeshType, int Topology = SINGLE_SHAPE>
void check_for_all_cell_nodes(int dimension)
{
  constexpr char* mesh_name = internal::mesh_type<MeshType, Topology>::name();
  SLIC_INFO("dimension=" << dimension
                         << ", policy=" << execution_space<ExecPolicy>::name()
                         << ", mesh_type=" << mesh_name);

  const IndexType Ni = 20;
  const IndexType Nj = (dimension >= 2) ? Ni : -1;
  const IndexType Nk = (dimension == 3) ? Ni : -1;

  const double lo[] = {-10, -10, -10};
  const double hi[] = {10, 10, 10};
  UniformMesh uniform_mesh(lo, hi, Ni, Nj, Nk);

  using MESH = typename internal::mesh_type<MeshType, Topology>::MeshType;
  MESH* test_mesh =
    dynamic_cast<MESH*>(internal::create_mesh<MeshType, Topology>(uniform_mesh));
  EXPECT_TRUE(test_mesh != nullptr);

  const IndexType numCells = test_mesh->getNumberOfCells();
  IndexType* conn = test_mesh->template createField<IndexType>("conn",
                                                               CELL_CENTERED,
                                                               MAX_CELL_NODES);

  for_all_cells<ExecPolicy, xargs::nodeids>(
    test_mesh,
    AXOM_LAMBDA(IndexType cellID, const IndexType* nodes, IndexType N) {
      for(int i = 0; i < N; ++i)
      {
        conn[cellID * MAX_CELL_NODES + i] = nodes[i];
      }  // END for all cell nodes
    });

  IndexType cellNodes[MAX_CELL_NODES];
  for(IndexType cellID = 0; cellID < numCells; ++cellID)
  {
    const IndexType N = test_mesh->getCellNodeIDs(cellID, cellNodes);
    for(int i = 0; i < N; ++i)
    {
      EXPECT_EQ(conn[cellID * MAX_CELL_NODES + i], cellNodes[i]);
    }
  }  // END for all cells

  /* clean up */
  delete test_mesh;
  test_mesh = nullptr;
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, int MeshType, int Topology = SINGLE_SHAPE>
void check_for_all_cell_coords(int dimension)
{
  constexpr char* mesh_name = internal::mesh_type<MeshType, Topology>::name();
  SLIC_INFO("dimension=" << dimension
                         << ", policy=" << execution_space<ExecPolicy>::name()
                         << ", mesh_type=" << mesh_name);

  const IndexType Ni = 20;
  const IndexType Nj = (dimension >= 2) ? Ni : -1;
  const IndexType Nk = (dimension == 3) ? Ni : -1;

  const double lo[] = {-10, -10, -10};
  const double hi[] = {10, 10, 10};
  UniformMesh uniform_mesh(lo, hi, Ni, Nj, Nk);

  using MESH = typename internal::mesh_type<MeshType, Topology>::MeshType;
  MESH* test_mesh =
    dynamic_cast<MESH*>(internal::create_mesh<MeshType, Topology>(uniform_mesh));
  EXPECT_TRUE(test_mesh != nullptr);

  const IndexType numCells = test_mesh->getNumberOfCells();
  IndexType* conn = test_mesh->template createField<IndexType>("conn",
                                                               CELL_CENTERED,
                                                               MAX_CELL_NODES);
  double* coords =
    test_mesh->template createField<double>("coords",
                                            CELL_CENTERED,
                                            dimension * MAX_CELL_NODES);

  for_all_cells<ExecPolicy, xargs::coords>(
    test_mesh,
    AXOM_LAMBDA(IndexType cellID,
                const numerics::Matrix<double>& coordsMatrix,
                const IndexType* nodes) {
      const IndexType numNodes = coordsMatrix.getNumColumns();
      for(int i = 0; i < numNodes; ++i)
      {
        conn[cellID * MAX_CELL_NODES + i] = nodes[i];

        for(int dim = 0; dim < dimension; ++dim)
        {
          coords[cellID * dimension * MAX_CELL_NODES + i * dimension + dim] =
            coordsMatrix(dim, i);
        }
      }  // END for all cell nodes
    });

  double nodeCoords[3];
  IndexType cellNodes[MAX_CELL_NODES];
  for(IndexType cellID = 0; cellID < numCells; ++cellID)
  {
    const IndexType numNodes = test_mesh->getCellNodeIDs(cellID, cellNodes);
    for(int i = 0; i < numNodes; ++i)
    {
      EXPECT_EQ(conn[cellID * MAX_CELL_NODES + i], cellNodes[i]);

      for(int dim = 0; dim < dimension; ++dim)
      {
        test_mesh->getNode(cellNodes[i], nodeCoords);
        EXPECT_NEAR(
          coords[cellID * dimension * MAX_CELL_NODES + i * dimension + dim],
          nodeCoords[dim],
          1e-8);
      }
    }
  }  // END for all cells

  /* clean up */
  delete test_mesh;
  test_mesh = nullptr;
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, int MeshType, int Topology = SINGLE_SHAPE>
void check_for_all_cell_faces(int dimension)
{
  constexpr char* mesh_name = internal::mesh_type<MeshType, Topology>::name();
  SLIC_INFO("dimension=" << dimension
                         << ", policy=" << execution_space<ExecPolicy>::name()
                         << ", mesh_type=" << mesh_name);

  const IndexType Ni = 20;
  const IndexType Nj = (dimension >= 2) ? Ni : -1;
  const IndexType Nk = (dimension == 3) ? Ni : -1;

  const double lo[] = {-10, -10, -10};
  const double hi[] = {10, 10, 10};
  UniformMesh uniform_mesh(lo, hi, Ni, Nj, Nk);

  using MESH = typename internal::mesh_type<MeshType, Topology>::MeshType;
  MESH* test_mesh =
    dynamic_cast<MESH*>(internal::create_mesh<MeshType, Topology>(uniform_mesh));
  EXPECT_TRUE(test_mesh != nullptr);

  const IndexType numCells = test_mesh->getNumberOfCells();
  IndexType* cellFaces =
    test_mesh->template createField<IndexType>("cellFaces",
                                               CELL_CENTERED,
                                               MAX_CELL_FACES);

  for_all_cells<ExecPolicy, xargs::faceids>(
    test_mesh,
    AXOM_LAMBDA(IndexType cellID, const IndexType* faces, IndexType N) {
      for(int i = 0; i < N; ++i)
      {
        cellFaces[cellID * MAX_CELL_FACES + i] = faces[i];
      }
    });

  IndexType faces[MAX_CELL_FACES];
  for(IndexType cellID = 0; cellID < numCells; ++cellID)
  {
    const IndexType N = test_mesh->getCellFaceIDs(cellID, faces);

    for(IndexType i = 0; i < N; ++i)
    {
      EXPECT_EQ(cellFaces[cellID * MAX_CELL_FACES + i], faces[i]);
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

AXOM_CUDA_TEST(mint_execution_cell_traversals, for_all_cells_nodeids)
{
  constexpr int NDIMS = 3;
  for(int i = 1; i <= NDIMS; ++i)
  {
    using seq_exec = axom::SEQ_EXEC;
    check_for_all_cell_nodes<seq_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_cell_nodes<seq_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_cell_nodes<seq_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_cell_nodes<seq_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_cell_nodes<seq_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && \
  defined(RAJA_ENABLE_OPENMP)

    using omp_exec = axom::OMP_EXEC;
    check_for_all_cell_nodes<omp_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_cell_nodes<omp_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_cell_nodes<omp_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_cell_nodes<omp_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_cell_nodes<omp_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && \
  defined(RAJA_ENABLE_CUDA) && defined(AXOM_USE_UMPIRE)

    using cuda_exec = axom::CUDA_EXEC<512>;

    // Use unified memory
    const int exec_space_id = axom::getUmpireResourceAllocatorID(
      umpire::resource::MemoryResourceType::Unified);
    const int prev_allocator = axom::getDefaultAllocatorID();
    axom::setDefaultAllocator(exec_space_id);

    check_for_all_cell_nodes<cuda_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_cell_nodes<cuda_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_cell_nodes<cuda_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_cell_nodes<cuda_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_cell_nodes<cuda_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

    setDefaultAllocator(prev_allocator);
#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_HIP) && \
  defined(RAJA_ENABLE_HIP) && defined(AXOM_USE_UMPIRE)

    using hip_exec = axom::HIP_EXEC<512>;

    // Use unified memory
    const int exec_space_id = axom::getUmpireResourceAllocatorID(
      umpire::resource::MemoryResourceType::Unified);
    const int prev_allocator = axom::getDefaultAllocatorID();
    axom::setDefaultAllocator(exec_space_id);

    check_for_all_cell_nodes<hip_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_cell_nodes<hip_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_cell_nodes<hip_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_cell_nodes<hip_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_cell_nodes<hip_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

    setDefaultAllocator(prev_allocator);
#endif

  }  // END for all dimensions
}

AXOM_CUDA_TEST(mint_execution_cell_traversals, for_all_cells_coords)
{
  constexpr int NDIMS = 3;
  for(int i = 1; i <= NDIMS; ++i)
  {
    using seq_exec = axom::SEQ_EXEC;
    check_for_all_cell_coords<seq_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_cell_coords<seq_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_cell_coords<seq_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_cell_coords<seq_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_cell_coords<seq_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && \
  defined(RAJA_ENABLE_OPENMP)

    using omp_exec = axom::OMP_EXEC;
    check_for_all_cell_coords<omp_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_cell_coords<omp_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_cell_coords<omp_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_cell_coords<omp_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_cell_coords<omp_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && \
  defined(RAJA_ENABLE_CUDA) && defined(AXOM_USE_UMPIRE)

    using cuda_exec = axom::CUDA_EXEC<512>;

    // Use unified memory
    const int exec_space_id = axom::getUmpireResourceAllocatorID(
      umpire::resource::MemoryResourceType::Unified);
    const int prev_allocator = axom::getDefaultAllocatorID();
    axom::setDefaultAllocator(exec_space_id);

    check_for_all_cell_coords<cuda_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_cell_coords<cuda_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_cell_coords<cuda_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_cell_coords<cuda_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_cell_coords<cuda_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

    setDefaultAllocator(prev_allocator);
#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_HIP) && \
  defined(RAJA_ENABLE_HIP) && defined(AXOM_USE_UMPIRE)

    using hip_exec = axom::HIP_EXEC<512>;

    // Use unified memory
    const int exec_space_id = axom::getUmpireResourceAllocatorID(
      umpire::resource::MemoryResourceType::Unified);
    const int prev_allocator = axom::getDefaultAllocatorID();
    axom::setDefaultAllocator(exec_space_id);

    check_for_all_cell_coords<hip_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_cell_coords<hip_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_cell_coords<hip_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_cell_coords<hip_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_cell_coords<hip_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

    setDefaultAllocator(prev_allocator);
#endif

  }  // END for all dimensions
}

AXOM_CUDA_TEST(mint_execution_cell_traversals, for_all_cells_faceids)
{
  constexpr int NDIMS = 3;
  for(int i = 2; i <= NDIMS; ++i)
  {
    using seq_exec = axom::SEQ_EXEC;
    check_for_all_cell_faces<seq_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_cell_faces<seq_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_cell_faces<seq_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_cell_faces<seq_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_cell_faces<seq_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && \
  defined(RAJA_ENABLE_OPENMP)

    using omp_exec = axom::OMP_EXEC;
    check_for_all_cell_faces<omp_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_cell_faces<omp_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_cell_faces<omp_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_cell_faces<seq_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_cell_faces<seq_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && \
  defined(RAJA_ENABLE_CUDA) && defined(AXOM_USE_UMPIRE)

    using cuda_exec = axom::CUDA_EXEC<512>;

    // Use unified memory
    const int exec_space_id = axom::getUmpireResourceAllocatorID(
      umpire::resource::MemoryResourceType::Unified);
    const int prev_allocator = axom::getDefaultAllocatorID();
    axom::setDefaultAllocator(exec_space_id);

    check_for_all_cell_faces<cuda_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_cell_faces<cuda_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_cell_faces<cuda_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_cell_faces<seq_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_cell_faces<seq_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

    setDefaultAllocator(prev_allocator);
#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_HIP) && \
  defined(RAJA_ENABLE_HIP) && defined(AXOM_USE_UMPIRE)

    using hip_exec = axom::HIP_EXEC<512>;

    // Use unified memory
    const int exec_space_id = axom::getUmpireResourceAllocatorID(
      umpire::resource::MemoryResourceType::Unified);
    const int prev_allocator = axom::getDefaultAllocatorID();
    axom::setDefaultAllocator(exec_space_id);

    check_for_all_cell_faces<hip_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_cell_faces<hip_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_cell_faces<hip_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_cell_faces<hip_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_cell_faces<hip_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

    setDefaultAllocator(prev_allocator);
#endif

  }  // END for all dimensions
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST(mint_execution_cell_traversals, for_all_cells_ij)
{
  using seq_exec = axom::SEQ_EXEC;
  check_for_all_cells_ij<seq_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_cells_ij<seq_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_cells_ij<seq_exec, STRUCTURED_RECTILINEAR_MESH>();

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && \
  defined(RAJA_ENABLE_OPENMP)

  using omp_exec = axom::OMP_EXEC;
  check_for_all_cells_ij<omp_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_cells_ij<omp_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_cells_ij<omp_exec, STRUCTURED_RECTILINEAR_MESH>();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && \
  defined(RAJA_ENABLE_CUDA) && defined(AXOM_USE_UMPIRE)

  using cuda_exec = axom::CUDA_EXEC<512>;

  // Use unified memory
  const int exec_space_id = axom::getUmpireResourceAllocatorID(
    umpire::resource::MemoryResourceType::Unified);
  const int prev_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(exec_space_id);

  check_for_all_cells_ij<cuda_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_cells_ij<cuda_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_cells_ij<cuda_exec, STRUCTURED_RECTILINEAR_MESH>();

  setDefaultAllocator(prev_allocator);
#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_HIP) && \
  defined(RAJA_ENABLE_HIP) && defined(AXOM_USE_UMPIRE)

  using hip_exec = axom::HIP_EXEC<512>;

  // Use unified memory
  const int exec_space_id = axom::getUmpireResourceAllocatorID(
    umpire::resource::MemoryResourceType::Unified);
  const int prev_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(exec_space_id);

  check_for_all_cells_ij<hip_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_cells_ij<hip_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_cells_ij<hip_exec, STRUCTURED_RECTILINEAR_MESH>();

  setDefaultAllocator(prev_allocator);
#endif
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST(mint_execution_cell_traversals, for_all_cells_ijk)
{
  using seq_exec = axom::SEQ_EXEC;
  check_for_all_cells_ijk<seq_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_cells_ijk<seq_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_cells_ijk<seq_exec, STRUCTURED_RECTILINEAR_MESH>();

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && \
  defined(RAJA_ENABLE_OPENMP)

  using omp_exec = axom::OMP_EXEC;
  check_for_all_cells_ijk<omp_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_cells_ijk<omp_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_cells_ijk<omp_exec, STRUCTURED_RECTILINEAR_MESH>();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && \
  defined(RAJA_ENABLE_CUDA) && defined(AXOM_USE_UMPIRE)

  using cuda_exec = axom::CUDA_EXEC<512>;

  // Use unified memory
  const int exec_space_id = axom::getUmpireResourceAllocatorID(
    umpire::resource::MemoryResourceType::Unified);
  const int prev_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(exec_space_id);

  check_for_all_cells_ijk<cuda_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_cells_ijk<cuda_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_cells_ijk<cuda_exec, STRUCTURED_RECTILINEAR_MESH>();

  setDefaultAllocator(prev_allocator);
#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_HIP) && \
  defined(RAJA_ENABLE_HIP) && defined(AXOM_USE_UMPIRE)

  using hip_exec = axom::HIP_EXEC<512>;

  // Use unified memory
  const int exec_space_id = axom::getUmpireResourceAllocatorID(
    umpire::resource::MemoryResourceType::Unified);
  const int prev_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(exec_space_id);

  check_for_all_cells_ijk<hip_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_cells_ijk<hip_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_cells_ijk<hip_exec, STRUCTURED_RECTILINEAR_MESH>();

  setDefaultAllocator(prev_allocator);
#endif
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST(mint_execution_cell_traversals, for_all_cells_index)
{
  constexpr int NDIMS = 3;
  for(int i = 1; i <= NDIMS; ++i)
  {
    using seq_exec = axom::SEQ_EXEC;
    check_for_all_cells_idx<seq_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_cells_idx<seq_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_cells_idx<seq_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_cells_idx<seq_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_cells_idx<seq_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && \
  defined(RAJA_ENABLE_OPENMP)

    using omp_exec = axom::OMP_EXEC;
    check_for_all_cells_idx<omp_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_cells_idx<omp_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_cells_idx<omp_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_cells_idx<omp_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_cells_idx<omp_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && \
  defined(RAJA_ENABLE_CUDA) && defined(AXOM_USE_UMPIRE)

    using cuda_exec = axom::CUDA_EXEC<512>;

    // Use unified memory
    const int exec_space_id = axom::getUmpireResourceAllocatorID(
      umpire::resource::MemoryResourceType::Unified);
    const int prev_allocator = axom::getDefaultAllocatorID();
    axom::setDefaultAllocator(exec_space_id);

    check_for_all_cells_idx<cuda_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_cells_idx<cuda_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_cells_idx<cuda_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_cells_idx<cuda_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_cells_idx<cuda_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

    setDefaultAllocator(prev_allocator);
#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_HIP) && \
  defined(RAJA_ENABLE_HIP) && defined(AXOM_USE_UMPIRE)

    using hip_exec = axom::HIP_EXEC<512>;

    // Use unified memory
    const int exec_space_id = axom::getUmpireResourceAllocatorID(
      umpire::resource::MemoryResourceType::Unified);

    const int prev_allocator = axom::getDefaultAllocatorID();
    axom::setDefaultAllocator(exec_space_id);

    check_for_all_cells_idx<hip_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_cells_idx<hip_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_cells_idx<hip_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_cells_idx<hip_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_cells_idx<hip_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

    setDefaultAllocator(prev_allocator);
#endif

  }  // END for all dimensions
}

} /* namespace mint */
} /* namespace axom */

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}

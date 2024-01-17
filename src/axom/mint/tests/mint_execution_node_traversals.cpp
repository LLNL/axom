// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/config.hpp"                          // compile-time definitions
#include "axom/core/execution/execution_space.hpp"  // for execution_space traits
#include "axom/core/memory_management.hpp"          // for alloc() /free()

// Mint includes
#include "axom/mint/config.hpp"               // mint compile-time definitions
#include "axom/mint/execution/interface.hpp"  // for_all()

// Slic includes
#include "axom/slic.hpp"  // for SLIC macros

#include "mint_test_utilities.hpp"

// gtest includes
#include "gtest/gtest.h"  // for gtest

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
void check_for_all_nodes_idx(int dimension)
{
  constexpr char* mesh_name = internal::mesh_type<MeshType, Topology>::name();
  SLIC_INFO("dimension=" << dimension
                         << ", policy=" << execution_space<ExecPolicy>::name()
                         << ", mesh_type=" << mesh_name);

  constexpr int MAGIC_VAL = 42;

  const IndexType Ni = 20;
  const IndexType Nj = (dimension >= 2) ? 20 : -1;
  const IndexType Nk = (dimension == 3) ? 20 : -1;

  const double lo[] = {-10, -10, -10};
  const double hi[] = {10, 10, 10};
  UniformMesh uniform_mesh(lo, hi, Ni, Nj, Nk);

  using MESH = typename internal::mesh_type<MeshType, Topology>::MeshType;
  MESH* test_mesh =
    dynamic_cast<MESH*>(internal::create_mesh<MeshType, Topology>(uniform_mesh));
  EXPECT_TRUE(test_mesh != nullptr);

  int* field = test_mesh->template createField<int>("n1", NODE_CENTERED);

  for_all_nodes<ExecPolicy>(
    test_mesh,
    AXOM_LAMBDA(IndexType nodeIdx) { field[nodeIdx] = MAGIC_VAL; });

  const IndexType numNodes = test_mesh->getNumberOfNodes();
  for(IndexType inode = 0; inode < numNodes; ++inode)
  {
    EXPECT_EQ(field[inode], MAGIC_VAL);
  }

  delete test_mesh;
  test_mesh = nullptr;
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, int MeshType>
void check_for_all_nodes_ij()
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

  const IndexType numNodes = test_mesh->getNumberOfNodes();

  IndexType* icoords = axom::allocate<IndexType>(numNodes);
  IndexType* jcoords = axom::allocate<IndexType>(numNodes);

  for_all_nodes<ExecPolicy, xargs::ij>(
    test_mesh,
    AXOM_LAMBDA(IndexType nodeIdx, IndexType i, IndexType j) {
      icoords[nodeIdx] = i;
      jcoords[nodeIdx] = j;
    });

  IndexType inode = 0;
  for(IndexType j = 0; j < N; ++j)
  {
    for(IndexType i = 0; i < N; ++i)
    {
      EXPECT_EQ(icoords[inode], i);
      EXPECT_EQ(jcoords[inode], j);
      ++inode;

    }  // END for all i
  }    // END for all j

  delete test_mesh;
  test_mesh = nullptr;

  axom::deallocate(icoords);
  axom::deallocate(jcoords);
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, int MeshType>
void check_for_all_nodes_ijk()
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

  const IndexType numNodes = test_mesh->getNumberOfNodes();

  IndexType* icoords = axom::allocate<IndexType>(numNodes);
  IndexType* jcoords = axom::allocate<IndexType>(numNodes);
  IndexType* kcoords = axom::allocate<IndexType>(numNodes);

  for_all_nodes<ExecPolicy, xargs::ijk>(
    test_mesh,
    AXOM_LAMBDA(IndexType nodeIdx, IndexType i, IndexType j, IndexType k) {
      icoords[nodeIdx] = i;
      jcoords[nodeIdx] = j;
      kcoords[nodeIdx] = k;
    });

  IndexType inode = 0;
  for(IndexType k = 0; k < N; ++k)
  {
    for(IndexType j = 0; j < N; ++j)
    {
      for(IndexType i = 0; i < N; ++i)
      {
        EXPECT_EQ(icoords[inode], i);
        EXPECT_EQ(jcoords[inode], j);
        EXPECT_EQ(kcoords[inode], k);
        ++inode;

      }  // END for all i
    }    // END for all j
  }      // END for all k

  delete test_mesh;
  test_mesh = nullptr;

  axom::deallocate(icoords);
  axom::deallocate(jcoords);
  axom::deallocate(kcoords);
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, int MeshType, int Topology = SINGLE_SHAPE>
void check_for_all_nodes_xyz()
{
  constexpr char* mesh_name = internal::mesh_type<MeshType, Topology>::name();
  SLIC_INFO("policy=" << execution_space<ExecPolicy>::name()
                      << ", mesh_type=" << mesh_name);

  constexpr IndexType N = 20;
  const double lo[] = {-10, -10, -10};
  const double hi[] = {10, 10, 10};
  UniformMesh uniform_mesh(lo, hi, N, N, N);

  // STEP 0: create the test mesh
  using MESH = typename internal::mesh_type<MeshType, Topology>::MeshType;
  MESH* test_mesh =
    dynamic_cast<MESH*>(internal::create_mesh<MeshType, Topology>(uniform_mesh));
  EXPECT_TRUE(test_mesh != nullptr);

  // STEP 1: generate test coordinate arrays
  const IndexType numNodes = test_mesh->getNumberOfNodes();
  double* x = axom::allocate<double>(numNodes);
  double* y = axom::allocate<double>(numNodes);
  double* z = axom::allocate<double>(numNodes);
  for_all_nodes<ExecPolicy, xargs::xyz>(
    test_mesh,
    AXOM_LAMBDA(IndexType idx, double xx, double yy, double zz) {
      x[idx] = xx;
      y[idx] = yy;
      z[idx] = zz;
    });

  // STEP 2:check coordinate arrays
  for(int inode = 0; inode < numNodes; ++inode)
  {
    double node[3];
    test_mesh->getNode(inode, node);
    EXPECT_DOUBLE_EQ(x[inode], node[X_COORDINATE]);
    EXPECT_DOUBLE_EQ(y[inode], node[Y_COORDINATE]);
    EXPECT_DOUBLE_EQ(z[inode], node[Z_COORDINATE]);
  }  // END for all nodes

  // STEP 3: clean up
  delete test_mesh;
  test_mesh = nullptr;

  axom::deallocate(x);
  axom::deallocate(y);
  axom::deallocate(z);
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, int MeshType, int Topology = SINGLE_SHAPE>
void check_for_all_nodes_xy()
{
  constexpr char* mesh_name = internal::mesh_type<MeshType, Topology>::name();
  SLIC_INFO("policy=" << execution_space<ExecPolicy>::name()
                      << ", mesh_type=" << mesh_name);

  constexpr IndexType N = 20;
  const double lo[] = {-10, -10};
  const double hi[] = {10, 10};
  UniformMesh uniform_mesh(lo, hi, N, N);

  // STEP 0: create the test mesh
  using MESH = typename internal::mesh_type<MeshType, Topology>::MeshType;
  MESH* test_mesh =
    dynamic_cast<MESH*>(internal::create_mesh<MeshType, Topology>(uniform_mesh));
  EXPECT_TRUE(test_mesh != nullptr);

  // STEP 1: generate test coordinate arrays
  const IndexType numNodes = test_mesh->getNumberOfNodes();
  double* x = axom::allocate<double>(numNodes);
  double* y = axom::allocate<double>(numNodes);
  for_all_nodes<ExecPolicy, xargs::xy>(
    test_mesh,
    AXOM_LAMBDA(IndexType idx, double xx, double yy) {
      x[idx] = xx;
      y[idx] = yy;
    });

  // STEP 2:check coordinate arrays
  for(int inode = 0; inode < numNodes; ++inode)
  {
    double node[2];
    test_mesh->getNode(inode, node);
    EXPECT_DOUBLE_EQ(x[inode], node[X_COORDINATE]);
    EXPECT_DOUBLE_EQ(y[inode], node[Y_COORDINATE]);
  }  // END for all nodes

  // STEP 3: clean up
  delete test_mesh;
  test_mesh = nullptr;

  axom::deallocate(x);
  axom::deallocate(y);
}

//------------------------------------------------------------------------------
template <typename ExecPolicy, int MeshType, int Topology = SINGLE_SHAPE>
void check_for_all_nodes_x()
{
  constexpr char* mesh_name = internal::mesh_type<MeshType, Topology>::name();
  SLIC_INFO("policy=" << execution_space<ExecPolicy>::name()
                      << ", mesh_type=" << mesh_name);

  constexpr IndexType Ni = 20;
  const double lo[] = {-10};
  const double hi[] = {10};
  UniformMesh uniform_mesh(lo, hi, Ni);

  // STEP 0: create the test mesh
  using MESH = typename internal::mesh_type<MeshType, Topology>::MeshType;
  MESH* test_mesh =
    dynamic_cast<MESH*>(internal::create_mesh<MeshType, Topology>(uniform_mesh));
  EXPECT_TRUE(test_mesh != nullptr);

  if(MeshType != PARTICLE_MESH)
  {
    EXPECT_EQ(test_mesh->getNumberOfCells(), uniform_mesh.getNumberOfCells());
  }

  // STEP 1: generate test coordinate arrays
  const IndexType numNodes = uniform_mesh.getNumberOfNodes();
  double* x = axom::allocate<double>(numNodes);
  for_all_nodes<ExecPolicy, xargs::x>(
    test_mesh,
    AXOM_LAMBDA(IndexType idx, double xx) { x[idx] = xx; });

  // STEP 2:check coordinate arrays
  for(int inode = 0; inode < numNodes; ++inode)
  {
    double node[1];
    uniform_mesh.getNode(inode, node);
    EXPECT_DOUBLE_EQ(x[inode], node[X_COORDINATE]);
  }  // END for all nodes

  // STEP 3: clean up
  delete test_mesh;
  test_mesh = nullptr;
  axom::deallocate(x);
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------

AXOM_CUDA_TEST(mint_execution_node_traversals, for_all_nodes_xyz)
{
  using seq_exec = axom::SEQ_EXEC;
  check_for_all_nodes_xyz<seq_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_xyz<seq_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_xyz<seq_exec, STRUCTURED_RECTILINEAR_MESH>();
  check_for_all_nodes_xyz<seq_exec, PARTICLE_MESH>();
  check_for_all_nodes_xyz<seq_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>();
  check_for_all_nodes_xyz<seq_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>();

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && \
  defined(RAJA_ENABLE_OPENMP)

  using openmp_exec = axom::OMP_EXEC;
  check_for_all_nodes_xyz<openmp_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_xyz<openmp_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_xyz<openmp_exec, STRUCTURED_RECTILINEAR_MESH>();
  check_for_all_nodes_xyz<openmp_exec, PARTICLE_MESH>();
  check_for_all_nodes_xyz<openmp_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>();
  check_for_all_nodes_xyz<openmp_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && \
  defined(RAJA_ENABLE_CUDA) && defined(AXOM_USE_UMPIRE)

  using cuda_exec = axom::CUDA_EXEC<512>;

  // Use unified memory
  const int exec_space_id = axom::getUmpireResourceAllocatorID(
    umpire::resource::MemoryResourceType::Unified);
  const int prev_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(exec_space_id);

  check_for_all_nodes_xyz<cuda_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_xyz<cuda_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_xyz<cuda_exec, STRUCTURED_RECTILINEAR_MESH>();
  check_for_all_nodes_xyz<cuda_exec, PARTICLE_MESH>();
  check_for_all_nodes_xyz<cuda_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>();
  check_for_all_nodes_xyz<cuda_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>();

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

  check_for_all_nodes_xyz<hip_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_xyz<hip_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_xyz<hip_exec, STRUCTURED_RECTILINEAR_MESH>();
  check_for_all_nodes_xyz<hip_exec, PARTICLE_MESH>();
  check_for_all_nodes_xyz<hip_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>();
  check_for_all_nodes_xyz<hip_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>();

  setDefaultAllocator(prev_allocator);
#endif
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST(mint_execution_node_traversals, for_all_nodes_xy)
{
  using seq_exec = axom::SEQ_EXEC;
  check_for_all_nodes_xy<seq_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_xy<seq_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_xy<seq_exec, STRUCTURED_RECTILINEAR_MESH>();
  check_for_all_nodes_xy<seq_exec, PARTICLE_MESH>();
  check_for_all_nodes_xy<seq_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>();
  check_for_all_nodes_xy<seq_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>();

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && \
  defined(RAJA_ENABLE_OPENMP)

  using openmp_exec = axom::OMP_EXEC;
  check_for_all_nodes_xy<openmp_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_xy<openmp_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_xy<openmp_exec, STRUCTURED_RECTILINEAR_MESH>();
  check_for_all_nodes_xy<openmp_exec, PARTICLE_MESH>();
  check_for_all_nodes_xy<openmp_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>();
  check_for_all_nodes_xy<openmp_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && \
  defined(RAJA_ENABLE_CUDA) && defined(AXOM_USE_UMPIRE)

  using cuda_exec = axom::CUDA_EXEC<512>;

  // Use unified memory
  const int exec_space_id = axom::getUmpireResourceAllocatorID(
    umpire::resource::MemoryResourceType::Unified);
  const int prev_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(exec_space_id);

  check_for_all_nodes_xy<cuda_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_xy<cuda_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_xy<cuda_exec, STRUCTURED_RECTILINEAR_MESH>();
  check_for_all_nodes_xy<cuda_exec, PARTICLE_MESH>();
  check_for_all_nodes_xy<cuda_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>();
  check_for_all_nodes_xy<cuda_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>();

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

  check_for_all_nodes_xy<hip_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_xy<hip_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_xy<hip_exec, STRUCTURED_RECTILINEAR_MESH>();
  check_for_all_nodes_xy<hip_exec, PARTICLE_MESH>();
  check_for_all_nodes_xy<hip_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>();
  check_for_all_nodes_xy<hip_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>();

  setDefaultAllocator(prev_allocator);
#endif
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST(mint_execution_node_traversals, for_all_nodes_x)
{
  using seq_exec = axom::SEQ_EXEC;
  check_for_all_nodes_x<seq_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_x<seq_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_x<seq_exec, STRUCTURED_RECTILINEAR_MESH>();
  check_for_all_nodes_x<seq_exec, PARTICLE_MESH>();
  check_for_all_nodes_x<seq_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>();
  check_for_all_nodes_x<seq_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>();

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && \
  defined(RAJA_ENABLE_OPENMP)

  using openmp_exec = axom::OMP_EXEC;
  check_for_all_nodes_x<openmp_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_x<openmp_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_x<openmp_exec, STRUCTURED_RECTILINEAR_MESH>();
  check_for_all_nodes_x<openmp_exec, PARTICLE_MESH>();
  check_for_all_nodes_x<openmp_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>();
  check_for_all_nodes_x<openmp_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && \
  defined(RAJA_ENABLE_CUDA) && defined(AXOM_USE_UMPIRE)

  using cuda_exec = axom::CUDA_EXEC<512>;

  // Use unified memory
  const int exec_space_id = axom::getUmpireResourceAllocatorID(
    umpire::resource::MemoryResourceType::Unified);
  const int prev_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(exec_space_id);

  check_for_all_nodes_x<cuda_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_x<cuda_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_x<cuda_exec, STRUCTURED_RECTILINEAR_MESH>();
  check_for_all_nodes_x<cuda_exec, PARTICLE_MESH>();
  check_for_all_nodes_x<cuda_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>();
  check_for_all_nodes_x<cuda_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>();

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

  check_for_all_nodes_x<hip_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_x<hip_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_x<hip_exec, STRUCTURED_RECTILINEAR_MESH>();
  check_for_all_nodes_x<hip_exec, PARTICLE_MESH>();
  check_for_all_nodes_x<hip_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>();
  check_for_all_nodes_x<hip_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>();

  setDefaultAllocator(prev_allocator);
#endif
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST(mint_execution_node_traversals, for_all_nodes_ijk)
{
  using seq_exec = axom::SEQ_EXEC;
  check_for_all_nodes_ijk<seq_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_ijk<seq_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_ijk<seq_exec, STRUCTURED_RECTILINEAR_MESH>();

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && \
  defined(RAJA_ENABLE_OPENMP)

  using openmp_exec = axom::OMP_EXEC;
  check_for_all_nodes_ijk<openmp_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_ijk<openmp_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_ijk<openmp_exec, STRUCTURED_RECTILINEAR_MESH>();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && \
  defined(RAJA_ENABLE_CUDA) && defined(AXOM_USE_UMPIRE)

  using cuda_exec = axom::CUDA_EXEC<512>;

  // Use unified memory
  const int exec_space_id = axom::getUmpireResourceAllocatorID(
    umpire::resource::MemoryResourceType::Unified);
  const int prev_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(exec_space_id);

  check_for_all_nodes_ijk<cuda_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_ijk<cuda_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_ijk<cuda_exec, STRUCTURED_RECTILINEAR_MESH>();

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

  check_for_all_nodes_ijk<hip_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_ijk<hip_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_ijk<hip_exec, STRUCTURED_RECTILINEAR_MESH>();

  setDefaultAllocator(prev_allocator);
#endif
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST(mint_execution_node_traversals, for_all_nodes_ij)
{
  using seq_exec = axom::SEQ_EXEC;
  check_for_all_nodes_ij<seq_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_ij<seq_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_ij<seq_exec, STRUCTURED_RECTILINEAR_MESH>();

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && \
  defined(RAJA_ENABLE_OPENMP)

  using openmp_exec = axom::OMP_EXEC;
  check_for_all_nodes_ij<openmp_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_ij<openmp_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_ij<openmp_exec, STRUCTURED_RECTILINEAR_MESH>();

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && \
  defined(RAJA_ENABLE_CUDA) && defined(AXOM_USE_UMPIRE)

  using cuda_exec = axom::CUDA_EXEC<512>;

  // Use unified memory
  const int exec_space_id = axom::getUmpireResourceAllocatorID(
    umpire::resource::MemoryResourceType::Unified);
  const int prev_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(exec_space_id);

  check_for_all_nodes_ij<cuda_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_ij<cuda_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_ij<cuda_exec, STRUCTURED_RECTILINEAR_MESH>();

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

  check_for_all_nodes_ij<hip_exec, STRUCTURED_UNIFORM_MESH>();
  check_for_all_nodes_ij<hip_exec, STRUCTURED_CURVILINEAR_MESH>();
  check_for_all_nodes_ij<hip_exec, STRUCTURED_RECTILINEAR_MESH>();

  setDefaultAllocator(prev_allocator);
#endif
}

//------------------------------------------------------------------------------
AXOM_CUDA_TEST(mint_execution_node_traversals, for_all_nodes_index)
{
  constexpr int NDIMS = 3;
  for(int i = 1; i <= NDIMS; ++i)
  {
    using seq_exec = axom::SEQ_EXEC;
    check_for_all_nodes_idx<seq_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_nodes_idx<seq_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_nodes_idx<seq_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_nodes_idx<seq_exec, PARTICLE_MESH>(i);
    check_for_all_nodes_idx<seq_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_nodes_idx<seq_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && \
  defined(RAJA_ENABLE_OPENMP)

    using omp_exec = axom::OMP_EXEC;
    check_for_all_nodes_idx<omp_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_nodes_idx<omp_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_nodes_idx<omp_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_nodes_idx<omp_exec, PARTICLE_MESH>(i);
    check_for_all_nodes_idx<omp_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_nodes_idx<omp_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && \
  defined(RAJA_ENABLE_CUDA) && defined(AXOM_USE_UMPIRE)

    using cuda_exec = axom::CUDA_EXEC<512>;

    // Use unified memory
    const int exec_space_id = axom::getUmpireResourceAllocatorID(
      umpire::resource::MemoryResourceType::Unified);
    const int prev_allocator = axom::getDefaultAllocatorID();
    axom::setDefaultAllocator(exec_space_id);

    check_for_all_nodes_idx<cuda_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_nodes_idx<cuda_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_nodes_idx<cuda_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_nodes_idx<cuda_exec, PARTICLE_MESH>(i);
    check_for_all_nodes_idx<cuda_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_nodes_idx<cuda_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

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

    check_for_all_nodes_idx<hip_exec, STRUCTURED_UNIFORM_MESH>(i);
    check_for_all_nodes_idx<hip_exec, STRUCTURED_CURVILINEAR_MESH>(i);
    check_for_all_nodes_idx<hip_exec, STRUCTURED_RECTILINEAR_MESH>(i);
    check_for_all_nodes_idx<hip_exec, PARTICLE_MESH>(i);
    check_for_all_nodes_idx<hip_exec, UNSTRUCTURED_MESH, SINGLE_SHAPE>(i);
    check_for_all_nodes_idx<hip_exec, UNSTRUCTURED_MESH, MIXED_SHAPE>(i);

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

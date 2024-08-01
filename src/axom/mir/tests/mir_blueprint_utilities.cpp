// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/mir.hpp"
#include "axom/mir/blueprint_utilities.hpp"

#include <iostream>

// RAJA
#if defined(AXOM_USE_RAJA)
  #include "RAJA/RAJA.hpp"
#endif

// clang-format off
#if defined (AXOM_USE_RAJA) && defined (AXOM_USE_UMPIRE)
  using seq_exec = axom::SEQ_EXEC;

  #if defined(AXOM_USE_OPENMP)
    using omp_exec = axom::OMP_EXEC;
  #else
    using omp_exec = seq_exec;
  #endif

  #if defined(AXOM_USE_CUDA) && defined(__CUDACC__)
    constexpr int CUDA_BLOCK_SIZE = 256;
    using cuda_exec = axom::CUDA_EXEC<CUDA_BLOCK_SIZE>;
  #else
    using cuda_exec = seq_exec;
  #endif

  #if defined(AXOM_USE_HIP)
    constexpr int HIP_BLOCK_SIZE = 64;
    using hip_exec = axom::HIP_EXEC<HIP_BLOCK_SIZE>;
  #else
    using hip_exec = seq_exec;
  #endif
#endif
// clang-format on


template <typename ExecSpace>
void test_conduit_allocate()
{
  axom::mir::utilities::blueprint::ConduitAllocateThroughAxom<ExecSpace> c2a;
  EXPECT_TRUE(c2a.getConduitAllocatorID() > 0);

  constexpr int nValues = 100;
  conduit::Node n;
  n.set_allocator(c2a.getConduitAllocatorID());
  n.set(conduit::DataType::int32(nValues));

  // Make sure we can store some values into the data that were allocated.
  int *ptr = static_cast<int *>(n.data_ptr());
  axom::for_all<ExecSpace>(nValues, AXOM_LAMBDA(auto index)
  {
    ptr[index] = index;
  });

  EXPECT_EQ(n.dtype().number_of_elements(), nValues);
}

TEST(mir_blueprint_utilities, allocate)
{
  test_conduit_allocate<seq_exec>(); 
#if defined(AXOM_USE_OPENMP)
  test_conduit_allocate<omp_exec>();
#endif
#if defined(AXOM_USE_CUDA) && defined(__CUDACC__)
  test_conduit_allocate<cuda_exec>();
#endif
#if defined(AXOM_USE_HIP)
  test_conduit_allocate<hip_exec>();
#endif
}

template <typename ExecSpace>
void test_copy_braid(const conduit::Node &mesh)
{
  // Copy the mesh to device.
  conduit::Node mesh_dev;
  axom::mir::utilities::blueprint::copy<ExecSpace>(mesh_dev, mesh);

  // Run some minmax operations on device (proves that the data was in the right place) and check the results.

  constexpr double eps = 1.e-7;

  auto x = axom::mir::utilities::blueprint::minmax<ExecSpace>(mesh["coordsets/coords/values/x"]);
  //std::cout << std::setw(16) << "x={" << x.first << ", " << x.second << "}\n";
  EXPECT_NEAR(x.first, -10., eps);
  EXPECT_NEAR(x.second, 10., eps);

  auto y = axom::mir::utilities::blueprint::minmax<ExecSpace>(mesh["coordsets/coords/values/y"]);
  //std::cout << std::setw(16) << "y={" << y.first << ", " << y.second << "}\n";
  EXPECT_NEAR(y.first, -10., eps);
  EXPECT_NEAR(y.second, 10., eps);

  auto c = axom::mir::utilities::blueprint::minmax<ExecSpace>(mesh["topologies/mesh/elements/connectivity"]);
  //std::cout << std::setw(16) << "conn={" << c.first << ", " << c.second << "}\n";
  EXPECT_NEAR(c.first, 0., eps);
  EXPECT_NEAR(c.second, 999., eps);

  auto r = axom::mir::utilities::blueprint::minmax<ExecSpace>(mesh["fields/radial/values"]);
  //std::cout << std::setw(16) << "radial={" << r.first << ", " << r.second << "}\n";
  EXPECT_NEAR(r.first, 19.2450089729875, eps);
  EXPECT_NEAR(r.second, 173.205080756888, eps);
}

TEST(mir_blueprint_utilities, copy)
{
  const int d[3] = {10, 10, 10};
  conduit::Node mesh;
  conduit::blueprint::mesh::examples::braid("hexs", d[0], d[1], d[2], mesh);

  test_copy_braid<seq_exec>(mesh);
#if defined(AXOM_USE_OPENMP)
  test_copy_braid<omp_exec>(mesh);
#endif
#if defined(AXOM_USE_CUDA) && defined(__CUDACC__)
  test_copy_braid<cuda_exec>(mesh);
#endif
#if defined(AXOM_USE_HIP)
  test_copy_braid<hip_exec>(mesh);
#endif
}

TEST(mir_blueprint_utilities, to_unstructured)
{
   // TODO: to_unstructured
}

//------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;  // create & initialize test logger,

  result = RUN_ALL_TESTS();
  return result;
}

// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/config.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/Array.hpp"
#include "axom/core/execution/execution_space.hpp"
#include "axom/core/execution/for_all.hpp"

// gtest includes
#include "gtest/gtest.h"

//------------------------------------------------------------------------------
//  HELPER METHODS
//------------------------------------------------------------------------------
namespace
{
//------------------------------------------------------------------------------
template <typename TheExecSpace>
class core_array_for_all : public ::testing::Test
{
public:
  using ExecSpace = TheExecSpace;

  // Define some memory spaces
#ifdef AXOM_USE_UMPIRE
  static constexpr axom::MemorySpace host_memory = axom::MemorySpace::Host;
#else
  static constexpr axom::MemorySpace host_memory = axom::MemorySpace::Dynamic;
#endif

  static constexpr axom::MemorySpace exec_space_memory =
    axom::execution_space<ExecSpace>::memory_space;

  // Define some Array type aliases
  using HostArray = axom::Array<int, 1, host_memory>;
  using DynamicArray = axom::Array<int, 1, axom::MemorySpace::Dynamic>;
  using KernelArray = axom::Array<int, 1, exec_space_memory>;
  using KernelArrayView = axom::ArrayView<int, 1, exec_space_memory>;
};

// Generate a list of available execution types
using MyTypes = ::testing::Types<
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP)
  axom::OMP_EXEC,
#endif
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && defined(AXOM_USE_UMPIRE)
  axom::CUDA_EXEC<256>,
#endif
  axom::SEQ_EXEC>;

TYPED_TEST_SUITE(core_array_for_all, MyTypes);

}  // end anonymous namespace

AXOM_TYPED_TEST(core_array_for_all, explicit_ArrayView)
{
  using ExecSpace = typename TestFixture::ExecSpace;
  using KernelArray = typename TestFixture::KernelArray;
  using KernelArrayView = typename TestFixture::KernelArrayView;
  using HostArray = typename TestFixture::HostArray;

  // Create an array of N items using default MemorySpace for ExecSpace
  constexpr int N = 374;
  KernelArray arr(N);

  // Modify array using mutable lambda and ArrayView
  KernelArrayView arr_view(arr);
  axom::for_all<ExecSpace>(
    N,
    AXOM_LAMBDA(axom::IndexType idx) mutable { arr_view[idx] = N - idx; });

  // Check array contents on device
  HostArray localArr = arr;
  for(int i = 0; i < N; ++i)
  {
    EXPECT_EQ(localArr[i], N - i);
  }
}

AXOM_TYPED_TEST(core_array_for_all, auto_ArrayView)
{
  using ExecSpace = typename TestFixture::ExecSpace;
  using KernelArray = typename TestFixture::KernelArray;
  using HostArray = typename TestFixture::HostArray;

  // Create an array of N items using default MemorySpace for ExecSpace
  constexpr int N = 374;
  KernelArray arr(N);

  // Modify array using mutable lambda and ArrayView
  auto arr_view = arr.view();
  axom::for_all<ExecSpace>(
    N,
    AXOM_LAMBDA(axom::IndexType idx) mutable { arr_view[idx] = N - idx; });

  // Check array contents on device
  HostArray localArr = arr;
  for(int i = 0; i < N; ++i)
  {
    EXPECT_EQ(localArr[i], N - i);
  }
}

AXOM_TYPED_TEST(core_array_for_all, dynamic_array)
{
  using ExecSpace = typename TestFixture::ExecSpace;
  using DynamicArray = typename TestFixture::DynamicArray;
  using HostArray = typename TestFixture::HostArray;

  int kernelAllocID = axom::execution_space<ExecSpace>::allocatorID();
  int hostAllocID = axom::execution_space<axom::SEQ_EXEC>::allocatorID();

  // Create an array of N items using default MemorySpace for ExecSpace
  constexpr axom::IndexType N = 374;
  DynamicArray arr(N, N, kernelAllocID);

  // Modify array using mutable lambda and ArrayView
  auto arr_view = arr.view();
  axom::for_all<ExecSpace>(
    N,
    AXOM_LAMBDA(axom::IndexType idx) mutable { arr_view[idx] = N - idx; });

  // Check array contents on device
  HostArray localArr(arr, hostAllocID);
  for(int i = 0; i < N; ++i)
  {
    EXPECT_EQ(localArr[i], N - i);
  }
}
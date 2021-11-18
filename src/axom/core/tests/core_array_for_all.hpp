// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/config.hpp"
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
template <typename ExecSpace>
void check_array_for_all_templated_memory()
{
  // Introduce some type aliases
#ifdef AXOM_USE_UMPIRE
  constexpr axom::MemorySpace host_memory = axom::MemorySpace::Host;
#else
  constexpr axom::MemorySpace host_memory = axom::MemorySpace::Dynamic;
#endif

  constexpr axom::MemorySpace exec_space_memory =
    axom::execution_space<ExecSpace>::memory_space;

  using HostArray = axom::Array<int, 1, host_memory>;
  using KernelArray = axom::Array<int, 1, exec_space_memory>;
  using KernelArrayView = axom::ArrayView<int, 1, exec_space_memory>;

  // Create an array of N items using default allocator ExecSpace
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

}  // namespace

//------------------------------------------------------------------------------
//  UNIT TESTS
//------------------------------------------------------------------------------
TEST(core_array_for_all, seq_exec)
{
  check_array_for_all_templated_memory<axom::SEQ_EXEC>();
}

//------------------------------------------------------------------------------

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP)
TEST(core_array_for_all, omp_exec)
{
  check_array_for_all_templated_memory<axom::OMP_EXEC>();
}
#endif

//------------------------------------------------------------------------------

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && defined(AXOM_USE_UMPIRE)
TEST(core_array_for_all, cuda_exec)
{
  constexpr int BLOCK_SIZE = 256;
  check_array_for_all_templated_memory<axom::CUDA_EXEC<BLOCK_SIZE>>();
}
#endif

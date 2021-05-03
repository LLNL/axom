// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/config.hpp"                         /* for compile time defs */
#include "axom/core/Macros.hpp"                    /* for axom macros */
#include "axom/core/execution/execution_space.hpp" /* execution_space traits */
#include "axom/core/execution/for_all.hpp"         /* for_all() traversals */
#include "axom/core/execution/synchronize.hpp"     /* synchronize() */
#include "axom/core/memory_management.hpp"         /* allocate/deallocate */

// gtest includes
#include "gtest/gtest.h"

//------------------------------------------------------------------------------
//  HELPER METHODS
//------------------------------------------------------------------------------
namespace
{
//------------------------------------------------------------------------------
template <typename ExecSpace>
void check_for_all()
{
  EXPECT_TRUE(axom::execution_space<ExecSpace>::valid());
  std::cout << "checking axom::for_all() with ["
            << axom::execution_space<ExecSpace>::name() << "]\n";

  // STEP 0: set some constants
  constexpr int VALUE_1 = -42;
  constexpr int VALUE_2 = 42;
  constexpr int N = 256;

  // STEP 1: set default allocator for the execution space
  const int currentAllocatorID = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(axom::execution_space<ExecSpace>::allocatorID());

  // STEP 0: allocate buffer
  int* a = axom::allocate<int>(N);

  // STEP 1: initialize to VALUE_1
  axom::for_all<ExecSpace>(
    N,
    AXOM_LAMBDA(axom::IndexType idx) { a[idx] = VALUE_1; });

  if(axom::execution_space<ExecSpace>::async())
  {
    axom::synchronize<ExecSpace>();
  }

  // STEP 2: check array
  for(int i = 0; i < N; ++i)
  {
    EXPECT_EQ(a[i], VALUE_1);
  }

  // STEP 3: add VALUE_2 to all entries resulting to zero
  axom::for_all<ExecSpace>(
    0,
    N,
    AXOM_LAMBDA(axom::IndexType idx) { a[idx] += VALUE_2; });

  if(axom::execution_space<ExecSpace>::async())
  {
    axom::synchronize<ExecSpace>();
  }

  // STEP 4: check result
  for(int i = 0; i < N; ++i)
  {
    EXPECT_EQ(a[i], 0);
  }

  // STEP 5: cleanup
  axom::deallocate(a);

  axom::setDefaultAllocator(currentAllocatorID);
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
//  UNIT TESTS
//------------------------------------------------------------------------------
TEST(core_execution_for_all, seq_exec) { check_for_all<axom::SEQ_EXEC>(); }

//------------------------------------------------------------------------------

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP)

TEST(core_execution_for_all, omp_exec) { check_for_all<axom::OMP_EXEC>(); }

#endif

//------------------------------------------------------------------------------

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && defined(AXOM_USE_UMPIRE)

TEST(core_execution_for_all, cuda_exec)
{
  constexpr int BLOCK_SIZE = 256;
  check_for_all<axom::CUDA_EXEC<BLOCK_SIZE>>();
}

//------------------------------------------------------------------------------
TEST(core_execution_for_all, cuda_exec_async)
{
  constexpr int BLOCK_SIZE = 256;
  check_for_all<axom::CUDA_EXEC<BLOCK_SIZE, axom::ASYNC>>();
}

#endif

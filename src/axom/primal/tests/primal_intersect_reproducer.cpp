// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"

#include "axom/core/execution/execution_space.hpp"
#include "axom/core/memory_management.hpp"

#include "axom/primal/geometry/Hexahedron.hpp"
#include "axom/primal/geometry/Tetrahedron.hpp"

#include "axom/primal/operators/intersect.hpp"

#include <cmath>

namespace primal = axom::primal;

namespace
{
// #if defined(AXOM_USE_UMPIRE) && defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA)
#if defined(AXOM_USE_CUDA) && defined(AXOM_USE_CUDA)

void reproducer() {
  constexpr int BLOCK_SIZE = 256;
  using exec = axom::CUDA_EXEC<BLOCK_SIZE>;

  //RAJA::ReduceMin<RAJA::cuda_reduce, int> reducer(0) ;

  axom::for_all<ExecSpace>(
    4,
    AXOM_LAMBDA(int i) {
//  RAJA::forall<RAJA::cuda_exec<256, true>>(RAJA::TypedRangeSegment<int>(0, 1), [=] __device__ (int i) {
      axom::primal::Hexahedron<double, 3> hex {};
      axom::primal::Tetrahedron<double, 3> tet {};
      axom::primal::intersection_volume(hex, tet, 1.0e-11);

      //reducer.min(i);
    });
}


AXOM_CUDA_TEST(primal_intersect_reproducer, tet_hex_intersect_reproducer)
{
  reproducer();
}


#endif /* AXOM_USE_RAJA && AXOM_USE_UMPIRE */

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger(axom::slic::message::Warning);

  int result = RUN_ALL_TESTS();
  return result;
}

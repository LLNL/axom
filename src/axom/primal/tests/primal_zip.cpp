// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/core/execution/for_all.hpp"
#include "axom/core/memory_management.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/utils/ZipIndexable.hpp"
#include "axom/primal/utils/ZipPoint.hpp"
#include "axom/primal/utils/ZipVector.hpp"

#include "gtest/gtest.h"

using namespace axom;

//------------------------------------------------------------------------------

template <typename ExecSpace, typename PrimitiveType>
void check_zip_points_3d()
{
  using ZipType = primal::ZipIndexable<PrimitiveType>;

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(axom::execution_space<ExecSpace>::allocatorID());

  // create arrays of data
  constexpr int N = 8;
  double* x = axom::allocate<double>(N);
  double* y = axom::allocate<double>(N);
  double* z = axom::allocate<double>(N);

  bool* valid = axom::allocate<bool>(N);

  axom::for_all<ExecSpace>(
    0,
    8,
    AXOM_LAMBDA(int idx) {
      x[idx] = (idx % 2) + 1;
      y[idx] = ((idx / 2) % 2) + 1;
      z[idx] = (idx / 4) + 1;
    });

  ZipType it {x, y, z};

  axom::for_all<ExecSpace>(
    0,
    8,
    AXOM_LAMBDA(int idx) {
      valid[idx] = (it[idx] == PrimitiveType {x[idx], y[idx], z[idx]});
    });

  for(int i = 0; i < N; i++)
  {
    EXPECT_EQ(valid[i], true);
  }

  axom::deallocate(x);
  axom::deallocate(y);
  axom::deallocate(z);
  axom::deallocate(valid);
  axom::setDefaultAllocator(current_allocator);
}

TEST(primal_zip, zip_points_3d)
{
  using PointType = primal::Point<double, 3>;
  using ExecSpace = axom::SEQ_EXEC;

  check_zip_points_3d<ExecSpace, PointType>();
}

TEST(primal_zip, zip_vectors_3d)
{
  using PointType = primal::Vector<double, 3>;
  using ExecSpace = axom::SEQ_EXEC;

  check_zip_points_3d<ExecSpace, PointType>();
}

#ifdef AXOM_USE_CUDA
AXOM_CUDA_TEST(primal_zip, zip_points_3d_gpu)
{
  using PointType = primal::Point<double, 3>;
  using ExecSpace = axom::CUDA_EXEC<256>;

  check_zip_points_3d<ExecSpace, PointType>();
}

AXOM_CUDA_TEST(primal_zip, zip_vectors_3d_gpu)
{
  using PointType = primal::Vector<double, 3>;
  using ExecSpace = axom::CUDA_EXEC<256>;

  check_zip_points_3d<ExecSpace, PointType>();
}
#endif

//------------------------------------------------------------------------------
#include "axom/slic/core/SimpleLogger.hpp"
using axom::slic::SimpleLogger;

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  SimpleLogger logger;  // create & initialize test logger,

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}

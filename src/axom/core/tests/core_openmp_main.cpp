// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"  // for compile-time definitions

#include "core_openmp_map.hpp"



int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);

  // add this line to avoid a warning in the output about thread safety
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

#ifdef EIGEN_SOLVE_TESTER_SHOULD_SEED
  std::srand(std::time(0));
#else
  std::srand(42);
#endif

  return RUN_ALL_TESTS();
}

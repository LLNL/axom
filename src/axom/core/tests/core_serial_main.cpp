// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"  // for compile-time definitions

#include "core_array.hpp"
#include "core_execution_for_all.hpp"
#include "core_execution_space.hpp"
#include "core_memory_management.hpp"
#include "core_Path.hpp"
#include "core_stack_array.hpp"

#ifndef AXOM_USE_MPI
  #include "core_types.hpp"
#endif

#include "numerics_determinants.hpp"
#include "numerics_eigen_solve.hpp"
#include "numerics_eigen_sort.hpp"
#include "numerics_floating_point_limits.hpp"
#include "numerics_jacobi_eigensolve.hpp"
#include "numerics_linear_solve.hpp"
#include "numerics_lu.hpp"
#include "numerics_matrix.hpp"
#include "numerics_polynomial_solvers.hpp"
#include "numerics_matvecops.hpp"

#include "utils_Timer.hpp"
#include "utils_endianness.hpp"
#include "utils_fileUtilities.hpp"
#include "utils_nvtx_settings.hpp"
#include "utils_utilities.hpp"
#include "utils_about.hpp"

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

// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core/utilities/Utilities.hpp"

// C/C++ includes
#include <type_traits>

//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------
template <typename T>
void check_random_real(int offset)
{
  AXOM_STATIC_ASSERT(std::is_floating_point<T>::value);

  constexpr T min = 0;
  constexpr T max = 1;

  const T curr_min = min - static_cast<T>(offset);
  const T curr_max = max - static_cast<T>(offset);

  constexpr int MAX_ITERS = 100;
  for(int i = 0; i < MAX_ITERS; ++i)
  {
    T val = axom::utilities::random_real(curr_min, curr_max);
    EXPECT_GE(val, curr_min);
    EXPECT_LT(val, curr_max);
  }
}

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST(utils_Utilities, log2)
{
  std::cout << "Testing log2 functions." << std::endl;

  // Test integer log2 value of type int
  {
    int val = 64;
    int exp = 6;
    EXPECT_EQ(exp, axom::utilities::log2(val));
  }

  // Test non-integer log2 value of type int
  {
    int val = 72;  // not a power of 2
    int exp = 6;
    EXPECT_EQ(exp, axom::utilities::log2(val));
  }

  // Test integer log2 value of type double
  {
    double val = 16.;
    double exp = 4.;
    EXPECT_EQ(exp, axom::utilities::log2(val));
  }

  // Test non-integer log2 value of type double
  {
    double val = 20.;  // not a power of 2
    double exp = 4.3219281;
    EXPECT_NEAR(exp, axom::utilities::log2(val), 1e-5);
  }
}

//------------------------------------------------------------------------------
TEST(utils_Utilities, random_real)
{
  std::cout << "Testing random_real functions (non-deterministic)." << std::endl;

  for(int offset = 0; offset < 10; ++offset)
  {
    check_random_real<float>(offset);
    check_random_real<double>(offset);
    check_random_real<long double>(offset);
  }
}

//------------------------------------------------------------------------------
TEST(utils_Utilities, random_real_with_seed)
{
  std::cout << "Testing random_real functions (deterministic)." << std::endl;
  constexpr unsigned int seed = 123456789;

  constexpr double a = -5.0;
  constexpr double b = 5.0;

  const double expected_reals[5] = {-1.5112829544380526,
                                    -2.3311429024686219,
                                    -3.6335370551231403,
                                    -4.714431326610093,
                                    3.6893326916732878};

  for(int i = 0; i < 5; ++i)
  {
    const double real = axom::utilities::random_real(a, b, seed);
    EXPECT_DOUBLE_EQ(real, expected_reals[i]);
    EXPECT_GE(real, a);
    EXPECT_LT(real, b);
  }
}

//------------------------------------------------------------------------------
TEST(utils_Utilities, minmax)
{
  std::cout << "Testing min and max functions." << std::endl;

  // Test simple min, max comparisons on ints
  {
    int a = 5;
    int b = 7;

    int temp_min = axom::utilities::min(a, b);
    int temp_max = axom::utilities::max(a, b);

    EXPECT_EQ(a, temp_min);
    EXPECT_EQ(b, temp_max);
  }

  // Test simple min, max comparisons on doubles
  {
    double a = 5.2;
    double b = -1.7;

    int temp_min = axom::utilities::min(a, b);
    int temp_max = axom::utilities::max(a, b);

    EXPECT_EQ(b, temp_min);
    EXPECT_EQ(a, temp_max);
  }
}

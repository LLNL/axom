// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core/utilities/Utilities.hpp"

#include <algorithm> // for std::is_permutation

TEST(core_Utilities,log2)
{
  std::cout<<"Testing log2 functions."<< std::endl;

  // Test integer log2 value of type int
  {
    int val = 64;
    int exp = 6;
    EXPECT_EQ(exp, axom::utilities::log2(val));
  }

  // Test non-integer log2 value of type int
  {
    int val = 72; // not a power of 2
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
    double val = 20.; // not a power of 2
    double exp = 4.3219281;
    EXPECT_NEAR(exp, axom::utilities::log2(val), 1e-5);
  }
}

TEST(core_Utilities,random_real)
{
  std::cout<<"Testing random_real functions (non-deterministic)."<< std::endl;

  int min = 0;
  int max = 1;
  for (int offset = 0 ; offset < 10 ; ++offset)
  {
    int cur_min = min - offset;
    int cur_max = max + offset;
    for (int i = 0 ; i < 100 ; ++i)
    {
      float f_val = axom::utilities::random_real<float>(cur_min, cur_max);
      EXPECT_GE(f_val, cur_min);
      EXPECT_LT(f_val, cur_max);

      double d_val = axom::utilities::random_real<double>(cur_min, cur_max);
      EXPECT_GE(d_val, cur_min);
      EXPECT_LT(d_val, cur_max);

      long double ld_val = axom::utilities::random_real<long double>(cur_min,
                                                                     cur_max);
      EXPECT_GE(ld_val, cur_min);
      EXPECT_LT(ld_val, cur_max);
    }
  }
}

TEST( core_Utilities,random_real_with_seed )
{
  std::cout<<"Testing random_real functions (deterministic)."<< std::endl;
  constexpr unsigned int seed = 123456789;

  constexpr double a = -5.0;
  constexpr double b =  5.0;

  const double expected_reals[ 5 ] = {
    -1.5112829544380526,
    -2.3311429024686219,
    -3.6335370551231403,
    -4.714431326610093,
    3.6893326916732878
  };

  for ( int i=0 ; i < 5 ; ++i )
  {
    const double real = axom::utilities::random_real( a, b, seed );
    EXPECT_DOUBLE_EQ( real, expected_reals[ i ] );
    EXPECT_GE(real, a);
    EXPECT_LT(real, b);
  }

}

TEST(core_Utilities,minmax)
{
  std::cout<<"Testing min and max functions."<< std::endl;

  // Test simple min, max comparisons on ints
  {
    int a = 5;
    int b = 7;

    EXPECT_EQ(a, axom::utilities::min(a,b));
    EXPECT_EQ(b, axom::utilities::max(a,b));
  }

  // Test simple min, max comparisons on doubles
  {
    double a = 5.2;
    double b = -1.7;

    EXPECT_EQ(b, axom::utilities::min(a,b));
    EXPECT_EQ(a, axom::utilities::max(a,b));
  }
}

TEST(core_Utilities, isValidPermutation)
{
  int const referenceValues[ 5 ] = { 0, 1, 2, 3, 4 };
  int values[ 5 ];

  for ( int i = 0; i < 5; ++i )
  {
    values[ 0 ] = i; 

    for ( int j = 0; j < 5; ++j )
    {
      values[ 1 ] = j; 

      for ( int k = 0; k < 5; ++k )
      {
        values[ 2 ] = k; 

        for ( int l = 0; l < 5; ++l )
        {
          values[ 3 ] = l;

          for ( int m = 0; m < 5; ++m )
          {
            values[ 4 ] = m;

            for ( int size = 0; size < 5; ++ size )
            {
              EXPECT_EQ(axom::utilities::isValidPermutation(values, size),
                        std::is_permutation(values, values + size, referenceValues));
            }
          }
        }
      }
    }
  }

  {
    constexpr int NVALS = 1;
    int badValues[NVALS] = {-5};
    EXPECT_FALSE(axom::utilities::isValidPermutation(badValues, NVALS));
  }

  {
    constexpr int NVALS = 1;
    int badValues[NVALS] = {1};
    EXPECT_FALSE(axom::utilities::isValidPermutation(badValues, NVALS));
  }
  
  {
    constexpr int NVALS = 2;
    int badValues[NVALS] = {-1, 0};
    EXPECT_FALSE(axom::utilities::isValidPermutation(badValues, NVALS));
  }

  {
    constexpr int NVALS = 2;
    int badValues[NVALS] = {3, 0};
    EXPECT_FALSE(axom::utilities::isValidPermutation(badValues, NVALS));
  }

  {
    constexpr int NVALS = 2;
    int badValues[NVALS] = {0, 0};
    EXPECT_FALSE(axom::utilities::isValidPermutation(badValues, NVALS));
  }

  {
    constexpr int NVALS = 3;
    int badValues[NVALS] = {0, 1, 55};
    EXPECT_FALSE(axom::utilities::isValidPermutation(badValues, NVALS));
  }

  {
    constexpr int NVALS = 3;
    int badValues[NVALS] = {0, 1, 0};
    EXPECT_FALSE(axom::utilities::isValidPermutation(badValues, NVALS));
  }

  {
    constexpr int NVALS = 3;
    int badValues[NVALS] = {0, -3, 1};
    EXPECT_FALSE(axom::utilities::isValidPermutation(badValues, NVALS));
  }
}

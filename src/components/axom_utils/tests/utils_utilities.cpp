/*
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


#include "gtest/gtest.h"

#include "axom_utils/Utilities.hpp"


TEST(axom_utils_Utilities,log2)
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

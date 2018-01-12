/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#include "gtest/gtest.h"

#include "axom_utils/Utilities.hpp"

TEST(axom_utils_Utilities, allocation)
{ 
  std::cout<<"Testing allocation functions."<< std::endl;  
  for ( int initial_size = 2; initial_size <= 1048576; initial_size *= 2 )
  {
    int buffer_size = initial_size;
    int * buffer = axom::utilities::alloc<int>(buffer_size);

    for (int i = 0; i < buffer_size; i++)
    {
      buffer[i] = i;
    }

    buffer_size *= 2;
    buffer = axom::utilities::realloc(buffer, buffer_size);
    for (int i = 0; i < buffer_size / 2; i++)
    {
      EXPECT_EQ(buffer[i],  i);
    }

    buffer_size /= 4;
    buffer = axom::utilities::realloc(buffer, buffer_size);
    for (int i = 0; i < buffer_size; i++)
    {
      EXPECT_EQ(buffer[i],  i);
    }

    axom::utilities::free(buffer);
  }
}


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

TEST(axom_utils_Utilities,minmax)
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

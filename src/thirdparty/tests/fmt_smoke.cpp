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

//-----------------------------------------------------------------------------
//
// file: fmt_smoke.cpp
// A simple test to see if the thirdparty fmt libary is working properly.
// Note: To ensure that fmt is configured properly, included the `fmt/fmt.hpp' helper file
// rather than directly including fmt's header files.
//
//-----------------------------------------------------------------------------

#include "fmt/fmt.hpp"

#include "gtest/gtest.h"

//-----------------------------------------------------------------------------
TEST(fmt_smoke, basic_use)
{
  // Test python style formatting -- both examples should produce 'hello world'
  std::string hw1 = fmt::format("{} {}", "hello", "world");
  EXPECT_EQ("hello world", hw1);

  std::string hw2 = fmt::format("{1} {0}", "world", "hello");
  EXPECT_EQ("hello world", hw2);

  // Test printf-style formatting.
  // It should produce two significant digits -- i.e. 1.23
  std::string fltFormat = fmt::sprintf("%.2f", 1.234567);
  EXPECT_EQ( "1.23", fltFormat);
}

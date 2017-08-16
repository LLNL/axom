/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

//-----------------------------------------------------------------------------
///
/// file: fmt_smoke.cpp
/// A simple test to see if the thirdparty fmt libary is working properly
///
//-----------------------------------------------------------------------------

#include "fmt/format.h"

#include "gtest/gtest.h"

//-----------------------------------------------------------------------------
TEST(fmt_smoke, basic_use)
{
    // Test python style formatting -- both examples should produce 'hello world'
    std::string hw1 = fmt::format("{} {}", "hello", "world");
    EXPECT_EQ("hello world", hw1);

    std::string hw2 = fmt::format("{1} {0}", "world", "hello");
    EXPECT_EQ("hello world", hw2);

    // Test printf-style formatting -- should produce two significant digits -- i.e. 1.23
    std::string fltFormat = fmt::sprintf("%.2f", 1.234567);
    EXPECT_EQ( "1.23", fltFormat);
}

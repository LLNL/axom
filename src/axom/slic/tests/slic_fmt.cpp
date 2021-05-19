// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file
 *
 * \brief A simple test to see if we can use the fmt library within slic macros.
 */

#include "axom/config.hpp"
#include "fmt/fmt.hpp"

#include "axom/slic/interface/slic.hpp"
#include "axom/slic/core/SimpleLogger.hpp"

using axom::slic::SimpleLogger;

#include "gtest/gtest.h"

//-----------------------------------------------------------------------------
TEST(slic_fmt, basic_use)
{
  SLIC_INFO("Formatting with C++ streams and with "
            << fmt::format("the '{}' string formatting library", "fmt")
            << " can be used together in slic macros.");

  // Style similar to Python format
  SLIC_INFO(fmt::format("{1} {0}", "world", "Hello"));
  EXPECT_EQ("Hello world", fmt::format("{1} {0}", "world", "Hello"));

#if !defined(__INTEL_COMPILER) && __INTEL_COMPILER != 1800
  // Note: This test seg faults on intel 18

  // Python style w/ 'dictionary'
  SLIC_INFO(
    fmt::format("Hello, {name}! Goodbye, {name}.", fmt::arg("name", "Axom")));
  EXPECT_EQ(
    "Hello, Axom! Goodbye, Axom.",
    fmt::format("Hello, {name}! Goodbye, {name}.", fmt::arg("name", "Axom")));
#endif

  // Python style with additional formatting
  SLIC_INFO(fmt::format("int:{0:d};  hex:{0:x};  oct:{0:o}; bin:{0:b}", 42));
  EXPECT_EQ("int:42;  hex:2a;  oct:52; bin:101010",
            fmt::format("int:{0:d};  hex:{0:x};  oct:{0:o}; bin:{0:b}", 42));

  // sprintf style
  SLIC_INFO(fmt::sprintf("Two significant digits: %.2f", 1.234567));
  EXPECT_EQ("1.23", fmt::sprintf("%.2f", 1.234567));
}

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  SimpleLogger logger;  // create & initialize test logger,
                        // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}

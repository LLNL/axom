// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core/utilities/System.hpp"

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST(utils_system, getUserName)
{
  // Returns empty string on failure
  EXPECT_TRUE(axom::utilities::getUserName() != "");
}

TEST(utils_system, getHostName)
{
  // Returns empty string on failure
  EXPECT_TRUE(axom::utilities::getHostName() != "");
}

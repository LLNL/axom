// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <iostream>

#include "sol/sol.hpp"
#include "gtest/gtest.h"

TEST(sol_smoke, basic_use)
{
  sol::state lua;
  lua.script(
    "table1={"
    "  table2={"
    "    some_bool = true,"
    "    some_double = 3.0 "
    "  }"
    "}");

  bool some_bool = lua["table1"]["table2"]["some_bool"];
  EXPECT_TRUE(some_bool);

  double some_double = lua["table1"]["table2"]["some_double"];
  EXPECT_NEAR(some_double, 3.0, 0.1);
}

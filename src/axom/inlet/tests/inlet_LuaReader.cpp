// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include <string>
#include <vector>

#include "axom/inlet/LuaReader.hpp"

TEST(inlet_LuaReader_getBool, getTopLevelBools)
{
  axom::inlet::LuaReader lr;
  lr.parseString("foo = true; bar = false");

  bool value, retValue;

  value = false;
  retValue = lr.getBool("foo", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, true);

  value = true;
  retValue = lr.getBool("bar", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, false);
}


TEST(inlet_LuaReader_getBool, getInsideBools)
{
  axom::inlet::LuaReader lr;
  lr.parseString("foo = { bar = false; baz = true }");

  bool value, retValue;

  value = true;
  retValue = lr.getBool("foo/bar", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, false);

  value = false;
  retValue = lr.getBool("foo/baz", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, true);
}


TEST(inlet_LuaReader_getString, getTopLevelStrings)
{
  axom::inlet::LuaReader lr;
  lr.parseString("foo = \"this is a test string\"; bar = \"TesT StrInG\"");

  bool retValue;
  std::string value;

  value = "";
  retValue = lr.getString("foo", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, "this is a test string");

  value = "";
  retValue = lr.getString("bar", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, "TesT StrInG");
}


TEST(inlet_LuaReader_getString, getInsideStrings)
{
  axom::inlet::LuaReader lr;
  lr.parseString(
    "foo = { bar = \"this is a test string\"; baz = \"TesT StrInG\" }");

  bool retValue;
  std::string value;

  value = "";
  retValue = lr.getString("foo/bar", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, "this is a test string");

  value = "";
  retValue = lr.getString("foo/baz", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, "TesT StrInG");
}

//------------------------------------------------------------------------------
#include "axom/slic/core/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}

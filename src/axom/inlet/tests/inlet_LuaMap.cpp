// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include <string>
#include <vector>

#include "axom/inlet/LuaMap.hpp"

TEST(inlet_LuaMap_getBool, getTopLevelBools)
{
  axom::inlet::LuaMap lm;
  lm.parseString("foo = true; bar = false");

  bool value, retValue;

  value = false;
  retValue = lm.getBool("foo", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, true);

  value = true;
  retValue = lm.getBool("bar", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, false);
}


TEST(inlet_LuaMap_getBool, getInsideBools)
{
  axom::inlet::LuaMap lm;
  lm.parseString("foo = { bar = false; baz = true }");

  bool value, retValue;

  value = true;
  retValue = lm.getBool("foo.bar", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, false);

  value = false;
  retValue = lm.getBool("foo.baz", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, true);
}


TEST(inlet_LuaMap_getString, getTopLevelStrings)
{
  axom::inlet::LuaMap lm;
  lm.parseString("foo = \"this is a test string\"; bar = \"TesT StrInG\"");

  bool retValue;
  std::string value;

  value = "";
  retValue = lm.getString("foo", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, "this is a test string");

  value = "";
  retValue = lm.getString("bar", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, "TesT StrInG");
}


TEST(inlet_LuaMap_getString, getInsideStrings)
{
  axom::inlet::LuaMap lm;
  lm.parseString(
    "foo = { bar = \"this is a test string\"; baz = \"TesT StrInG\" }");

  bool retValue;
  std::string value;

  value = "";
  retValue = lm.getString("foo.bar", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, "this is a test string");

  value = "";
  retValue = lm.getString("foo.baz", value);
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

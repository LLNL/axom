// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include <string>
#include <vector>

#include "axom/inlet/LuaReader.hpp"

TEST(inlet_LuaReader, getTopLevelBools)
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

TEST(inlet_LuaReader, getInsideBools)
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

TEST(inlet_LuaReader, getTopLevelStrings)
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

TEST(inlet_LuaReader, getInsideStrings)
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

TEST(inlet_LuaReader, mixLevelTables)
{
  axom::inlet::LuaReader lr;
  lr.parseString("t = { innerT = { foo = 1 }, anotherInnerT = {baz = 3}}");

  bool retValue;
  int value;

  value = 0;
  retValue = lr.getInt("t/innerT/foo", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, 1);

  value = 0;
  retValue = lr.getInt("t/doesntexist", value);
  EXPECT_EQ(retValue, false);
  EXPECT_EQ(value, 0);

  value = 0;
  retValue = lr.getInt("t/anotherInnerT/baz", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, 3);
}

// Checks that LuaReader parses array information as expected
TEST(inlet_LuaReader, getMap)
{
  std::string testString =
    "luaArray = { [1] = 4, [2] = 5, [3] = 6 , [4] = true, [8] = false, [12] = "
    "2.4, [33] = 'hello', [200] = 'bye' }";
  axom::inlet::LuaReader lr;
  lr.parseString(testString);

  std::unordered_map<int, int> ints;
  bool found = lr.getIntMap("luaArray", ints);
  EXPECT_TRUE(found);
  std::unordered_map<int, int> expectedInts {{1, 4}, {2, 5}, {3, 6}, {12, 2}};
  EXPECT_EQ(expectedInts, ints);

  std::unordered_map<int, double> doubles;
  found = lr.getDoubleMap("luaArray", doubles);
  EXPECT_TRUE(found);
  std::unordered_map<int, double> expectedDoubles {{1, 4},
                                                   {2, 5},
                                                   {3, 6},
                                                   {12, 2.4}};
  EXPECT_EQ(expectedDoubles, doubles);

  std::unordered_map<int, bool> bools;
  found = lr.getBoolMap("luaArray", bools);
  EXPECT_TRUE(found);
  std::unordered_map<int, bool> expectedBools {{4, true}, {8, false}};
  EXPECT_EQ(expectedBools, bools);

  std::unordered_map<int, std::string> strs;
  found = lr.getStringMap("luaArray", strs);
  EXPECT_TRUE(found);
  std::unordered_map<int, std::string> expectedStrs {{33, "hello"},
                                                     {200, "bye"}};
  EXPECT_EQ(expectedStrs, strs);
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

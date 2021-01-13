// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include <string>
#include <vector>
#include <memory>

#include "axom/config.hpp"

#include "axom/inlet/tests/inlet_test_utils.hpp"

template <typename InletReader>
class inlet_Reader : public testing::Test
{ };

TYPED_TEST_SUITE(inlet_Reader, axom::inlet::detail::ReaderTypes);

using axom::inlet::detail::fromLuaTo;

TYPED_TEST(inlet_Reader, getTopLevelBools)
{
  TypeParam reader;
  reader.parseString(fromLuaTo<TypeParam>("foo = true; bar = false"));

  bool value, retValue;

  value = false;
  retValue = reader.getBool("foo", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, true);

  value = true;
  retValue = reader.getBool("bar", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, false);
}

TYPED_TEST(inlet_Reader, getInsideBools)
{
  TypeParam reader;
  reader.parseString(fromLuaTo<TypeParam>("foo = { bar = false; baz = true }"));

  bool value, retValue;

  value = true;
  retValue = reader.getBool("foo/bar", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, false);

  value = false;
  retValue = reader.getBool("foo/baz", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, true);
}

TYPED_TEST(inlet_Reader, getTopLevelStrings)
{
  TypeParam reader;
  reader.parseString(fromLuaTo<TypeParam>(
    "foo = \"this is a test string\"; bar = \"TesT StrInG\""));

  bool retValue;
  std::string value;

  value = "";
  retValue = reader.getString("foo", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, "this is a test string");

  value = "";
  retValue = reader.getString("bar", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, "TesT StrInG");
}

TYPED_TEST(inlet_Reader, getInsideStrings)
{
  TypeParam reader;
  reader.parseString(fromLuaTo<TypeParam>(
    "foo = { bar = \"this is a test string\"; baz = \"TesT StrInG\" }"));

  bool retValue;
  std::string value;

  value = "";
  retValue = reader.getString("foo/bar", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, "this is a test string");

  value = "";
  retValue = reader.getString("foo/baz", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, "TesT StrInG");
}

TYPED_TEST(inlet_Reader, mixLevelTables)
{
  TypeParam reader;
  reader.parseString(fromLuaTo<TypeParam>(
    "t = { innerT = { foo = 1 }, anotherInnerT = {baz = 3}}"));

  bool retValue;
  int value;

  value = 0;
  retValue = reader.getInt("t/innerT/foo", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, 1);

  value = 0;
  retValue = reader.getInt("t/doesntexist", value);
  EXPECT_EQ(retValue, false);
  EXPECT_EQ(value, 0);

  value = 0;
  retValue = reader.getInt("t/anotherInnerT/baz", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, 3);
}

TYPED_TEST(inlet_Reader, getMap)
{
  // Keep this contiguous in order to test all supported input languages
  std::string testString =
    "luaArray = { [0] = 4, [1] = 5, [2] = 6 , [3] = true, [4] = false, [5] = "
    "2.4, [6] = 'hello', [7] = 'bye' }";
  TypeParam reader;
  reader.parseString(fromLuaTo<TypeParam>(testString));

  std::unordered_map<int, int> ints;
  bool found = reader.getIntMap("luaArray", ints);
  EXPECT_TRUE(found);
  std::unordered_map<int, int> expectedInts {{0, 4}, {1, 5}, {2, 6}, {5, 2}};
  EXPECT_EQ(expectedInts, ints);

  std::unordered_map<int, double> doubles;
  found = reader.getDoubleMap("luaArray", doubles);
  EXPECT_TRUE(found);
  std::unordered_map<int, double> expectedDoubles {{0, 4},
                                                   {1, 5},
                                                   {2, 6},
                                                   {5, 2.4}};
  EXPECT_EQ(expectedDoubles, doubles);

  std::unordered_map<int, bool> bools;
  found = reader.getBoolMap("luaArray", bools);
  EXPECT_TRUE(found);
  std::unordered_map<int, bool> expectedBools {{3, true}, {4, false}};
  EXPECT_EQ(expectedBools, bools);

  // Conduit's YAML parser doesn't distinguish boolean literals from strings
  // so the YAML version will extract the "true" and "false" here
  std::unordered_map<int, std::string> strs;
  found = reader.getStringMap("luaArray", strs);
  EXPECT_TRUE(found);
  // std::unordered_map<int, std::string> expectedStrs {{6, "hello"},
  //                                                    {7, "bye"}};
  // EXPECT_EQ(expectedStrs, strs);
}

TEST(inlet_Reader_YAML, getInsideBools)
{
  axom::inlet::YAMLReader reader;
  bool result = reader.parseString(
    "foo:\n"
    "  bar: false\n"
    "  baz: true");
  EXPECT_TRUE(result);
  bool value, retValue;

  value = true;
  retValue = reader.getBool("foo/bar", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, false);

  value = false;
  retValue = reader.getBool("foo/baz", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, true);
}

TEST(inlet_Reader_YAML, mixLevelTables)
{
  axom::inlet::YAMLReader reader;
  bool result = reader.parseString(
    "t:\n"
    "  innerT:\n"
    "    foo: 1\n"
    "  anotherInnerT:\n"
    "    baz: 3");
  EXPECT_TRUE(result);

  bool retValue;
  int value;

  value = 0;
  retValue = reader.getInt("t/innerT/foo", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, 1);

  value = 0;
  retValue = reader.getInt("t/doesntexist", value);
  EXPECT_EQ(retValue, false);
  EXPECT_EQ(value, 0);

  value = 0;
  retValue = reader.getInt("t/anotherInnerT/baz", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, 3);
}

TEST(inlet_Reader_YAML, mixLevelTables_invalid)
{
  axom::inlet::YAMLReader reader;
  bool result = reader.parseString(
    "t:\n"
    "  innerT: foo: 1\n"
    "  anotherInnerT:\n"
    "    baz: 3");

  EXPECT_FALSE(result);
}

TEST(inlet_Reader_JSON, getInsideBools)
{
  axom::inlet::JSONReader reader;
  bool result = reader.parseString(
    "{\n"
    "  foo: {\n"
    "    bar: false,\n"
    "    baz: true\n"
    "  }\n"
    "}");
  EXPECT_TRUE(result);

  bool value, retValue;

  value = true;
  retValue = reader.getBool("foo/bar", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, false);

  value = false;
  retValue = reader.getBool("foo/baz", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, true);
}

TEST(inlet_Reader_JSON, mixLevelTables)
{
  axom::inlet::JSONReader reader;
  bool result = reader.parseString(
    "{\n"
    "  t: {\n"
    "    innerT: {\n"
    "      foo: 1\n"
    "    },\n"
    "    anotherInnerT: {\n"
    "      baz: 3\n"
    "    }\n"
    "  }\n"
    "}");
  EXPECT_TRUE(result);

  bool retValue;
  int value;

  value = 0;
  retValue = reader.getInt("t/innerT/foo", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, 1);

  value = 0;
  retValue = reader.getInt("t/doesntexist", value);
  EXPECT_EQ(retValue, false);
  EXPECT_EQ(value, 0);

  value = 0;
  retValue = reader.getInt("t/anotherInnerT/baz", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, 3);
}

TEST(inlet_Reader_JSON, mixLevelTables_invalid)
{
  axom::inlet::JSONReader reader;
  bool result = reader.parseString(
    "{\n"
    "  t: {\n"
    "    innerT: {\n"
    "      foo: 1\n"
    "    }\n"
    "    anotherInnerT: {\n"
    "      baz: 3\n"
    "    }\n"
    "  }\n"
    "}");

  EXPECT_FALSE(result);
}

#ifdef AXOM_USE_SOL
// Checks that LuaReader parses array information as expected
// Discontiguous arrays are lua-specific
TEST(inlet_Reader_lua, getDiscontiguousMap)
{
  std::string testString =
    "luaArray = { [1] = 4, [2] = 5, [3] = 6 , [4] = true, [8] = false, [12] = "
    "2.4, [33] = 'hello', [200] = 'bye' }";
  axom::inlet::LuaReader reader;
  reader.parseString(testString);

  std::unordered_map<int, int> ints;
  bool found = reader.getIntMap("luaArray", ints);
  EXPECT_TRUE(found);
  std::unordered_map<int, int> expectedInts {{1, 4}, {2, 5}, {3, 6}, {12, 2}};
  EXPECT_EQ(expectedInts, ints);

  std::unordered_map<int, double> doubles;
  found = reader.getDoubleMap("luaArray", doubles);
  EXPECT_TRUE(found);
  std::unordered_map<int, double> expectedDoubles {{1, 4},
                                                   {2, 5},
                                                   {3, 6},
                                                   {12, 2.4}};
  EXPECT_EQ(expectedDoubles, doubles);

  std::unordered_map<int, bool> bools;
  found = reader.getBoolMap("luaArray", bools);
  EXPECT_TRUE(found);
  std::unordered_map<int, bool> expectedBools {{4, true}, {8, false}};
  EXPECT_EQ(expectedBools, bools);

  std::unordered_map<int, std::string> strs;
  found = reader.getStringMap("luaArray", strs);
  EXPECT_TRUE(found);
  std::unordered_map<int, std::string> expectedStrs {{33, "hello"},
                                                     {200, "bye"}};
  EXPECT_EQ(expectedStrs, strs);
}
#endif

//------------------------------------------------------------------------------
#include "axom/slic/core/SimpleLogger.hpp"
using axom::slic::SimpleLogger;

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  SimpleLogger logger;  // create & initialize test logger,

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}

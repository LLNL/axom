// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include <array>
#include <string>
#include <vector>
#include <unordered_map>

#include <iostream>

#include "axom/sidre.hpp"

#include "axom/inlet/Inlet.hpp"
#include "axom/inlet/tests/inlet_test_utils.hpp"

using axom::inlet::Inlet;
using axom::inlet::InletType;
using axom::inlet::LuaReader;
using axom::sidre::DataStore;

template <typename InletReader>
Inlet createBasicInlet(DataStore* ds,
                       const std::string& luaString,
                       bool enableDocs = true)
{
  auto lr = std::make_unique<InletReader>();
  lr->parseString(axom::inlet::detail::fromLuaTo<InletReader>(luaString));

  return Inlet(std::move(lr), ds->getRoot(), enableDocs);
}

struct Foo
{
  bool bar;
  bool baz;

  bool operator==(const Foo& other) const
  {
    return bar == other.bar && baz == other.baz;
  }
};

template <>
struct FromInlet<Foo>
{
  Foo operator()(const axom::inlet::Table& base)
  {
    Foo f {base["bar"], base["baz"]};
    return f;
  }
};

template <typename InletReader>
class InletObjectTest : public ::testing::Test
{ };

TYPED_TEST_SUITE_P(InletObjectTest);

TYPED_TEST_P(InletObjectTest, simple_struct_by_value)
{
  std::string testString = "foo = { bar = true; baz = false }";
  DataStore ds;
  Inlet inlet = createBasicInlet<TypeParam>(&ds, testString);

  // Define schema
  inlet.addBool("foo/bar", "bar's description");
  inlet.addBool("foo/baz", "baz's description");

  Foo foo;

  foo = inlet.get<Foo>("foo");
  EXPECT_TRUE(foo.bar);
  EXPECT_FALSE(foo.baz);
}

TYPED_TEST_P(InletObjectTest, simple_array_of_struct_by_value)
{
  std::string testString =
    "foo = { [0] = { bar = true; baz = false}, "
    "        [1] = { bar = false; baz = true} }";
  DataStore ds;
  Inlet inlet = createBasicInlet<TypeParam>(&ds, testString);

  auto& arr_table = inlet.addGenericArray("foo");

  arr_table.addBool("bar", "bar's description");
  arr_table.addBool("baz", "baz's description");
  std::unordered_map<int, Foo> expected_foos = {{0, {true, false}},
                                                {1, {false, true}}};
  auto foos = inlet["foo"].get<std::unordered_map<int, Foo>>();
  EXPECT_EQ(foos, expected_foos);
}

TYPED_TEST_P(InletObjectTest, simple_array_of_struct_implicit_idx)
{
  std::string testString =
    "foo = { { bar = true; baz = false}, "
    "        { bar = false; baz = true} }";
  DataStore ds;
  Inlet inlet = createBasicInlet<TypeParam>(&ds, testString);

  auto& arr_table = inlet.addGenericArray("foo");

  arr_table.addBool("bar", "bar's description");
  arr_table.addBool("baz", "baz's description");
  // Lua is 1-indexed
  const int base_idx =
    std::is_same<TypeParam, axom::inlet::LuaReader>::value ? 1 : 0;
  std::unordered_map<int, Foo> expected_foos = {{base_idx, {true, false}},
                                                {base_idx + 1, {false, true}}};
  auto foos = inlet["foo"].get<std::unordered_map<int, Foo>>();
  EXPECT_EQ(foos, expected_foos);
}

TYPED_TEST_P(InletObjectTest, simple_array_of_struct_verify_optional)
{
  std::string testString =
    "foo = { [4] = { bar = true;}, "
    "        [7] = { bar = false;} }";
  DataStore ds;
  Inlet inlet = createBasicInlet<TypeParam>(&ds, testString);

  auto& arr_table = inlet.addGenericArray("foo");

  arr_table.addBool("bar", "bar's description").required(true);
  arr_table.addBool("baz", "baz's description").required(false);

  EXPECT_TRUE(inlet.verify());
}

TYPED_TEST_P(InletObjectTest, simple_array_of_struct_verify_reqd)
{
  std::string testString =
    "foo = { [4] = { bar = true;}, "
    "        [7] = { bar = false;} }";
  DataStore ds;
  Inlet inlet = createBasicInlet<TypeParam>(&ds, testString);

  auto& arr_table = inlet.addGenericArray("foo");

  arr_table.addBool("bar", "bar's description").required(true);
  arr_table.addBool("baz", "baz's description").required(true);

  EXPECT_FALSE(inlet.verify());
}

struct FooWithArray
{
  std::unordered_map<int, int> arr;
  bool operator==(const FooWithArray& other) const { return arr == other.arr; }
};

template <>
struct FromInlet<FooWithArray>
{
  FooWithArray operator()(const axom::inlet::Table& base)
  {
    FooWithArray f = {base["arr"]};
    return f;
  }
};

TYPED_TEST_P(InletObjectTest, array_of_struct_containing_array)
{
  std::string testString =
    "foo = { [0] = { arr = { [0] = 3 }; }, "
    "        [1] = { arr = { [0] = 2 }; } }";
  DataStore ds;
  Inlet inlet = createBasicInlet<TypeParam>(&ds, testString);

  auto& arr_table = inlet.addGenericArray("foo");

  arr_table.addIntArray("arr", "arr's description");
  // Contiguous indexing for generality
  std::unordered_map<int, FooWithArray> expected_foos = {{0, {{{0, 3}}}},
                                                         {1, {{{0, 2}}}}};
  std::unordered_map<int, FooWithArray> foos_with_arr;
  foos_with_arr = inlet["foo"].get<std::unordered_map<int, FooWithArray>>();
  EXPECT_EQ(foos_with_arr, expected_foos);
}

struct MoveOnlyFoo
{
  MoveOnlyFoo() = delete;
  MoveOnlyFoo(const MoveOnlyFoo&) = delete;
  MoveOnlyFoo(MoveOnlyFoo&&) = default;
  MoveOnlyFoo(bool first, bool second) : bar(first), baz(second) { }
  bool bar;
  bool baz;
};

template <>
struct FromInlet<MoveOnlyFoo>
{
  MoveOnlyFoo operator()(const axom::inlet::Table& base)
  {
    MoveOnlyFoo f(base["bar"], base["baz"]);
    return f;
  }
};

TYPED_TEST_P(InletObjectTest, simple_moveonly_struct_by_value)
{
  std::string testString = "foo = { bar = true; baz = false }";
  DataStore ds;
  Inlet inlet = createBasicInlet<TypeParam>(&ds, testString);

  // Define schema
  // Check for existing fields
  inlet.addBool("foo/bar", "bar's description");

  inlet.addBool("foo/baz", "baz's description");

  MoveOnlyFoo foo = inlet.get<MoveOnlyFoo>("foo");
  EXPECT_TRUE(foo.bar);
  EXPECT_FALSE(foo.baz);
}

TYPED_TEST_P(InletObjectTest, simple_value_from_bracket)
{
  std::string testString = "foo = true";
  DataStore ds;
  Inlet inlet = createBasicInlet<TypeParam>(&ds, testString);

  // Define schema
  // Check for existing fields
  inlet.addBool("foo", "foo's description");

  bool foo = inlet["foo"];
  EXPECT_TRUE(foo);
}

TYPED_TEST_P(InletObjectTest, simple_struct_from_bracket)
{
  std::string testString = "foo = { bar = true; baz = false }";
  DataStore ds;
  Inlet inlet = createBasicInlet<TypeParam>(&ds, testString);

  // Define schema
  // Check for existing fields
  inlet.addBool("foo/bar", "bar's description");

  inlet.addBool("foo/baz", "baz's description");

  auto foo = inlet["foo"].get<Foo>();
  EXPECT_TRUE(foo.bar);
  EXPECT_FALSE(foo.baz);
}

TYPED_TEST_P(InletObjectTest, contains_from_table)
{
  std::string testString = "foo = { bar = true; baz = false }";
  DataStore ds;
  Inlet inlet = createBasicInlet<TypeParam>(&ds, testString);

  // Define schema
  // Check for existing fields
  inlet.addBool("foo/bar", "bar's description");

  inlet.addBool("foo/baz", "baz's description");

  EXPECT_TRUE(inlet.contains("foo/bar"));
  EXPECT_TRUE(inlet.contains("foo/baz"));

  auto& foo_table = inlet.getTable("foo");
  EXPECT_TRUE(foo_table.contains("bar"));
  EXPECT_TRUE(foo_table.contains("baz"));
}

TYPED_TEST_P(InletObjectTest, contains_from_bracket)
{
  std::string testString = "foo = { bar = true; baz = false }";
  DataStore ds;
  Inlet inlet = createBasicInlet<TypeParam>(&ds, testString);

  // Define schema
  // Check for existing fields
  inlet.addBool("foo/bar", "bar's description");

  inlet.addBool("foo/baz", "baz's description");

  EXPECT_TRUE(inlet["foo"].contains("bar"));
  EXPECT_TRUE(inlet["foo"].contains("baz"));
}

TYPED_TEST_P(InletObjectTest, array_from_bracket)
{
  DataStore ds;
  std::string testString =
    "luaArrays = { arr1 = { [0] = 4}, "
    "              arr2 = {[0] = true, [1] = false}, "
    "              arr3 = {[0] = 'hello', [1] = 'bye'}, "
    "              arr4 = { [0] = 2.4 } }";
  Inlet inlet = createBasicInlet<TypeParam>(&ds, testString);

  std::unordered_map<int, int> intMap;
  std::unordered_map<int, bool> boolMap;
  std::unordered_map<int, std::string> strMap;
  std::unordered_map<int, double> doubleMap;
  inlet.addIntArray("luaArrays/arr1");
  inlet.addBoolArray("luaArrays/arr2");
  inlet.addStringArray("luaArrays/arr3");
  inlet.addDoubleArray("luaArrays/arr4");

  std::unordered_map<int, int> expectedInts {{0, 4}};
  std::unordered_map<int, bool> expectedBools {{0, true}, {1, false}};
  std::unordered_map<int, std::string> expectedStrs {{0, "hello"}, {1, "bye"}};
  std::unordered_map<int, double> expectedDoubles {{0, 2.4}};

  intMap = inlet["luaArrays/arr1"].get<std::unordered_map<int, int>>();
  EXPECT_EQ(intMap, expectedInts);

  boolMap = inlet["luaArrays/arr2"].get<std::unordered_map<int, bool>>();
  EXPECT_EQ(boolMap, expectedBools);

  strMap = inlet["luaArrays/arr3"].get<std::unordered_map<int, std::string>>();
  EXPECT_EQ(strMap, expectedStrs);

  doubleMap = inlet["luaArrays/arr4"].get<std::unordered_map<int, double>>();
  EXPECT_EQ(doubleMap, expectedDoubles);
}

TYPED_TEST_P(InletObjectTest, primitive_type_checks)
{
  std::string testString =
    " bar = true; baz = 12; quux = 2.5; corge = 'hello' ";
  DataStore ds;
  Inlet inlet = createBasicInlet<TypeParam>(&ds, testString);

  // Define schema
  // Check for existing fields
  inlet.addBool("bar", "bar's description");

  inlet.addInt("baz", "baz's description");

  inlet.addDouble("quux", "quux's description");

  inlet.addString("corge", "corge's description");

  EXPECT_EQ(inlet["bar"].type(), InletType::Bool);
  bool bar = inlet["bar"];
  EXPECT_EQ(bar, true);

  EXPECT_EQ(inlet["baz"].type(), InletType::Integer);
  int baz = inlet["baz"];
  EXPECT_EQ(baz, 12);

  EXPECT_EQ(inlet["quux"].type(), InletType::Double);
  double quux = inlet["quux"];
  EXPECT_DOUBLE_EQ(quux, 2.5);

  EXPECT_EQ(inlet["corge"].type(), InletType::String);
  std::string corge = inlet["corge"];
  EXPECT_EQ(corge, "hello");
}

TYPED_TEST_P(InletObjectTest, composite_type_checks)
{
  std::string testString =
    "luaArrays = { arr1 = { [1] = 4}, "
    "              arr2 = {[4] = true, [8] = false} }; "
    "foo = { bar = true; baz = false }";
  DataStore ds;
  Inlet inlet = createBasicInlet<TypeParam>(&ds, testString);

  // Define schema
  // Check for existing fields
  inlet.addBool("foo/bar", "bar's description");

  inlet.addBool("foo/baz", "baz's description");

  inlet.addIntArray("luaArrays/arr1");
  inlet.addBoolArray("luaArrays/arr2");

  auto arr_table = inlet["luaArrays"];
  // The table containing the two arrays is not an array, but an object
  EXPECT_EQ(arr_table.type(), InletType::Object);

  // But the things it contains are arrays
  EXPECT_EQ(arr_table["arr1"].type(), InletType::Container);
  EXPECT_EQ(arr_table["arr2"].type(), InletType::Container);

  auto foo_table = inlet["foo"];
  // Similarly, the table containing the two bools is an object
  EXPECT_EQ(foo_table.type(), InletType::Object);

  // But the things it contains are booleans
  EXPECT_EQ(foo_table["bar"].type(), InletType::Bool);
  EXPECT_EQ(foo_table["baz"].type(), InletType::Bool);
}

TYPED_TEST_P(InletObjectTest, implicit_conversion_primitives)
{
  std::string testString =
    " bar = true; baz = 12; quux = 2.5; corge = 'hello'; arr = { [0] = 4, [1] "
    "= 6, [2] = 10}";
  DataStore ds;
  Inlet inlet = createBasicInlet<TypeParam>(&ds, testString);

  // Define schema
  inlet.addBool("bar", "bar's description");
  inlet.addInt("baz", "baz's description");
  inlet.addDouble("quux", "quux's description");
  inlet.addString("corge", "corge's description");

  inlet.addIntArray("arr");

  // Attempt both construction and assignment
  bool bar = inlet["bar"];
  EXPECT_EQ(bar, true);
  bar = inlet["bar"];
  EXPECT_EQ(bar, true);

  int baz = inlet["baz"];
  EXPECT_EQ(baz, 12);
  baz = inlet["baz"];
  EXPECT_EQ(baz, 12);

  double quux = inlet["quux"];
  EXPECT_DOUBLE_EQ(quux, 2.5);
  quux = inlet["quux"];
  EXPECT_DOUBLE_EQ(quux, 2.5);

  std::string corge = inlet["corge"];
  EXPECT_EQ(corge, "hello");
  corge = inlet["corge"];
  EXPECT_EQ(corge, "hello");

  std::unordered_map<int, int> expected_arr {{0, 4}, {1, 6}, {2, 10}};
  std::unordered_map<int, int> arr = inlet["arr"];
  EXPECT_EQ(arr, expected_arr);
  arr = inlet["arr"];
  EXPECT_EQ(arr, expected_arr);
}

REGISTER_TYPED_TEST_SUITE_P(InletObjectTest,
                            simple_struct_by_value,
                            simple_array_of_struct_by_value,
                            simple_array_of_struct_implicit_idx,
                            simple_array_of_struct_verify_optional,
                            simple_array_of_struct_verify_reqd,
                            array_of_struct_containing_array,
                            simple_moveonly_struct_by_value,
                            simple_value_from_bracket,
                            simple_struct_from_bracket,
                            contains_from_table,
                            contains_from_bracket,
                            array_from_bracket,
                            primitive_type_checks,
                            composite_type_checks,
                            implicit_conversion_primitives);

INSTANTIATE_TYPED_TEST_SUITE_P(ObjectTests,
                               InletObjectTest,
                               axom::inlet::detail::ReaderTypes);

template <typename InletReader>
class InletDictTest : public ::testing::Test
{ };

TYPED_TEST_SUITE_P(InletDictTest);

TYPED_TEST_P(InletDictTest, basic_dicts)
{
  std::string testString = "foo = { ['key1'] = 4, ['key3'] = 6, ['key2'] = 10}";
  DataStore ds;
  Inlet inlet = createBasicInlet<TypeParam>(&ds, testString);

  inlet.addIntDictionary("foo", "foo's description");
  std::unordered_map<std::string, int> dict = inlet["foo"];
  std::unordered_map<std::string, int> correct = {{"key1", 4},
                                                  {"key3", 6},
                                                  {"key2", 10}};
  EXPECT_EQ(dict, correct);
}

TYPED_TEST_P(InletDictTest, simple_dict_of_struct_by_value)
{
  std::string testString =
    "foo = { ['key1'] = { bar = true; baz = false}, "
    "        ['key2'] = { bar = false; baz = true} }";
  DataStore ds;
  Inlet inlet = createBasicInlet<TypeParam>(&ds, testString);

  auto& dict_table = inlet.addGenericDictionary("foo");

  dict_table.addBool("bar", "bar's description");
  dict_table.addBool("baz", "baz's description");
  std::unordered_map<std::string, Foo> expected_foos = {{"key1", {true, false}},
                                                        {"key2", {false, true}}};
  std::unordered_map<std::string, Foo> foos;
  foos = inlet["foo"].get<std::unordered_map<std::string, Foo>>();
  EXPECT_EQ(foos, expected_foos);
}

struct FooWithDict
{
  std::unordered_map<std::string, int> arr;
  bool operator==(const FooWithDict& other) const { return arr == other.arr; }
};

template <>
struct FromInlet<FooWithDict>
{
  FooWithDict operator()(const axom::inlet::Table& base)
  {
    FooWithDict f = {base["arr"]};
    return f;
  }
};

TYPED_TEST_P(InletDictTest, dict_of_struct_containing_dict)
{
  std::string testString =
    "foo = { ['key3'] = { arr = { ['key1'] = 3 }; }, "
    "        ['key4'] = { arr = { ['key2'] = 2 }; } }";
  DataStore ds;
  Inlet inlet = createBasicInlet<TypeParam>(&ds, testString);

  auto& dict_table = inlet.addGenericDictionary("foo");

  dict_table.addIntDictionary("arr", "arr's description");
  std::unordered_map<std::string, FooWithDict> expected_foos = {
    {"key3", {{{"key1", 3}}}},
    {"key4", {{{"key2", 2}}}}};
  std::unordered_map<std::string, FooWithDict> foos_with_dict;
  foos_with_dict =
    inlet["foo"].get<std::unordered_map<std::string, FooWithDict>>();
  EXPECT_EQ(foos_with_dict, expected_foos);
}

TYPED_TEST_P(InletDictTest, dict_of_struct_containing_array)
{
  std::string testString =
    "foo = { ['key3'] = { arr = { [0] = 3 }; }, "
    "        ['key4'] = { arr = { [0] = 2 }; } }";
  DataStore ds;
  Inlet inlet = createBasicInlet<TypeParam>(&ds, testString);

  auto& dict_table = inlet.addGenericDictionary("foo");

  dict_table.addIntArray("arr", "arr's description");
  std::unordered_map<std::string, FooWithArray> expected_foos = {
    {"key3", {{{0, 3}}}},
    {"key4", {{{0, 2}}}}};
  std::unordered_map<std::string, FooWithArray> foos_with_array;
  foos_with_array =
    inlet["foo"].get<std::unordered_map<std::string, FooWithArray>>();
  EXPECT_EQ(foos_with_array, expected_foos);
}

TYPED_TEST_P(InletDictTest, array_of_struct_containing_dict)
{
  std::string testString =
    "foo = { [0] = { arr = { ['key1'] = 3 }; }, "
    "        [1] = { arr = { ['key2'] = 2 }; } }";
  DataStore ds;
  Inlet inlet = createBasicInlet<TypeParam>(&ds, testString);

  auto& arr_table = inlet.addGenericArray("foo");

  arr_table.addIntDictionary("arr", "arr's description");
  std::unordered_map<int, FooWithDict> expected_foos = {{0, {{{"key1", 3}}}},
                                                        {1, {{{"key2", 2}}}}};
  std::unordered_map<int, FooWithDict> foos_with_dict;
  foos_with_dict = inlet["foo"].get<std::unordered_map<int, FooWithDict>>();
  EXPECT_EQ(foos_with_dict, expected_foos);
}

/*
FIXME: These are currently error conditions.  If these should be supported
or handled differently these tests can be re-enabled.

TEST(inlet_dict, mixed_keys)
{
  std::string testString = "foo = { ['key1'] = 4, [1] = 6 }";
  DataStore ds;
  Inlet inlet = createBasicInlet(&ds, testString);

  inlet.addIntDictionary("foo", "foo's description");
  std::unordered_map<std::string, int> dict = inlet["foo"];
  std::unordered_map<std::string, int> correct_dict = {{"key1", 4}};
  EXPECT_EQ(dict, correct_dict);

  inlet.addIntArray("foo", "foo's description");
  std::unordered_map<int, int> array = inlet["foo"];
  std::unordered_map<int, int> correct_array = {{1, 6}};
  EXPECT_EQ(array, correct_array);
}

TEST(inlet_dict, key_with_slash)
{
  std::string testString =
    "foo = { ['key1/subkey1'] = 4, ['key3'] = 6, ['key2'] = 10}";
  DataStore ds;
  Inlet inlet = createBasicInlet(&ds, testString);

  inlet.addIntDictionary("foo", "foo's description");
  std::unordered_map<std::string, int> dict = inlet["foo"];
  std::unordered_map<std::string, int> correct = {{"key1/subkey1", 4},
                                                  {"key3", 6},
                                                  {"key2", 10}};
  EXPECT_EQ(dict, correct);
}
*/

REGISTER_TYPED_TEST_SUITE_P(InletDictTest,
                            basic_dicts,
                            simple_dict_of_struct_by_value,
                            dict_of_struct_containing_dict,
                            dict_of_struct_containing_array,
                            array_of_struct_containing_dict);

INSTANTIATE_TYPED_TEST_SUITE_P(DictTests,
                               InletDictTest,
                               axom::inlet::detail::ReaderTypes);

// Noncontiguous array tests that are only valid in Lua
#ifdef AXOM_USE_SOL

TEST(LuaInletObjectTest, array_of_struct_containing_array)
{
  std::string testString =
    "foo = { [4] = { arr = { [1] = 3 }; }, "
    "        [7] = { arr = { [6] = 2 }; } }";
  DataStore ds;
  Inlet inlet = createBasicInlet<axom::inlet::LuaReader>(&ds, testString);

  auto& arr_table = inlet.addGenericArray("foo");

  arr_table.addIntArray("arr", "arr's description");
  std::unordered_map<int, FooWithArray> expected_foos = {{4, {{{1, 3}}}},
                                                         {7, {{{6, 2}}}}};
  std::unordered_map<int, FooWithArray> foos_with_arr;
  foos_with_arr = inlet["foo"].get<std::unordered_map<int, FooWithArray>>();
  EXPECT_EQ(foos_with_arr, expected_foos);
}

TEST(LuaInletObjectTest, array_from_bracket)
{
  DataStore ds;
  std::string testString =
    "luaArrays = { arr1 = { [1] = 4}, "
    "              arr2 = {[4] = true, [8] = false}, "
    "              arr3 = {[33] = 'hello', [2] = 'bye'}, "
    "              arr4 = { [12] = 2.4 } }";
  Inlet inlet = createBasicInlet<axom::inlet::LuaReader>(&ds, testString);

  std::unordered_map<int, int> intMap;
  std::unordered_map<int, bool> boolMap;
  std::unordered_map<int, std::string> strMap;
  std::unordered_map<int, double> doubleMap;
  inlet.addIntArray("luaArrays/arr1");
  inlet.addBoolArray("luaArrays/arr2");
  inlet.addStringArray("luaArrays/arr3");
  inlet.addDoubleArray("luaArrays/arr4");

  std::unordered_map<int, int> expectedInts {{1, 4}};
  std::unordered_map<int, bool> expectedBools {{4, true}, {8, false}};
  std::unordered_map<int, std::string> expectedStrs {{33, "hello"}, {2, "bye"}};
  std::unordered_map<int, double> expectedDoubles {{12, 2.4}};

  intMap = inlet["luaArrays/arr1"].get<std::unordered_map<int, int>>();
  EXPECT_EQ(intMap, expectedInts);

  boolMap = inlet["luaArrays/arr2"].get<std::unordered_map<int, bool>>();
  EXPECT_EQ(boolMap, expectedBools);

  strMap = inlet["luaArrays/arr3"].get<std::unordered_map<int, std::string>>();
  EXPECT_EQ(strMap, expectedStrs);

  doubleMap = inlet["luaArrays/arr4"].get<std::unordered_map<int, double>>();
  EXPECT_EQ(doubleMap, expectedDoubles);
}

TEST(LuaInletDictTest, dict_of_struct_containing_array)
{
  std::string testString =
    "foo = { ['key3'] = { arr = { [1] = 3 }; }, "
    "        ['key4'] = { arr = { [6] = 2 }; } }";
  DataStore ds;
  Inlet inlet = createBasicInlet<axom::inlet::LuaReader>(&ds, testString);

  auto& dict_table = inlet.addGenericDictionary("foo");

  dict_table.addIntArray("arr", "arr's description");
  std::unordered_map<std::string, FooWithArray> expected_foos = {
    {"key3", {{{1, 3}}}},
    {"key4", {{{6, 2}}}}};
  std::unordered_map<std::string, FooWithArray> foos_with_array;
  foos_with_array =
    inlet["foo"].get<std::unordered_map<std::string, FooWithArray>>();
  EXPECT_EQ(foos_with_array, expected_foos);
}

TEST(LuaInletDictTest, array_of_struct_containing_dict)
{
  std::string testString =
    "foo = { [7] = { arr = { ['key1'] = 3 }; }, "
    "        [4] = { arr = { ['key2'] = 2 }; } }";
  DataStore ds;
  Inlet inlet = createBasicInlet<axom::inlet::LuaReader>(&ds, testString);

  auto& arr_table = inlet.addGenericArray("foo");

  arr_table.addIntDictionary("arr", "arr's description");
  std::unordered_map<int, FooWithDict> expected_foos = {{7, {{{"key1", 3}}}},
                                                        {4, {{{"key2", 2}}}}};
  std::unordered_map<int, FooWithDict> foos_with_dict;
  foos_with_dict = inlet["foo"].get<std::unordered_map<int, FooWithDict>>();
  EXPECT_EQ(foos_with_dict, expected_foos);
}

#endif

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

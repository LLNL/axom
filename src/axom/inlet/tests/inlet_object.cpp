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
using axom::inlet::VariantKey;
using axom::sidre::DataStore;

template <typename InletReader>
Inlet createBasicInlet(DataStore* ds,
                       const std::string& luaString,
                       bool enableDocs = true)
{
  std::unique_ptr<InletReader> reader(new InletReader());
  reader->parseString(axom::inlet::detail::fromLuaTo<InletReader>(luaString));

  return Inlet(std::move(reader), ds->getRoot(), enableDocs);
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
class inlet_object : public ::testing::Test
{ };

TYPED_TEST_SUITE(inlet_object, axom::inlet::detail::ReaderTypes);

TYPED_TEST(inlet_object, simple_struct_by_value)
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

TYPED_TEST(inlet_object, simple_array_of_struct_by_value)
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

TYPED_TEST(inlet_object, simple_array_of_struct_implicit_idx)
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
  const int base_idx = TypeParam::baseIndex;
  std::unordered_map<int, Foo> expected_foos = {{base_idx, {true, false}},
                                                {base_idx + 1, {false, true}}};
  auto foos = inlet["foo"].get<std::unordered_map<int, Foo>>();
  EXPECT_EQ(foos, expected_foos);
}

TYPED_TEST(inlet_object, simple_array_of_struct_verify_optional)
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

TYPED_TEST(inlet_object, simple_array_of_struct_verify_reqd)
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

TYPED_TEST(inlet_object, simple_array_of_struct_verify_empty)
{
  std::string testString = "foo = { }";
  DataStore ds;
  Inlet inlet = createBasicInlet<TypeParam>(&ds, testString);

  auto& arr_table = inlet.addGenericArray("foo");
  // Even though these are required, the input file should still
  // "verify" as the array is empty
  arr_table.addBool("bar", "bar's description").required(true);
  arr_table.addBool("baz", "baz's description").required(true);

  EXPECT_TRUE(inlet.verify());
}

TYPED_TEST(inlet_object, simple_array_of_struct_verify_lambda)
{
  std::string testString =
    "foo = { [4] = { bar = true;}, "
    "        [7] = { bar = false;} }";
  DataStore ds;
  Inlet inlet = createBasicInlet<TypeParam>(&ds, testString);

  auto& arr_table = inlet.addGenericArray("foo");

  arr_table.addBool("bar", "bar's description").required();

  // Ensuring that foo is the element of the container
  arr_table.registerVerifier(
    [](const axom::inlet::Table& foo) { return foo.contains("bar"); });

  EXPECT_TRUE(inlet.verify());
}

TYPED_TEST(inlet_object, simple_array_of_struct_verify_lambda_pass)
{
  std::string testString =
    "foo = { [4] = { bar = true;}, "
    "        [7] = { baz = false;} }";
  DataStore ds;
  Inlet inlet = createBasicInlet<TypeParam>(&ds, testString);

  auto& arr_table = inlet.addGenericArray("foo");

  arr_table.addBool("bar", "bar's description");
  arr_table.addBool("baz", "baz's description");

  // Check for mutual exclusivity
  arr_table.registerVerifier([](const axom::inlet::Table& foo) {
    return !(foo.contains("bar") && foo.contains("baz"));
  });

  EXPECT_TRUE(inlet.verify());
}

TYPED_TEST(inlet_object, simple_array_of_struct_verify_lambda_fail)
{
  std::string testString =
    "foo = { [4] = { bar = true;}, "
    "        [7] = { bar = false; baz = true} }";
  DataStore ds;
  Inlet inlet = createBasicInlet<TypeParam>(&ds, testString);

  auto& arr_table = inlet.addGenericArray("foo");

  arr_table.addBool("bar", "bar's description");
  arr_table.addBool("baz", "baz's description");

  // Check for mutual exclusivity
  arr_table.registerVerifier([](const axom::inlet::Table& foo) {
    return !(foo.contains("bar") && foo.contains("baz"));
  });

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

TYPED_TEST(inlet_object, array_of_struct_containing_array)
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

TYPED_TEST(inlet_object, simple_moveonly_struct_by_value)
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

TYPED_TEST(inlet_object, simple_value_from_bracket)
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

TYPED_TEST(inlet_object, simple_struct_from_bracket)
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

TYPED_TEST(inlet_object, contains_from_table)
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

TYPED_TEST(inlet_object, contains_from_bracket)
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

TYPED_TEST(inlet_object, array_from_bracket)
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

TYPED_TEST(inlet_object, primitive_type_checks)
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

TYPED_TEST(inlet_object, composite_type_checks)
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

TYPED_TEST(inlet_object, implicit_conversion_primitives)
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

template <typename InletReader>
class inlet_object_dict : public ::testing::Test
{ };

TYPED_TEST_SUITE(inlet_object_dict, axom::inlet::detail::ReaderTypes);

TYPED_TEST(inlet_object_dict, basic_dicts)
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

TYPED_TEST(inlet_object_dict, simple_dict_of_struct_by_value)
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

TYPED_TEST(inlet_object_dict, dict_of_struct_containing_dict)
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

TYPED_TEST(inlet_object_dict, dict_of_struct_containing_array)
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

TYPED_TEST(inlet_object_dict, array_of_struct_containing_dict)
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

TEST(inlet_dict, mixed_keys_primitive_duplicated)
{
  std::string testString = "foo = { ['1'] = 4, [1] = 6 }";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  inlet.addIntDictionary("foo", "foo's description");
  std::unordered_map<VariantKey, int> dict = inlet["foo"];
  std::unordered_map<VariantKey, int> correct_dict = {{"1", 4}, {1, 6}};
  EXPECT_EQ(dict, correct_dict);
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

// Noncontiguous array tests that are only valid in Lua
#ifdef AXOM_USE_SOL

TEST(inlet_object_lua, array_of_struct_containing_array)
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

TEST(inlet_object_lua, array_from_bracket)
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

TEST(inlet_object_lua_dict, dict_of_struct_containing_array)
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

TEST(inlet_object_lua_dict, array_of_struct_containing_dict)
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

TEST(inlet_object_lua_dict, mixed_keys_primitive)
{
  std::string testString = "foo = { ['key1'] = 4, [1] = 6 }";
  DataStore ds;
  Inlet inlet = createBasicInlet<axom::inlet::LuaReader>(&ds, testString);

  inlet.addIntDictionary("foo", "foo's description");
  std::unordered_map<VariantKey, int> dict = inlet["foo"];
  std::unordered_map<VariantKey, int> correct_dict = {{"key1", 4}, {1, 6}};
  EXPECT_EQ(dict, correct_dict);
}

TEST(inlet_object_lua_dict, mixed_keys_primitive_ignore_string_only)
{
  std::string testString = "foo = { ['key1'] = 4, [1] = 6 }";
  DataStore ds;
  Inlet inlet = createBasicInlet<axom::inlet::LuaReader>(&ds, testString);

  inlet.addIntDictionary("foo", "foo's description");
  std::unordered_map<std::string, int> dict = inlet["foo"];
  std::unordered_map<std::string, int> correct_dict = {{"key1", 4}};
  EXPECT_EQ(dict, correct_dict);
}

TEST(inlet_object_lua_dict, mixed_keys_primitive_ignore_int_only)
{
  std::string testString = "foo = { ['key1'] = 4, [1] = 6 }";
  DataStore ds;
  Inlet inlet = createBasicInlet<axom::inlet::LuaReader>(&ds, testString);

  inlet.addIntArray("foo", "foo's description");
  std::unordered_map<int, int> array = inlet["foo"];
  std::unordered_map<int, int> correct_array = {{1, 6}};
  EXPECT_EQ(array, correct_array);
}

TEST(inlet_object_lua_dict, mixed_keys_object)
{
  std::string testString =
    "foo = { ['key1'] = { bar = true; baz = false}, "
    "        [1] = { bar = false; baz = true} }";
  DataStore ds;
  Inlet inlet = createBasicInlet<axom::inlet::LuaReader>(&ds, testString);

  auto& dict_table = inlet.addGenericDictionary("foo");

  dict_table.addBool("bar", "bar's description");
  dict_table.addBool("baz", "baz's description");
  std::unordered_map<VariantKey, Foo> expected_foos = {{"key1", {true, false}},
                                                       {1, {false, true}}};
  std::unordered_map<VariantKey, Foo> foos;
  foos = inlet["foo"].get<std::unordered_map<VariantKey, Foo>>();
  EXPECT_EQ(foos, expected_foos);
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

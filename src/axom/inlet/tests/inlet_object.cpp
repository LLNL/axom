// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include <string>
#include <vector>
#include <unordered_map>

#include <iostream>

#include "axom/sidre.hpp"

#include "axom/inlet/LuaReader.hpp"
#include "axom/inlet/Inlet.hpp"

using axom::inlet::Inlet;
using axom::inlet::LuaReader;
using axom::sidre::DataStore;

std::shared_ptr<Inlet> createBasicInlet(DataStore* ds,
                                        const std::string& luaString,
                                        bool enableDocs = true)
{
  auto lr = std::make_shared<LuaReader>();
  lr->parseString(luaString);

  return std::make_shared<Inlet>(lr, ds->getRoot(), enableDocs);
}

struct Foo
{
  bool bar;
  bool baz;
};

void from_inlet(axom::inlet::Table& base, Foo& f)
{
  f.bar = base["bar"];
  f.baz = base["baz"];
}

TEST(inlet_object, simple_struct_by_ref)
{
  std::string testString = "foo = { bar = true; baz = false }";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  // Define schema
  std::shared_ptr<axom::inlet::Field> currField;

  // Check for existing fields
  currField = inlet->addBool("foo/bar", "bar's description");
  EXPECT_TRUE(currField);

  currField = inlet->addBool("foo/baz", "baz's description");
  EXPECT_TRUE(currField);

  Foo foo;
  bool found = inlet->get("foo", foo);

  EXPECT_TRUE(found);
  EXPECT_TRUE(foo.bar);
  EXPECT_FALSE(foo.baz);
}

template <>
Foo from_inlet<Foo>(axom::inlet::Table& base)
{
  Foo f {base["bar"], base["baz"]};
  return f;
}

TEST(inlet_object, simple_struct_by_value)
{
  std::string testString = "foo = { bar = true; baz = false }";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  // Define schema
  std::shared_ptr<axom::inlet::Field> currField;

  // Check for existing fields
  currField = inlet->addBool("foo/bar", "bar's description");
  EXPECT_TRUE(currField);

  currField = inlet->addBool("foo/baz", "baz's description");
  EXPECT_TRUE(currField);

  Foo foo;

  foo = inlet->get<Foo>("foo");
  EXPECT_TRUE(foo.bar);
  EXPECT_FALSE(foo.baz);
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
MoveOnlyFoo from_inlet<MoveOnlyFoo>(axom::inlet::Table& base)
{
  MoveOnlyFoo f(base["bar"], base["baz"]);
  return f;
}

TEST(inlet_object, simple_moveonly_struct_by_value)
{
  std::string testString = "foo = { bar = true; baz = false }";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  // Define schema
  std::shared_ptr<axom::inlet::Field> currField;

  // Check for existing fields
  currField = inlet->addBool("foo/bar", "bar's description");
  EXPECT_TRUE(currField);

  currField = inlet->addBool("foo/baz", "baz's description");
  EXPECT_TRUE(currField);

  MoveOnlyFoo foo = inlet->get<MoveOnlyFoo>("foo");
  EXPECT_TRUE(foo.bar);
  EXPECT_FALSE(foo.baz);
}

TEST(inlet_object, simple_value_from_bracket)
{
  std::string testString = "foo = true";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  // Define schema
  std::shared_ptr<axom::inlet::Field> currField;

  // Check for existing fields
  currField = inlet->addBool("foo", "foo's description");
  EXPECT_TRUE(currField);

  bool foo = (*inlet)["foo"];
  EXPECT_TRUE(foo);
}

TEST(inlet_object, simple_struct_from_bracket)
{
  std::string testString = "foo = { bar = true; baz = false }";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  // Define schema
  std::shared_ptr<axom::inlet::Field> currField;

  // Check for existing fields
  currField = inlet->addBool("foo/bar", "bar's description");
  EXPECT_TRUE(currField);

  currField = inlet->addBool("foo/baz", "baz's description");
  EXPECT_TRUE(currField);

  Foo foo = (*inlet)["foo"];
  EXPECT_TRUE(foo.bar);
  EXPECT_FALSE(foo.baz);
}

TEST(inlet_object, array_from_bracket)
{
  DataStore ds;
  std::string testString =
    "luaArrays = { arr1 = { [1] = 4}, "
    "              arr2 = {[4] = true, [8] = false}, "
    "              arr3 = {[33] = 'hello', [2] = 'bye'}, "
    "              arr4 = { [12] = 2.4 } }";
  auto inlet = createBasicInlet(&ds, testString);

  std::unordered_map<int, int> intMap;
  std::unordered_map<int, bool> boolMap;
  std::unordered_map<int, std::string> strMap;
  std::unordered_map<int, double> doubleMap;
  auto arr1_field = inlet->getGlobalTable()->addIntArray("luaArrays/arr1");
  auto arr2_field = inlet->getGlobalTable()->addBoolArray("luaArrays/arr2");
  auto arr3_field = inlet->getGlobalTable()->addStringArray("luaArrays/arr3");
  auto arr4_field = inlet->getGlobalTable()->addDoubleArray("luaArrays/arr4");

  std::unordered_map<int, int> expectedInts {{1, 4}};
  std::unordered_map<int, bool> expectedBools {{4, true}, {8, false}};
  std::unordered_map<int, std::string> expectedStrs {{33, "hello"}, {2, "bye"}};
  std::unordered_map<int, double> expectedDoubles {{12, 2.4}};

  intMap = (*inlet)["luaArrays/arr1"].get<std::unordered_map<int, int>>();
  EXPECT_EQ(intMap, expectedInts);

  boolMap = (*inlet)["luaArrays/arr2"].get<std::unordered_map<int, bool>>();
  EXPECT_EQ(boolMap, expectedBools);

  strMap = (*inlet)["luaArrays/arr3"].get<std::unordered_map<int, std::string>>();
  EXPECT_EQ(strMap, expectedStrs);

  doubleMap = (*inlet)["luaArrays/arr4"].get<std::unordered_map<int, double>>();
  EXPECT_EQ(doubleMap, expectedDoubles);
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

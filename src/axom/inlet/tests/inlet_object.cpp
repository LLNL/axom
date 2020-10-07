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

bool from_inlet(axom::inlet::Table& base, Foo& f)
{
  if(!base.get("bar", f.bar))
  {
    return false;
  }

  if(!base.get("baz", f.baz))
  {
    return false;
  }

  return true;
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
  Foo f;
  base.get("bar", f.bar);
  base.get("baz", f.baz);
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
  bool bar;
  bool baz;
  base.get("bar", bar);
  base.get("baz", baz);
  MoveOnlyFoo f {bar, baz};
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

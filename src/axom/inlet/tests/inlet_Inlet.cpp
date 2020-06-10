// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include <string>
#include <vector>

#include <iostream>

#include "axom/sidre.hpp"

#include "axom/inlet/LuaReader.hpp"
#include "axom/inlet/Inlet.hpp"

TEST(inlet_Inlet_basic, getTopLevelBools)
{
  axom::inlet::LuaReader lr;
  lr.parseString("foo = true; bar = false");

  axom::inlet::Inlet inlet;
  inlet.reader(&lr);
  axom::sidre::DataStore ds;
  inlet.sidreGroup(ds.getRoot());

  axom::sidre::Group* group = nullptr;

  //
  // Define schema
  //

  // Check for existing fields
  group = inlet.addBool("foo", "foo's description");
  EXPECT_TRUE(group != nullptr);

  group = inlet.addBool("bar", "bar's description");
  EXPECT_TRUE(group != nullptr);

  // Check one that doesn't exist and doesn't have a default value
  group = inlet.addBool("nonexistant", "nothing");
  EXPECT_FALSE(group == nullptr);

  //
  // Check stored values from get
  //

  bool value = false;
  bool found = false;

  // Check for existing fields
  found = inlet.get("foo", value);
  EXPECT_TRUE(found);
  EXPECT_TRUE(value);

  found = inlet.get("bar", value);
  EXPECT_TRUE(found);
  EXPECT_FALSE(value);

  // Check one that doesn't exist and doesn't have a default value
  value = true;
  found = inlet.get("nonexistant", value);
  EXPECT_FALSE(found);
  EXPECT_TRUE(value);
}

TEST(inlet_Inlet_basic, getNestedBools)
{
  axom::inlet::LuaReader lr;
  lr.parseString("foo = { bar = true; baz = false }");

  axom::inlet::Inlet inlet;
  inlet.reader(&lr);
  axom::sidre::DataStore ds;
  inlet.sidreGroup(ds.getRoot());

  axom::sidre::Group* group = nullptr;

  //
  // Define schema
  //

  // Check for existing fields
  group = inlet.addBool("foo/bar", "bar's description");
  EXPECT_TRUE(group != nullptr);

  group = inlet.addBool("foo/baz", "baz's description");
  EXPECT_TRUE(group != nullptr);

  // Check one that doesn't exist and doesn't have a default value
  group = inlet.addBool("foo/nonexistant", "nothing");
  EXPECT_FALSE(group == nullptr);

  //
  // Check stored values from get
  //

  bool value = false;
  bool found = false;

  // Check for existing fields
  found = inlet.get("foo/bar", value);
  EXPECT_TRUE(found);
  EXPECT_TRUE(value);

  found = inlet.get("foo/baz", value);
  EXPECT_TRUE(found);
  EXPECT_FALSE(value);

  // Check one that doesn't exist and doesn't have a default value
  value = true;
  found = inlet.get("foo/nonexistant", value);
  EXPECT_FALSE(found);
  EXPECT_TRUE(value);
}

TEST(inlet_Inlet_basic, getTopLevelDoubles)
{
  axom::inlet::LuaReader lr;
  lr.parseString("foo = 5.05; bar = 15.1");

  axom::inlet::Inlet inlet;
  inlet.reader(&lr);
  axom::sidre::DataStore ds;
  inlet.sidreGroup(ds.getRoot());

  axom::sidre::Group* group = nullptr;

  //
  // Define schema
  //

  // Check for existing fields
  group = inlet.addDouble("foo", "foo's description");
  EXPECT_TRUE(group != nullptr);

  group = inlet.addDouble("bar", "bar's description");
  EXPECT_TRUE(group != nullptr);

  // Check one that doesn't exist and doesn't have a default value
  group = inlet.addDouble("nonexistant", "nothing");
  EXPECT_FALSE(group == nullptr);

  //
  // Check stored values from get
  //

  double value = -1;
  bool found = false;

  // Check for existing fields
  found = inlet.get("foo", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, 5.05);

  found = inlet.get("bar", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, 15.1);

  // Check one that doesn't exist and doesn't have a default value
  value = 1000;
  found = inlet.get("nonexistant", value);
  EXPECT_FALSE(found);
  EXPECT_EQ(value, 1000);
}

TEST(inlet_Inlet_basic, getNestedDoubles)
{
  axom::inlet::LuaReader lr;
  lr.parseString("foo = { bar = 200.5; baz = 100.987654321 }");

  axom::inlet::Inlet inlet;
  inlet.reader(&lr);
  axom::sidre::DataStore ds;
  inlet.sidreGroup(ds.getRoot());

  axom::sidre::Group* group = nullptr;

  //
  // Define schema
  //

  // Check for existing fields
  group = inlet.addDouble("foo/bar", "bar's description");
  EXPECT_TRUE(group != nullptr);

  group = inlet.addDouble("foo/baz", "baz's description");
  EXPECT_TRUE(group != nullptr);

  // Check one that doesn't exist and doesn't have a default value
  group = inlet.addDouble("foo/nonexistant", "nothing");
  EXPECT_FALSE(group == nullptr);

  //
  // Check stored values from get
  //

  double value = -1;
  bool found = false;

  // Check for existing fields
  found = inlet.get("foo/bar", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, 200.5);

  found = inlet.get("foo/baz", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, 100.987654321);

  // Check one that doesn't exist and doesn't have a default value
  value = 1000;
  found = inlet.get("foo/nonexistant", value);
  EXPECT_FALSE(found);
  EXPECT_EQ(value, 1000);
}

TEST(inlet_Inlet_basic, getTopLevelInts)
{
  axom::inlet::LuaReader lr;
  lr.parseString("foo = 5; bar = 15");

  axom::inlet::Inlet inlet;
  inlet.reader(&lr);
  axom::sidre::DataStore ds;
  inlet.sidreGroup(ds.getRoot());

  axom::sidre::Group* group = nullptr;

  //
  // Define schema
  //

  // Check for existing fields
  group = inlet.addInt("foo", "foo's description");
  EXPECT_TRUE(group != nullptr);

  group = inlet.addInt("bar", "bar's description");
  EXPECT_TRUE(group != nullptr);

  // Check one that doesn't exist and doesn't have a default value
  group = inlet.addInt("nonexistant", "nothing");
  EXPECT_FALSE(group == nullptr);

  //
  // Check stored values from get
  //

  int value = -1;
  bool found = false;

  // Check for existing fields
  found = inlet.get("foo", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, 5);

  found = inlet.get("bar", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, 15);

  // Check one that doesn't exist and doesn't have a default value
  value = 1000;
  found = inlet.get("nonexistant", value);
  EXPECT_FALSE(found);
  EXPECT_EQ(value, 1000);
}

TEST(inlet_Inlet_basic, getNestedInts)
{
  axom::inlet::LuaReader lr;
  lr.parseString("foo = { bar = 200; baz = 100 }");

  axom::inlet::Inlet inlet;
  inlet.reader(&lr);
  axom::sidre::DataStore ds;
  inlet.sidreGroup(ds.getRoot());

  axom::sidre::Group* group = nullptr;

  //
  // Define schema
  //

  // Check for existing fields
  group = inlet.addInt("foo/bar", "bar's description");
  EXPECT_TRUE(group != nullptr);

  group = inlet.addInt("foo/baz", "baz's description");
  EXPECT_TRUE(group != nullptr);

  // Check one that doesn't exist and doesn't have a default value
  group = inlet.addInt("foo/nonexistant", "nothing");
  EXPECT_FALSE(group == nullptr);

  //
  // Check stored values from get
  //

  int value = -1;
  bool found = false;

  // Check for existing fields
  found = inlet.get("foo/bar", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, 200);

  found = inlet.get("foo/baz", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, 100);

  // Check one that doesn't exist and doesn't have a default value
  value = 1000;
  found = inlet.get("foo/nonexistant", value);
  EXPECT_FALSE(found);
  EXPECT_EQ(value, 1000);
}

TEST(inlet_Inlet_basic, getTopLevelStrings)
{
  axom::inlet::LuaReader lr;
  lr.parseString("foo = 'test string'; bar = '15'");

  axom::inlet::Inlet inlet;
  inlet.reader(&lr);
  axom::sidre::DataStore ds;
  inlet.sidreGroup(ds.getRoot());

  axom::sidre::Group* group = nullptr;

  //
  // Define schema
  //

  // Check for existing fields
  group = inlet.addString("foo", "foo's description");
  EXPECT_TRUE(group != nullptr);

  group = inlet.addString("bar", "bar's description");
  EXPECT_TRUE(group != nullptr);

  // Check one that doesn't exist and doesn't have a default value
  group = inlet.addString("nonexistant", "nothing");
  EXPECT_FALSE(group == nullptr);

  //
  // Check stored values from get
  //

  std::string value = "";
  bool found = false;

  // Check for existing fields
  found = inlet.get("foo", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, "test string");

  found = inlet.get("bar", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, "15");

  // Check one that doesn't exist and doesn't have a default value
  value = "don't change";
  found = inlet.get("nonexistant", value);
  EXPECT_FALSE(found);
  EXPECT_EQ(value, "don't change");
}

TEST(inlet_Inlet_basic, getNestedStrings)
{
  axom::inlet::LuaReader lr;
  lr.parseString("foo = { bar = 'yet another string'; baz = '' }");

  axom::inlet::Inlet inlet;
  inlet.reader(&lr);
  axom::sidre::DataStore ds;
  inlet.sidreGroup(ds.getRoot());

  axom::sidre::Group* group = nullptr;

  //
  // Define schema
  //

  // Check for existing fields
  group = inlet.addString("foo/bar", "bar's description");
  EXPECT_TRUE(group != nullptr);

  group = inlet.addString("foo/baz", "baz's description");
  EXPECT_TRUE(group != nullptr);

  // Check one that doesn't exist and doesn't have a default value
  group = inlet.addString("foo/nonexistant", "nothing");
  EXPECT_FALSE(group == nullptr);

  //
  // Check stored values from get
  //

  std::string value = "";
  bool found = false;

  // Check for existing fields
  found = inlet.get("foo/bar", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, "yet another string");

  found = inlet.get("foo/baz", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, "");

  // Check one that doesn't exist and doesn't have a default value
  value = "1000";
  found = inlet.get("foo/nonexistant", value);
  EXPECT_FALSE(found);
  EXPECT_EQ(value, "1000");
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

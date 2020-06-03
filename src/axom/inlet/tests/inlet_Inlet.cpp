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

TEST(inlet_Inlet_basic, getTopLevelInts)
{
  axom::inlet::LuaReader lr;
  lr.parseString("foo = 5; bar = 15");

  axom::inlet::Inlet inlet;
  inlet.reader(&lr);
  axom::sidre::DataStore ds;
  inlet.datastore(&ds);

  bool found = false;

  //
  // Define schema
  //

  // Check for existing fields
  found = inlet.addInt("foo", "foo's description");
  EXPECT_TRUE(found);

  found = inlet.addInt("bar", "bar's description");
  EXPECT_TRUE(found);

  // Check one that doesn't exist but has a default value
  found = inlet.addInt("baz", "baz's description", 100);
  EXPECT_TRUE(found);

  // Check one that doesn't exist and doesn't have a default value
  found = inlet.addInt("nonexistant", "nothing");
  EXPECT_FALSE(found);

  //
  // Check stored values from get
  //

  int value = -1;

  // Check for existing fields
  found = inlet.get("foo", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, 5);

  found = inlet.get("bar", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, 15);

  // Check one that doesn't exist but has a default value
  found = inlet.get("baz", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, 100);

  // Check one that doesn't exist and doesn't have a default value
  found = inlet.get("nonexistant", value);
  EXPECT_FALSE(found);
}

TEST(inlet_Inlet_basic, getNestedInts)
{
  axom::inlet::LuaReader lr;
  lr.parseString("foo = { bar = 200; baz = 100 }");

  axom::inlet::Inlet inlet;
  inlet.reader(&lr);
  axom::sidre::DataStore ds;
  inlet.datastore(&ds);

  bool found = false;

  //
  // Define schema
  //

  // Check for existing fields
  found = inlet.addInt("foo/bar", "bar's description");
  EXPECT_TRUE(found);

  found = inlet.addInt("foo/baz", "baz's description");
  EXPECT_TRUE(found);

  // Check one that doesn't exist but has a default value
  found = inlet.addInt("foo/default", "defaults's description", 900);
  EXPECT_TRUE(found);

  // Check one that doesn't exist and doesn't have a default value
  found = inlet.addInt("foo/nonexistant", "nothing");
  EXPECT_FALSE(found);

  //
  // Check stored values from get
  //

  int value = -1;

  // Check for existing fields
  found = inlet.get("foo/bar", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, 200);

  found = inlet.get("foo/baz", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, 100);

  // Check one that doesn't exist but has a default value
  found = inlet.get("foo/default", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, 900);

  // Check one that doesn't exist and doesn't have a default value
  found = inlet.get("foo/nonexistant", value);
  EXPECT_FALSE(found);
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

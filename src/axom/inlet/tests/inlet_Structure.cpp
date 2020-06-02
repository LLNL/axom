// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include <string>
#include <vector>

#include <iostream>

#include "axom/inlet/LuaMap.hpp"
#include "axom/inlet/MapBackend.hpp"
#include "axom/inlet/Structure.hpp"

TEST(inlet_Structure_getBool, getTopLevelInts)
{
  axom::inlet::LuaMap lm;
  lm.parseString("foo = 5; bar = 15");

  axom::inlet::MapBackend mb;

  axom::inlet::Structure s;
  s.map(&lm);
  s.backend(&mb);

  axom::inlet::IntField* intfield = nullptr;

  //
  // Check return values from add
  //

  // Check for existing fields
  intfield = s.addIntField("foo", "foo's description");
  EXPECT_TRUE(intfield != nullptr);
  EXPECT_EQ(intfield->value(), 5);

  intfield = s.addIntField("bar", "bar's description");
  EXPECT_TRUE(intfield != nullptr);
  EXPECT_EQ(intfield->value(), 15);

  // Check one that doesn't exist but has a default value
  intfield = s.addIntField("baz", "baz's description", 100);
  EXPECT_TRUE(intfield != nullptr);
  EXPECT_EQ(intfield->value(), 100);

  // Check one that doesn't exist and doesn't have a default value
  intfield = s.addIntField("nonexistant", "nothing");
  EXPECT_TRUE(intfield == nullptr);

  //
  // Check stored values from get
  //
  
  // Check for existing fields
  intfield = s.getIntField("foo");
  EXPECT_TRUE(intfield != nullptr);
  EXPECT_EQ(intfield->value(), 5);

  intfield = s.getIntField("bar");
  EXPECT_TRUE(intfield != nullptr);
  EXPECT_EQ(intfield->value(), 15);

  // Check one that doesn't exist but has a default value
  intfield = s.getIntField("baz");
  EXPECT_TRUE(intfield != nullptr);
  EXPECT_EQ(intfield->value(), 100);

  // Check one that doesn't exist and doesn't have a default value
  intfield = s.getIntField("nonexistant");
  EXPECT_TRUE(intfield == nullptr);
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

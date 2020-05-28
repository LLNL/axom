// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include <string>
#include <vector>

#include <iostream>

#include "axom/slim/LuaMap.hpp"
#include "axom/slim/MapBackend.hpp"
#include "axom/slim/Structure.hpp"

TEST(slim_Structure_getBool, getTopLevelInts)
{
  axom::slim::LuaMap lm;
  lm.parseString("foo = 5; bar = 15");

  axom::slim::MapBackend mb;

  axom::slim::Structure s;
  s.map(&lm);
  s.backend(&mb);

  axom::slim::IntField* intfield = nullptr;

  std::cout << "foo" << std::endl;

  // Check for existing fields
  intfield = s.addIntField("foo", "foo's description");
  EXPECT_TRUE(intfield != nullptr);
  EXPECT_EQ(intfield->value(), 5);

  std::cout << "bar" << std::endl;

  intfield = s.addIntField("bar", "bar's description");
  EXPECT_TRUE(intfield != nullptr);
  EXPECT_EQ(intfield->value(), 15);

  std::cout << "baz" << std::endl;

  // Check one that doesn't exist but has a default value
  intfield = s.addIntField("baz", "baz's description", 100);
  EXPECT_TRUE(intfield != nullptr);
  EXPECT_EQ(intfield->value(), 100);

  std::cout << "nonexistant" << std::endl;

  // Check one that doesn't exist and doesn't have a default value
  intfield = s.addIntField("nonexistant", "nothing");
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

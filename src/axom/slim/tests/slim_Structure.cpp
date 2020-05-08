// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include <string>
#include <vector>

#include "axom/slim/LuaMap.hpp"
#include "axom/slim/Structure.hpp"

TEST(slim_Structure_getBool, getTopLevelInts)
{
  axom::slim::LuaMap lm;
  lm.parseString("foo = 5; bar = 15");

  axom::slim::Structure s;
  s.map(&lm);

  axom::slim::IntField* intfield = nullptr;

  intfield = s.addIntField("foo", "foo's description");
  EXPECT_TRUE(intfield != nullptr);
  EXPECT_EQ(intfield->value(), 5);

  intfield = s.addIntField("bar", "bar's description");
  EXPECT_TRUE(intfield != nullptr);
  EXPECT_EQ(intfield->value(), 15);

  intfield = s.addIntField("baz", "baz's description", 100);
  EXPECT_TRUE(intfield != nullptr);
  EXPECT_EQ(intfield->value(), 100);
}

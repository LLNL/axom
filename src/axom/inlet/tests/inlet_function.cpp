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

#include "axom/inlet/LuaReader.hpp"
#include "axom/inlet/Inlet.hpp"

using axom::inlet::Inlet;
using axom::inlet::InletType;
using axom::inlet::LuaReader;
using axom::sidre::DataStore;

Inlet createBasicInlet(DataStore* ds,
                       const std::string& luaString,
                       bool enableDocs = true)
{
  auto lr = std::make_unique<LuaReader>();
  lr->parseString(luaString);

  return Inlet(std::move(lr), ds->getRoot(), enableDocs);
}

TEST(inlet_function, simple_vec_to_double_raw)
{
  std::string testString = "function foo (x, y, z) return x + y + z end";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  auto func = inlet.reader().getFunction("foo",
                                         axom::inlet::InletFunctionType::Double,
                                         axom::inlet::InletFunctionType::Vec3D);

  auto result = func.call<double>(axom::primal::Vector3D {1, 2, 3});
  EXPECT_FLOAT_EQ(result, 6);
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

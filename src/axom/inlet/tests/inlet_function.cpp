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
using axom::inlet::InletFunctionType;
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

TEST(inlet_function, simple_vec3_to_double_raw)
{
  std::string testString = "function foo (x, y, z) return x + y + z end";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  auto func = inlet.reader().getFunction("foo",
                                         InletFunctionType::Double,
                                         InletFunctionType::Vec3D);

  EXPECT_TRUE(func);
  auto result = func.call<double>(axom::primal::Vector3D {1, 2, 3});
  EXPECT_FLOAT_EQ(result, 6);
}

TEST(inlet_function, simple_vec3_to_vec3_raw)
{
  std::string testString = "function foo (x, y, z) return 2*x, 2*y, 2*z end";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  auto func = inlet.reader().getFunction("foo",
                                         InletFunctionType::Vec3D,
                                         InletFunctionType::Vec3D);

  EXPECT_TRUE(func);
  auto result =
    func.call<axom::primal::Vector3D>(axom::primal::Vector3D {1, 2, 3});
  EXPECT_FLOAT_EQ(result[0], 2);
  EXPECT_FLOAT_EQ(result[1], 4);
  EXPECT_FLOAT_EQ(result[2], 6);
}

TEST(inlet_function, simple_vec3_to_double_through_table)
{
  std::string testString = "function foo (x, y, z) return x + y + z end";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  inlet.addFunction("foo",
                    InletFunctionType::Double,
                    InletFunctionType::Vec3D,
                    "foo's description");

  auto callable =
    inlet["foo"].get<std::function<double(axom::primal::Vector3D)>>();
  auto result = callable({1, 2, 3});
  EXPECT_FLOAT_EQ(result, 6);
}

TEST(inlet_function, simple_vec3_to_vec3_through_table)
{
  std::string testString = "function foo (x, y, z) return 2*x, 2*y, 2*z end";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  inlet.addFunction("foo",
                    InletFunctionType::Vec3D,
                    InletFunctionType::Vec3D,
                    "foo's description");

  auto callable =
    inlet["foo"]
      .get<std::function<axom::primal::Vector3D(axom::primal::Vector3D)>>();
  auto result = callable({1, 2, 3});
  EXPECT_FLOAT_EQ(result[0], 2);
  EXPECT_FLOAT_EQ(result[1], 4);
  EXPECT_FLOAT_EQ(result[2], 6);
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

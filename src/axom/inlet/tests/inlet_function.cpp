// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
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

using axom::inlet::FunctionTag;
using axom::inlet::FunctionType;
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

TEST(inlet_function, simple_vec3_to_double_raw)
{
  std::string testString = "function foo (x, y, z) return x + y + z end";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  auto func =
    inlet.reader().getFunction("foo", FunctionTag::Double, {FunctionTag::Vec3D});

  EXPECT_TRUE(func);
  auto result = func.call<double>(FunctionType::Vec3D {1, 2, 3});
  EXPECT_FLOAT_EQ(result, 6);
}

TEST(inlet_function, simple_vec3_to_vec3_raw)
{
  std::string testString = "function foo (x, y, z) return 2*x, 2*y, 2*z end";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  auto func =
    inlet.reader().getFunction("foo", FunctionTag::Vec3D, {FunctionTag::Vec3D});

  EXPECT_TRUE(func);
  auto result = func.call<FunctionType::Vec3D>(FunctionType::Vec3D {1, 2, 3});
  EXPECT_FLOAT_EQ(result[0], 2);
  EXPECT_FLOAT_EQ(result[1], 4);
  EXPECT_FLOAT_EQ(result[2], 6);
}

TEST(inlet_function, simple_vec3_to_vec3_raw_partial_init)
{
  std::string testString = "function foo (x, y, z) return 2*x, 2*y, 2*z end";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  auto func =
    inlet.reader().getFunction("foo", FunctionTag::Vec3D, {FunctionTag::Vec3D});

  EXPECT_TRUE(func);

  auto result = func.call<FunctionType::Vec3D>(FunctionType::Vec3D {1, 2});
  EXPECT_FLOAT_EQ(result[0], 2);
  EXPECT_FLOAT_EQ(result[1], 4);
  EXPECT_FLOAT_EQ(result[2], 0);

  result = func.call<FunctionType::Vec3D>(FunctionType::Vec3D {1});
  EXPECT_FLOAT_EQ(result[0], 2);
  EXPECT_FLOAT_EQ(result[1], 0);
  EXPECT_FLOAT_EQ(result[2], 0);
}

TEST(inlet_function, simple_vec3_to_double_through_table)
{
  std::string testString = "function foo (x, y, z) return x + y + z end";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  inlet.addFunction("foo",
                    FunctionTag::Double,
                    {FunctionTag::Vec3D},
                    "foo's description");

  auto callable = inlet["foo"].get<std::function<double(FunctionType::Vec3D)>>();
  auto result = callable({1, 2, 3});
  EXPECT_FLOAT_EQ(result, 6);
}

TEST(inlet_function, simple_vec3_to_vec3_through_table)
{
  std::string testString = "function foo (x, y, z) return 2*x, 2*y, 2*z end";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  inlet.addFunction("foo",
                    FunctionTag::Vec3D,
                    {FunctionTag::Vec3D},
                    "foo's description");

  auto callable =
    inlet["foo"].get<std::function<FunctionType::Vec3D(FunctionType::Vec3D)>>();
  auto result = callable({1, 2, 3});
  EXPECT_FLOAT_EQ(result[0], 2);
  EXPECT_FLOAT_EQ(result[1], 4);
  EXPECT_FLOAT_EQ(result[2], 6);
}

TEST(inlet_function, simple_double_to_double_through_table)
{
  std::string testString = "function foo (a) return (a * 3.4) + 9.64 end";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  inlet.addFunction("foo",
                    FunctionTag::Double,
                    {FunctionTag::Double},
                    "foo's description");

  auto callable =
    inlet["foo"].get<std::function<FunctionType::Double(FunctionType::Double)>>();
  double arg = -6.37;
  double result = callable(arg);
  EXPECT_FLOAT_EQ(result, (arg * 3.4) + 9.64);
}

TEST(inlet_function, simple_vec3_to_double_through_table_call)
{
  std::string testString = "function foo (x, y, z) return x + y + z end";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  inlet.addFunction("foo",
                    FunctionTag::Double,
                    {FunctionTag::Vec3D},
                    "foo's description");

  auto result = inlet["foo"].call<double>(FunctionType::Vec3D {1, 2, 3});
  EXPECT_FLOAT_EQ(result, 6);
}

TEST(inlet_function, simple_vec3_to_vec3_through_table_call)
{
  std::string testString = "function foo (x, y, z) return 2*x, 2*y, 2*z end";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  inlet.addFunction("foo",
                    FunctionTag::Vec3D,
                    {FunctionTag::Vec3D},
                    "foo's description");

  auto result =
    inlet["foo"].call<FunctionType::Vec3D>(FunctionType::Vec3D {1, 2, 3});
  EXPECT_FLOAT_EQ(result[0], 2);
  EXPECT_FLOAT_EQ(result[1], 4);
  EXPECT_FLOAT_EQ(result[2], 6);
}

TEST(inlet_function, simple_vec3_double_to_double_through_table_call)
{
  std::string testString =
    "function foo (x, y, z, t) return t * (x + y + z) end";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  inlet.addFunction("foo",
                    FunctionTag::Double,
                    {FunctionTag::Vec3D, FunctionTag::Double},
                    "foo's description");

  auto result = inlet["foo"].call<double>(FunctionType::Vec3D {1, 2, 3}, 2.0);
  EXPECT_FLOAT_EQ(result, 12);
}

TEST(inlet_function, simple_vec3_double_to_vec3_through_table_call)
{
  std::string testString = "function foo (x, y, z, t) return t*x, t*y, t*z end";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  inlet.addFunction("foo",
                    FunctionTag::Vec3D,
                    {FunctionTag::Vec3D, FunctionTag::Double},
                    "foo's description");

  auto result =
    inlet["foo"].call<FunctionType::Vec3D>(FunctionType::Vec3D {1, 2, 3}, 2.0);
  EXPECT_FLOAT_EQ(result[0], 2);
  EXPECT_FLOAT_EQ(result[1], 4);
  EXPECT_FLOAT_EQ(result[2], 6);
}

TEST(inlet_function, simple_vec3_to_vec3_verify_lambda_pass)
{
  std::string testString = "function foo (x, y, z) return 2*x, 2*y, 2*z end";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  auto& func = inlet
                 .addFunction("foo",
                              FunctionTag::Vec3D,
                              {FunctionTag::Vec3D},
                              "foo's description")
                 .required();
  func.registerVerifier([](const axom::inlet::Function& func) {
    auto result = func.call<FunctionType::Vec3D>(FunctionType::Vec3D {1, 0, 0});
    return std::abs(result[0] - 2) < 1e-5;
  });

  EXPECT_TRUE(inlet.verify());
}

TEST(inlet_function, simple_vec3_to_vec3_verify_lambda_fail)
{
  std::string testString = "function foo (x, y, z) return 2*x, 2*y, 2*z end";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  auto& func = inlet
                 .addFunction("foo",
                              FunctionTag::Vec3D,
                              {FunctionTag::Vec3D},
                              "foo's description")
                 .required();
  func.registerVerifier([](const axom::inlet::Function& func) {
    auto result = func.call<FunctionType::Vec3D>(FunctionType::Vec3D {2, 0, 0});
    return std::abs(result[0] - 2) < 1e-5;
  });

  EXPECT_FALSE(inlet.verify());
}

struct Foo
{
  bool bar;
  std::function<FunctionType::Vec3D(FunctionType::Vec3D)> baz;
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

TEST(inlet_function, simple_vec3_to_vec3_struct)
{
  std::string testString =
    "foo = { bar = true; baz = function (x, y, z) return 2*x, 2*y, 2*z end }";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  // Define schema
  inlet.addBool("foo/bar", "bar's description");
  inlet
    .addFunction("foo/baz",
                 FunctionTag::Vec3D,
                 {FunctionTag::Vec3D},
                 "baz's description")
    .required();
  Foo foo = inlet["foo"].get<Foo>();
  EXPECT_TRUE(foo.bar);
  auto result = foo.baz({4, 5, 6});
  EXPECT_FLOAT_EQ(result[0], 8);
  EXPECT_FLOAT_EQ(result[1], 10);
  EXPECT_FLOAT_EQ(result[2], 12);
}

TEST(inlet_function, simple_vec3_to_vec3_array_of_struct)
{
  std::string testString =
    "foo = { [7] = { bar = true, "
    "                baz = function (x, y, z) return 2*x, 2*y, 2*z end }, "
    "       [12] = { bar = false, "
    "                baz = function (x, y, z) return 3*x, 3*y, 3*z end } "
    "}";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  auto& arr_table = inlet.addGenericArray("foo");

  // Define schema
  arr_table.addBool("bar", "bar's description");
  arr_table
    .addFunction("baz",
                 FunctionTag::Vec3D,
                 {FunctionTag::Vec3D},
                 "baz's description")
    .required();

  auto foos = inlet["foo"].get<std::unordered_map<int, Foo>>();
  EXPECT_TRUE(foos[7].bar);
  auto first_result = foos[7].baz({4, 5, 6});
  EXPECT_FLOAT_EQ(first_result[0], 8);
  EXPECT_FLOAT_EQ(first_result[1], 10);
  EXPECT_FLOAT_EQ(first_result[2], 12);

  EXPECT_FALSE(foos[12].bar);
  auto second_result = foos[12].baz({4, 5, 6});
  EXPECT_FLOAT_EQ(second_result[0], 12);
  EXPECT_FLOAT_EQ(second_result[1], 15);
  EXPECT_FLOAT_EQ(second_result[2], 18);
}

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

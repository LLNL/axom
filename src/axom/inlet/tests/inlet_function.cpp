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
  std::string testString = "function foo (v) return v.x + v.y + v.z end";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  auto func =
    inlet.reader().getFunction("foo", FunctionTag::Double, {FunctionTag::Vector});

  EXPECT_TRUE(func);
  auto result = func.call<double>(FunctionType::Vector {1, 2, 3});
  EXPECT_FLOAT_EQ(result, 6);
}

TEST(inlet_function, simple_vec3_to_vec3_raw)
{
  std::string testString = "function foo (v) return 2*v end";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  auto func =
    inlet.reader().getFunction("foo", FunctionTag::Vector, {FunctionTag::Vector});

  EXPECT_TRUE(func);
  auto result = func.call<FunctionType::Vector>(FunctionType::Vector {1, 2, 3});
  EXPECT_FLOAT_EQ(result[0], 2);
  EXPECT_FLOAT_EQ(result[1], 4);
  EXPECT_FLOAT_EQ(result[2], 6);
}

TEST(inlet_function, simple_vec3_to_vec3_raw_partial_init)
{
  std::string testString = "function foo (v) return 2*v end";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  auto func =
    inlet.reader().getFunction("foo", FunctionTag::Vector, {FunctionTag::Vector});

  EXPECT_TRUE(func);

  auto result = func.call<FunctionType::Vector>(FunctionType::Vector {1, 2});
  EXPECT_FLOAT_EQ(result[0], 2);
  EXPECT_FLOAT_EQ(result[1], 4);
  EXPECT_FLOAT_EQ(result[2], 0);

  result = func.call<FunctionType::Vector>(FunctionType::Vector {1});
  EXPECT_FLOAT_EQ(result[0], 2);
  EXPECT_FLOAT_EQ(result[1], 0);
  EXPECT_FLOAT_EQ(result[2], 0);
}

TEST(inlet_function, simple_vec3_to_double_through_table)
{
  std::string testString = "function foo (v) return v.x + v.y + v.z end";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  inlet.addFunction("foo",
                    FunctionTag::Double,
                    {FunctionTag::Vector},
                    "foo's description");

  auto callable = inlet["foo"].get<std::function<double(FunctionType::Vector)>>();
  auto result = callable({1, 2, 3});
  EXPECT_FLOAT_EQ(result, 6);
}

TEST(inlet_function, simple_vec3_to_vec3_through_table)
{
  std::string testString = "function foo (v) return 2*v end";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  inlet.addFunction("foo",
                    FunctionTag::Vector,
                    {FunctionTag::Vector},
                    "foo's description");

  auto callable =
    inlet["foo"].get<std::function<FunctionType::Vector(FunctionType::Vector)>>();
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
  std::string testString = "function foo (v) return v.x + v.y + v.z end";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  inlet.addFunction("foo",
                    FunctionTag::Double,
                    {FunctionTag::Vector},
                    "foo's description");

  auto result = inlet["foo"].call<double>(FunctionType::Vector {1, 2, 3});
  EXPECT_FLOAT_EQ(result, 6);
}

TEST(inlet_function, simple_vec3_to_vec3_through_table_call)
{
  std::string testString = "function foo (v) return 2*v end";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  inlet.addFunction("foo",
                    FunctionTag::Vector,
                    {FunctionTag::Vector},
                    "foo's description");

  auto result =
    inlet["foo"].call<FunctionType::Vector>(FunctionType::Vector {1, 2, 3});
  EXPECT_FLOAT_EQ(result[0], 2);
  EXPECT_FLOAT_EQ(result[1], 4);
  EXPECT_FLOAT_EQ(result[2], 6);
}

TEST(inlet_function, simple_vec3_double_to_double_through_table_call)
{
  std::string testString =
    "function foo (v, t) return t * (v.x + v.y + v.z) end";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  inlet.addFunction("foo",
                    FunctionTag::Double,
                    {FunctionTag::Vector, FunctionTag::Double},
                    "foo's description");

  auto result = inlet["foo"].call<double>(FunctionType::Vector {1, 2, 3}, 2.0);
  EXPECT_FLOAT_EQ(result, 12);
}

TEST(inlet_function, simple_vec3_double_to_vec3_through_table_call)
{
  std::string testString = "function foo (v, t) return t*v end";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  inlet.addFunction("foo",
                    FunctionTag::Vector,
                    {FunctionTag::Vector, FunctionTag::Double},
                    "foo's description");

  auto result =
    inlet["foo"].call<FunctionType::Vector>(FunctionType::Vector {1, 2, 3}, 2.0);
  EXPECT_FLOAT_EQ(result[0], 2);
  EXPECT_FLOAT_EQ(result[1], 4);
  EXPECT_FLOAT_EQ(result[2], 6);
}

TEST(inlet_function, simple_vec3_to_vec3_verify_lambda_pass)
{
  std::string testString = "function foo (v) return 2*v end";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  auto& func = inlet
                 .addFunction("foo",
                              FunctionTag::Vector,
                              {FunctionTag::Vector},
                              "foo's description")
                 .required();
  func.registerVerifier([](const axom::inlet::Function& func) {
    auto result = func.call<FunctionType::Vector>(FunctionType::Vector {1, 0, 0});
    return std::abs(result[0] - 2) < 1e-5;
  });

  EXPECT_TRUE(inlet.verify());
}

TEST(inlet_function, simple_vec3_to_vec3_verify_lambda_fail)
{
  std::string testString = "function foo (v) return 2*v end";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  auto& func = inlet
                 .addFunction("foo",
                              FunctionTag::Vector,
                              {FunctionTag::Vector},
                              "foo's description")
                 .required();
  func.registerVerifier([](const axom::inlet::Function& func) {
    auto result = func.call<FunctionType::Vector>(FunctionType::Vector {2, 0, 0});
    return std::abs(result[0] - 2) < 1e-5;
  });

  EXPECT_FALSE(inlet.verify());
}

struct Foo
{
  bool bar;
  std::function<FunctionType::Vector(FunctionType::Vector)> baz;
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
    "foo = { bar = true; baz = function (v) return 2*v end }";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  // Define schema
  inlet.addBool("foo/bar", "bar's description");
  inlet
    .addFunction("foo/baz",
                 FunctionTag::Vector,
                 {FunctionTag::Vector},
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
    "                baz = function (v) return 2*v end }, "
    "       [12] = { bar = false, "
    "                baz = function (v) return 3*v end } "
    "}";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  auto& arr_table = inlet.addGenericArray("foo");

  // Define schema
  arr_table.addBool("bar", "bar's description");
  arr_table
    .addFunction("baz",
                 FunctionTag::Vector,
                 {FunctionTag::Vector},
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

TEST(inlet_function, dimension_dependent_result)
{
  std::string testString =
    "foo = function (v) "
    "first = 2 * v.x "
    "last = 2 * v.y "
    "if v:dim() == 2 then "
    "return Vector.new(first, last) "
    "else "
    "return Vector.new(first, 0, last) "
    "end "
    "end";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  inlet.addFunction("foo",
                    FunctionTag::Vector,
                    {FunctionTag::Vector},
                    "foo's description");

  auto callable =
    inlet["foo"].get<std::function<FunctionType::Vector(FunctionType::Vector)>>();

  FunctionType::Vector input_3d({3.5, 0.5, 7.5});
  auto result = callable(input_3d);

  EXPECT_EQ(result.dim, 3);
  EXPECT_FLOAT_EQ(result.vec[0], 7);
  EXPECT_FLOAT_EQ(result.vec[1], 0);
  EXPECT_FLOAT_EQ(result.vec[2], 1);

  FunctionType::Vector input_2d({3.5, 0.5});
  result = callable(input_2d);

  EXPECT_EQ(result.dim, 2);
  EXPECT_FLOAT_EQ(result.vec[0], 7);
  EXPECT_FLOAT_EQ(result.vec[1], 1);
}

template <typename Ret, typename... Args>
Ret checkedCall(const sol::protected_function& func, Args&&... args)
{
  auto tentative_result = func(std::forward<Args>(args)...);
  EXPECT_TRUE(tentative_result.valid());
  return tentative_result;
}

TEST(inlet_function, lua_usertype_basic)
{
  std::string testString = "function func(vec) return 7 end";
  LuaReader lr;
  lr.parseString(testString);
  sol::protected_function func = lr.solState()["func"];
  axom::inlet::FunctionType::Vector vec {1, 2, 3};
  int result = checkedCall<int>(func, vec);
  EXPECT_EQ(result, 7);
}

TEST(inlet_function, lua_usertype_basic_ret)
{
  std::string testString =
    "function func(x, y, z) return Vector.new(x, y, z) end";
  LuaReader lr;
  lr.parseString(testString);
  sol::protected_function func = lr.solState()["func"];
  axom::inlet::FunctionType::Vector vec {1, 2, 3};
  auto result = checkedCall<axom::inlet::FunctionType::Vector>(func, 1, 2, 3);
  EXPECT_EQ(vec, result);
}

TEST(inlet_function, lua_usertype_basic_ret_2d)
{
  std::string testString = "function func(x, y, z) return Vector.new(x, y) end";
  LuaReader lr;
  lr.parseString(testString);
  sol::protected_function func = lr.solState()["func"];
  axom::inlet::FunctionType::Vector vec {1, 2};
  auto result = checkedCall<axom::inlet::FunctionType::Vector>(func, 1, 2, 3);
  EXPECT_EQ(vec, result);
}

TEST(inlet_function, lua_usertype_basic_ret_default)
{
  std::string testString = "function func(x, y, z) return Vector.new() end";
  LuaReader lr;
  lr.parseString(testString);
  sol::protected_function func = lr.solState()["func"];
  axom::inlet::FunctionType::Vector vec {0, 0, 0};
  auto result = checkedCall<axom::inlet::FunctionType::Vector>(func, 1, 2, 3);
  EXPECT_EQ(vec, result);
}

TEST(inlet_function, lua_usertype_basic_add)
{
  std::string testString = "function func(vec1, vec2) return vec1 + vec2 end";
  LuaReader lr;
  lr.parseString(testString);
  sol::protected_function func = lr.solState()["func"];
  axom::inlet::FunctionType::Vector vec1 {1, 2, 3};
  axom::inlet::FunctionType::Vector vec2 {4, 5, 6};
  auto result = checkedCall<axom::inlet::FunctionType::Vector>(func, vec1, vec2);
  EXPECT_EQ(result, vec1.vec + vec2.vec);
}

TEST(inlet_function, lua_usertype_basic_negate)
{
  std::string testString = "function func(vec) return -vec end";
  LuaReader lr;
  lr.parseString(testString);
  sol::protected_function func = lr.solState()["func"];
  axom::inlet::FunctionType::Vector vec {1, 2, 3};
  auto result = checkedCall<axom::inlet::FunctionType::Vector>(func, vec);
  EXPECT_EQ(result, -vec.vec);
}

TEST(inlet_function, lua_usertype_basic_scalar_mult)
{
  std::string testString =
    "function func1(vec, x) return vec * x end; function func2(vec, x) return "
    "x * vec end";
  LuaReader lr;
  lr.parseString(testString);
  sol::protected_function func1 = lr.solState()["func1"];
  sol::protected_function func2 = lr.solState()["func2"];
  axom::inlet::FunctionType::Vector vec {1, 2, 3};
  auto result = checkedCall<axom::inlet::FunctionType::Vector>(func1, vec, 2.0);
  EXPECT_EQ(result, 2.0 * vec.vec);
  result = checkedCall<axom::inlet::FunctionType::Vector>(func2, vec, 3.0);
  EXPECT_EQ(result, 3.0 * vec.vec);
}

TEST(inlet_function, lua_usertype_basic_index_get)
{
  std::string testString = "function func(vec, idx) return vec[idx] end";
  LuaReader lr;
  lr.parseString(testString);
  sol::protected_function func = lr.solState()["func"];
  axom::inlet::FunctionType::Vector vec {1, 2, 3};
  // Use 1-based indexing in these tests as lua is 1-indexed
  auto result = checkedCall<double>(func, vec, 1);
  EXPECT_FLOAT_EQ(1, result);
  result = checkedCall<double>(func, vec, 2);
  EXPECT_FLOAT_EQ(2, result);
  result = checkedCall<double>(func, vec, 3);
  EXPECT_FLOAT_EQ(3, result);
}

TEST(inlet_function, lua_usertype_basic_index_set)
{
  std::string testString =
    "function func(idx) vec = Vector.new(1,1,1); vec[idx] = -1; return vec end";
  LuaReader lr;
  lr.parseString(testString);
  sol::protected_function func = lr.solState()["func"];
  auto result = checkedCall<axom::inlet::FunctionType::Vector>(func, 1);
  EXPECT_FLOAT_EQ(-1, result[0]);
  result = checkedCall<axom::inlet::FunctionType::Vector>(func, 2);
  EXPECT_FLOAT_EQ(-1, result[1]);
  result = checkedCall<axom::inlet::FunctionType::Vector>(func, 3);
  EXPECT_FLOAT_EQ(-1, result[2]);
}

TEST(inlet_function, lua_usertype_basic_norm)
{
  std::string testString = "function func(vec) return vec:norm() end";
  LuaReader lr;
  lr.parseString(testString);
  sol::protected_function func = lr.solState()["func"];
  axom::inlet::FunctionType::Vector vec {1, 2, 3};
  auto result = checkedCall<double>(func, vec);
  EXPECT_FLOAT_EQ(vec.vec.norm(), result);
}

TEST(inlet_function, lua_usertype_basic_squared_norm)
{
  std::string testString = "function func(vec) return vec:squared_norm() end";
  LuaReader lr;
  lr.parseString(testString);
  sol::protected_function func = lr.solState()["func"];
  axom::inlet::FunctionType::Vector vec {1, 2, 3};
  auto result = checkedCall<double>(func, vec);
  EXPECT_FLOAT_EQ(vec.vec.squared_norm(), result);
}

TEST(inlet_function, lua_usertype_basic_unit_vec)
{
  std::string testString = "function func(vec) return vec:unitVector() end";
  LuaReader lr;
  lr.parseString(testString);
  sol::protected_function func = lr.solState()["func"];
  axom::inlet::FunctionType::Vector vec {1, 2, 3};
  auto result = checkedCall<axom::inlet::FunctionType::Vector>(func, vec);
  EXPECT_EQ(vec.vec.unitVector(), result);
}

TEST(inlet_function, lua_usertype_basic_dot)
{
  std::string testString =
    "function func(vec1, vec2) return vec1:dot(vec2) end";
  LuaReader lr;
  lr.parseString(testString);
  sol::protected_function func = lr.solState()["func"];
  axom::inlet::FunctionType::Vector vec1 {1, 2, 3};
  axom::inlet::FunctionType::Vector vec2 {4, 5, 6};
  auto result = checkedCall<double>(func, vec1, vec2);
  EXPECT_EQ(vec1.vec.dot(vec2), result);
}

TEST(inlet_function, lua_usertype_basic_cross)
{
  std::string testString =
    "function func(vec1, vec2) return vec1:cross(vec2) end";
  LuaReader lr;
  lr.parseString(testString);
  sol::protected_function func = lr.solState()["func"];
  axom::inlet::FunctionType::Vector vec1 {1, 2, 3};
  axom::inlet::FunctionType::Vector vec2 {4, 5, 6};
  auto result = checkedCall<axom::inlet::FunctionType::Vector>(func, vec1, vec2);
  EXPECT_EQ(axom::primal::Vector3D::cross_product(vec1.vec, vec2.vec), result);
}

TEST(inlet_function, lua_usertype_check_dim)
{
  std::string testString = "function func(vec) return vec:dim() end";
  LuaReader lr;
  lr.parseString(testString);
  sol::protected_function func = lr.solState()["func"];
  axom::inlet::FunctionType::Vector vec1 {1, 2, 3};
  axom::inlet::FunctionType::Vector vec2 {4, 5};
  auto result = checkedCall<double>(func, vec1);
  EXPECT_EQ(result, 3);

  result = checkedCall<double>(func, vec2);
  EXPECT_EQ(result, 2);
}

TEST(inlet_function, lua_usertype_named_access)
{
  std::string testString =
    "function func(vec, comp) if comp == 1 then return vec.x elseif comp == 2 "
    "then return vec.y else return vec.z end end";
  LuaReader lr;
  lr.parseString(testString);
  sol::protected_function func = lr.solState()["func"];
  axom::inlet::FunctionType::Vector vec1 {4, 5, 6};
  auto result = checkedCall<double>(func, vec1, 1);
  EXPECT_EQ(result, 4);

  result = checkedCall<double>(func, vec1, 2);
  EXPECT_EQ(result, 5);

  result = checkedCall<double>(func, vec1, 3);
  EXPECT_EQ(result, 6);
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

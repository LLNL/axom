// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/slic.hpp"

#include "axom/core/utilities/FileUtilities.hpp"

#include "axom/inlet/LuaReader.hpp"
#include "axom/inlet/YAMLReader.hpp"
#include "axom/inlet/Inlet.hpp"
#include "axom/inlet/JSONSchemaWriter.hpp"
#include "axom/inlet/SphinxWriter.hpp"

#include "gtest/gtest.h"

#include <array>
#include <functional>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
#include <iostream>

using axom::inlet::FunctionTag;
using axom::inlet::FunctionType;
using axom::inlet::Inlet;
using axom::inlet::InletType;
using axom::inlet::LuaReader;
using axom::inlet::VerificationError;

#include "axom/sol.hpp"

class SolStateReader : public LuaReader
{
public:
  using LuaReader::solState;
};

Inlet createBasicInlet(const std::string& luaString, bool enableDocs = true)
{
  auto lr = std::make_unique<LuaReader>();
  lr->parseString(luaString);
  return Inlet(std::move(lr), enableDocs);
}

/// Collects SLIC error messages; registered with the logger in main()
std::ostringstream& slicErrors()
{
  static std::ostringstream errors;
  return errors;
}

/// Checks whether a SLIC error containing \a expected was logged, then clears the log
::testing::AssertionResult slicErrorContains(const std::string& expected)
{
  const std::string errors = slicErrors().str();
  slicErrors().str("");
  if(errors.find(expected) != std::string::npos)
  {
    return ::testing::AssertionSuccess();
  }
  return ::testing::AssertionFailure()
    << "Expected a SLIC error containing '" << expected << "', got: '" << errors << "'";
}

TEST(inlet_function, simple_vec3_to_double_raw)
{
  std::string testString = "function foo (v) return v.x + v.y + v.z end";
  auto inlet = createBasicInlet(testString);

  auto func = inlet.reader().getFunction("foo", FunctionTag::Double, {FunctionTag::Vector});

  EXPECT_TRUE(func);
  auto result = func.call<double>(FunctionType::Vector {1, 2, 3});
  EXPECT_FLOAT_EQ(result, 6);
}

TEST(inlet_function, simple_vec3_to_vec3_raw)
{
  std::string testString = "function foo (v) return 2*v end";
  auto inlet = createBasicInlet(testString);

  auto func = inlet.reader().getFunction("foo", FunctionTag::Vector, {FunctionTag::Vector});

  EXPECT_TRUE(func);
  auto result = func.call<FunctionType::Vector>(FunctionType::Vector {1, 2, 3});
  EXPECT_FLOAT_EQ(result[0], 2);
  EXPECT_FLOAT_EQ(result[1], 4);
  EXPECT_FLOAT_EQ(result[2], 6);
}

TEST(inlet_function, vector_function_accepts_lua_table_returns)
{
  auto inlet = createBasicInlet(R"(
    function make_vector (dim)
      local result = {}
      for i = 1, dim do result[i] = i end
      return result
    end
  )");

  auto function =
    inlet.reader().getFunction("make_vector", FunctionTag::Vector, {FunctionTag::Double});
  ASSERT_TRUE(function);
  for(int dim = 1; dim <= 3; ++dim)
  {
    const auto result = function.call<FunctionType::Vector>(static_cast<double>(dim));
    ASSERT_EQ(result.dim, dim);
    for(int component = 0; component < result.dim; ++component)
    {
      EXPECT_DOUBLE_EQ(result[component], component + 1.0);
    }
  }
}

TEST(inlet_function, vector_function_accepts_equivalent_table_forms)
{
  // Explicit keys, in any order, and float keys with integral values are the same vector
  for(const std::string table : {"{4.5, -2}", "{[2] = -2, [1] = 4.5}", "{[1.0] = 4.5, [2.0] = -2}"})
  {
    auto inlet = createBasicInlet("function foo () return " + table + " end");
    inlet.addFunction("foo", FunctionTag::Vector, {});
    ASSERT_TRUE(inlet.verify()) << table;

    // Through Proxy::call and through a std::function copy
    const auto viaCall = inlet["foo"].call<FunctionType::Vector>();
    const auto viaCopy = inlet["foo"].get<std::function<FunctionType::Vector()>>()();
    for(const auto& result : {viaCall, viaCopy})
    {
      ASSERT_EQ(result.dim, 2) << table;
      EXPECT_DOUBLE_EQ(result[0], 4.5) << table;
      EXPECT_DOUBLE_EQ(result[1], -2.0) << table;
    }
  }
}

TEST(inlet_function, vector_function_rejects_malformed_table_returns)
{
  axom::slic::ScopedAbortToThrow abortGuard;

  // Lua vectors must be dense numeric sequences with a supported dimension.
  // Each case fails exactly one check; match its message to confirm which one.
  const std::array<std::pair<std::string, std::string>, 13> cases {{
    {"2.0", "return types possibly incorrect"},
    {"'text'", "return types possibly incorrect"},
    {"{}", "empty table"},
    {"{1, 2, 3, 4}", "integers between 1 and 3"},
    {"{[1] = 1, [3] = 3}", "contiguous starting at 1"},
    {"{1, 'two'}", "components must be numeric"},
    {"{1, true}", "components must be numeric"},
    {"{1, {2}}", "components must be numeric"},
    {"{1, 2, label = 3}", "only contain numeric indices"},
    {"{[0] = 0, 1, 2}", "integers between 1 and 3"},
    {"{[1.5] = 1}", "integers between 1 and 3"},
    {"{[-1] = 1}", "integers between 1 and 3"},
    {"{[math.huge] = 1}", "integers between 1 and 3"},
  }};

  for(const auto& [result, message] : cases)
  {
    auto inlet = createBasicInlet("function foo () return " + result + " end");
    auto func = inlet.reader().getFunction("foo", FunctionTag::Vector, {});
    ASSERT_TRUE(func);
    EXPECT_THROW(func.call<FunctionType::Vector>(), axom::slic::SlicAbortException) << result;
    EXPECT_TRUE(slicErrorContains(message)) << result;
  }
}

TEST(inlet_function, lua_callback_failures_report_slic_errors)
{
  axom::slic::ScopedAbortToThrow abortGuard;

  auto inlet = createBasicInlet(R"(
    function runtime_error () error('callback failed') end
    function wrong_type () return 'not a number' end
  )");

  auto wrongType = inlet.reader().getFunction("wrong_type", FunctionTag::Double, {});
  ASSERT_TRUE(wrongType);
  EXPECT_THROW(wrongType.call<FunctionType::Double>(), axom::slic::SlicAbortException);
  EXPECT_TRUE(slicErrorContains("return types possibly incorrect"));

  // The SLIC error includes Lua's error message
  auto runtimeError = inlet.reader().getFunction("runtime_error", FunctionTag::Double, {});
  ASSERT_TRUE(runtimeError);
  EXPECT_THROW(runtimeError.call<FunctionType::Double>(), axom::slic::SlicAbortException);
  EXPECT_TRUE(slicErrorContains("callback failed"));
}

TEST(inlet_function, lua_callback_failures_report_slic_errors_through_every_access_path)
{
  axom::slic::ScopedAbortToThrow abortGuard;

  auto inlet = createBasicInlet(R"(
    function fail_void (x) error('void failed') end
    function add_one (x) return x + 1 end
  )");
  inlet.addFunction("fail_void", FunctionTag::Void, {FunctionTag::Double});
  inlet.addFunction("add_one", FunctionTag::Double, {FunctionTag::String});
  ASSERT_TRUE(inlet.verify());

  // A void callback still reports an execution error
  EXPECT_THROW(inlet["fail_void"].call<void>(1.0), axom::slic::SlicAbortException);
  EXPECT_TRUE(slicErrorContains("void failed"));

  // An invalid input reports an execution error through Proxy::call.
  EXPECT_THROW(inlet["add_one"].call<double>(std::string {"text"}), axom::slic::SlicAbortException);
  EXPECT_TRUE(slicErrorContains("Lua function call failed:"));

  // Copies made with get<std::function> report errors the same way
  auto addOne = inlet["add_one"].get<std::function<double(std::string)>>();
  EXPECT_DOUBLE_EQ(addOne("2"), 3.0);
  EXPECT_THROW(addOne("text"), axom::slic::SlicAbortException);
  EXPECT_TRUE(slicErrorContains("Lua function call failed:"));
}

TEST(inlet_function, lua_callback_failure_in_verifier_reports_slic_error)
{
  axom::slic::ScopedAbortToThrow abortGuard;

  auto inlet = createBasicInlet("function scale () error('scale failed') end");
  inlet.addFunction("scale", FunctionTag::Double, {})
    .registerVerifier([](const axom::inlet::Function& func) { return func.call<double>() > 0.0; });

  EXPECT_THROW(inlet.verify(), axom::slic::SlicAbortException);
  EXPECT_TRUE(slicErrorContains("scale failed"));
}

TEST(inlet_function, simple_vec3_to_vec3_raw_partial_init)
{
  std::string testString = "function foo (v) return 2*v end";
  auto inlet = createBasicInlet(testString);

  auto func = inlet.reader().getFunction("foo", FunctionTag::Vector, {FunctionTag::Vector});

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

TEST(inlet_function, simple_vec3_to_double_through_container)
{
  std::string testString = "function foo (v) return v.x + v.y + v.z end";
  auto inlet = createBasicInlet(testString);

  inlet.addFunction("foo", FunctionTag::Double, {FunctionTag::Vector}, "foo's description");

  auto callable = inlet["foo"].get<std::function<double(FunctionType::Vector)>>();
  auto result = callable({1, 2, 3});
  EXPECT_FLOAT_EQ(result, 6);
}

TEST(inlet_function, simple_vec3_to_vec3_through_container)
{
  std::string testString = "function foo (v) return 2*v end";
  auto inlet = createBasicInlet(testString);

  inlet.addFunction("foo", FunctionTag::Vector, {FunctionTag::Vector}, "foo's description");

  auto callable = inlet["foo"].get<std::function<FunctionType::Vector(FunctionType::Vector)>>();
  auto result = callable({1, 2, 3});
  EXPECT_FLOAT_EQ(result[0], 2);
  EXPECT_FLOAT_EQ(result[1], 4);
  EXPECT_FLOAT_EQ(result[2], 6);
}

TEST(inlet_function, simple_double_to_double_through_container)
{
  std::string testString = "function foo (a) return (a * 3.4) + 9.64 end";
  auto inlet = createBasicInlet(testString);

  inlet.addFunction("foo", FunctionTag::Double, {FunctionTag::Double}, "foo's description");

  auto callable = inlet["foo"].get<std::function<FunctionType::Double(FunctionType::Double)>>();
  double arg = -6.37;
  double result = callable(arg);
  EXPECT_FLOAT_EQ(result, (arg * 3.4) + 9.64);
}

TEST(inlet_function, returned_function_keeps_lua_state_alive)
{
  // An extracted callback must retain its Lua state after Inlet is destroyed.
  std::function<double(double)> callback;
  {
    auto inlet = createBasicInlet("offset = 3.0; function foo (value) return value + offset end");
    inlet.addFunction("foo", FunctionTag::Double, {FunctionTag::Double});
    callback = inlet["foo"].get<std::function<double(double)>>();
  }

  EXPECT_DOUBLE_EQ(callback(4.0), 7.0);
}

TEST(inlet_function, returned_functions_share_lua_state_after_inlet_is_destroyed)
{
  // Callbacks from one input share its Lua state, including after the Inlet is gone.
  // Copying, calling, and destroying them in any order must leave the state usable.
  std::function<double()> increment;
  std::function<double()> current;
  {
    auto inlet = createBasicInlet(R"(
      count = 0
      function increment () count = count + 1; return count end
      function current () return count end
    )");
    inlet.addFunction("increment", FunctionTag::Double, {});
    inlet.addFunction("current", FunctionTag::Double, {});
    increment = inlet["increment"].get<std::function<double()>>();
    current = inlet["current"].get<std::function<double()>>();
  }

  EXPECT_DOUBLE_EQ(increment(), 1.0);
  auto incrementCopy = increment;
  increment = nullptr;
  EXPECT_DOUBLE_EQ(incrementCopy(), 2.0);
  incrementCopy = nullptr;
  EXPECT_DOUBLE_EQ(current(), 2.0);
}

TEST(inlet_function, reader_function_outlives_reader)
{
  axom::inlet::FunctionVariant function;
  {
    LuaReader reader;
    reader.parseString("function foo (value) return 2 * value end");
    function = reader.getFunction("foo", FunctionTag::Double, {FunctionTag::Double});
  }

  ASSERT_TRUE(function);
  EXPECT_DOUBLE_EQ(function.call<double>(4.0), 8.0);
}

TEST(inlet_function, function_declaration_is_idempotent_across_containers)
{
  // Redeclaring a function returns the existing one, whichever Container declared it first
  for(const bool firstOnRoot : {true, false})
  {
    for(const bool secondOnRoot : {true, false})
    {
      auto inlet = createBasicInlet("group = { f = function() return 2.0 end }");
      auto& group = inlet.addStruct("group");

      auto declare = [&](bool onRoot) -> axom::inlet::Verifiable<axom::inlet::Function>& {
        return onRoot ? inlet.addFunction("group/f", FunctionTag::Double, {})
                      : group.addFunction("f", FunctionTag::Double, {});
      };

      auto& first = declare(firstOnRoot);
      auto& second = declare(secondOnRoot);
      EXPECT_EQ(&first, &second) << firstOnRoot << secondOnRoot;
      EXPECT_EQ(1u, group.getChildFunctions().size());
      EXPECT_TRUE(inlet.getGlobalContainer().getChildFunctions().empty());
      EXPECT_TRUE(inlet.verify());
      EXPECT_DOUBLE_EQ(group["f"].call<double>(), 2.0);
    }
  }
}

TEST(inlet_function, function_value_alternative_selects_supplied_representation)
{
  auto inlet = createBasicInlet(R"(
    function label () return 'computed' end
    scale = 4.0
  )");

  inlet.addFunctionAsValueAlternative("label", FunctionTag::String, {});
  inlet.addString("label");

  inlet.addFunctionAsValueAlternative("scale", FunctionTag::Double, {});
  inlet.addDouble("scale");

  EXPECT_TRUE(inlet.verify());
  EXPECT_FALSE(inlet.contains("label"));
  ASSERT_TRUE(inlet.containsFunctionValueAlternative("label"));
  EXPECT_EQ(inlet.getFunctionValueAlternative("label").call<std::string>(), "computed");
  EXPECT_EQ(inlet.getGlobalContainer().getFunctionValueAlternativeNames(),
            std::vector<std::string> {"label"});

  EXPECT_TRUE(inlet.contains("scale"));
  EXPECT_FALSE(inlet.containsFunctionValueAlternative("scale"));
  EXPECT_DOUBLE_EQ(inlet["scale"].get<double>(), 4.0);
}

TEST(inlet_function, function_value_alternative_rejects_unrelated_wrong_type)
{
  auto inlet = createBasicInlet("foo = 'not a number or function'");
  inlet.addFunctionAsValueAlternative("foo", FunctionTag::Double, {});
  inlet.addDouble("foo");

  EXPECT_FALSE(inlet.verify());
  EXPECT_FALSE(inlet.contains("foo"));
  EXPECT_FALSE(inlet.getGlobalContainer().containsFunctionValueAlternative("foo"));
  // The input exists even though neither schema entry accepts its type.
  EXPECT_TRUE(inlet.isUserProvided("foo"));
  EXPECT_FALSE(inlet.getGlobalContainer().exists());
}

TEST(inlet_function, function_value_alternative_is_container_independent)
{
  for(const bool functionOnRoot : {true, false})
  {
    auto inlet = createBasicInlet("group = { value = function() return 2.0 end }");
    inlet.getGlobalContainer().strict();
    auto& group = inlet.addStruct("group");

    // The alternative and the concrete entry may be declared through different Containers,
    // so long as the alternative comes first.
    if(functionOnRoot)
    {
      inlet.addFunctionAsValueAlternative("group/value", FunctionTag::Double, {});
      group.addDouble("value");
    }
    else
    {
      group.addFunctionAsValueAlternative("value", FunctionTag::Double, {});
      inlet.addDouble("group/value");
    }

    EXPECT_TRUE(inlet.verify());
    EXPECT_TRUE(inlet.unexpectedNames().empty());
    EXPECT_TRUE(inlet.getGlobalContainer().containsFunctionValueAlternative("group/value"));
    EXPECT_TRUE(group.containsFunctionValueAlternative("value"));
    EXPECT_FALSE(inlet.contains("group/value"));
    EXPECT_TRUE(inlet.isUserProvided("group/value"));
    EXPECT_TRUE(group.isUserProvided("value"));
    EXPECT_TRUE(group.exists());
    EXPECT_DOUBLE_EQ(group.getFunctionValueAlternative("value").call<double>(), 2.0);
    EXPECT_EQ(group.getFunctionValueAlternativeNames(), std::vector<std::string> {"value"});
    EXPECT_TRUE(inlet.getGlobalContainer().getFunctionValueAlternativeNames().empty());
  }
}

TEST(inlet_function, function_value_alternative_names_are_sorted)
{
  auto inlet = createBasicInlet(R"(
    zeta = function() return 1.0 end
    alpha = function() return 2.0 end
    plain = 3.0
  )");

  inlet.addFunctionAsValueAlternative("zeta", FunctionTag::Double, {});
  inlet.addFunctionAsValueAlternative("alpha", FunctionTag::Double, {});
  // Declared but not supplied by the input, so it is not reported
  inlet.addFunctionAsValueAlternative("plain", FunctionTag::Double, {});
  inlet.addDouble("plain");

  EXPECT_TRUE(inlet.verify());
  EXPECT_EQ(inlet.getGlobalContainer().getFunctionValueAlternativeNames(),
            (std::vector<std::string> {"alpha", "zeta"}));
}

TEST(inlet_function, ordinary_function_suffix_does_not_register_an_alternative)
{
  for(const bool enableDocs : {true, false})
  {
    auto inlet = createBasicInlet(R"(
      scale = function() return 2.0 end
      scale_inlet_function_alternative = function() return 9.0 end
      offset = function() return 3.0 end
    )",
                                  enableDocs);
    inlet.addFunction("scale_inlet_function_alternative", FunctionTag::Double, {});
    inlet.addFunctionAsValueAlternative("offset", FunctionTag::Double, {});
    inlet.addDouble("offset");

    auto& container = inlet.getGlobalContainer();
    EXPECT_FALSE(container.containsFunctionValueAlternative("scale"));
    EXPECT_EQ(container.getFunctionValueAlternativeNames(), std::vector<std::string> {"offset"});
    EXPECT_DOUBLE_EQ(inlet["scale_inlet_function_alternative"].call<double>(), 9.0);
    EXPECT_DOUBLE_EQ(container.getFunctionValueAlternative("offset").call<double>(), 3.0);
    {
      axom::slic::ScopedAbortToThrow abortGuard;
      EXPECT_THROW(container.getFunctionValueAlternative("scale"), axom::slic::SlicAbortException);
    }

    // An unrelated ordinary function must not suppress a concrete field's type error.
    inlet.addDouble("scale");
    EXPECT_FALSE(inlet.verify());
  }
}

TEST(inlet_function, function_value_alternative_rejects_function_name_collisions)
{
  axom::slic::ScopedAbortToThrow abortGuard;
  for(const bool alternativeFirst : {true, false})
  {
    for(const bool alternativeOnRoot : {true, false})
    {
      auto inlet = createBasicInlet(R"(
        group = {
          scale = function() return 2.0 end,
          scale_inlet_function_alternative = function() return 9.0 end
        }
      )");
      auto& group = inlet.addStruct("group");
      auto& root = inlet.getGlobalContainer();
      auto& alternativeContainer = alternativeOnRoot ? root : group;
      auto& ordinaryContainer = alternativeOnRoot ? group : root;
      const std::string valueName = alternativeOnRoot ? "group/scale" : "scale";
      const std::string functionName = alternativeOnRoot ? "scale_inlet_function_alternative"
                                                         : "group/scale_inlet_function_alternative";
      if(alternativeFirst)
      {
        alternativeContainer.addFunctionAsValueAlternative(valueName, FunctionTag::Double, {});
        EXPECT_THROW(ordinaryContainer.addFunction(functionName, FunctionTag::Double, {}),
                     axom::slic::SlicAbortException);
        EXPECT_DOUBLE_EQ(group.getFunctionValueAlternative("scale").call<double>(), 2.0);
      }
      else
      {
        ordinaryContainer.addFunction(functionName, FunctionTag::Double, {});
        EXPECT_THROW(
          alternativeContainer.addFunctionAsValueAlternative(valueName, FunctionTag::Double, {}),
          axom::slic::SlicAbortException);
        EXPECT_FALSE(group.containsFunctionValueAlternative("scale"));
        EXPECT_DOUBLE_EQ(group["scale_inlet_function_alternative"].call<double>(), 9.0);
      }
    }
  }
}

TEST(inlet_function, function_value_alternative_array_is_container_independent)
{
  for(const bool functionOnRoot : {true, false})
  {
    auto inlet = createBasicInlet("group = { values = function() return {1.0, 2.0, 3.0} end }");
    auto& group = inlet.addStruct("group");

    if(functionOnRoot)
    {
      inlet.addFunctionAsValueAlternative("group/values", FunctionTag::Vector, {});
    }
    else
    {
      group.addFunctionAsValueAlternative("values", FunctionTag::Vector, {});
    }
    inlet.addDoubleArray("group/values");

    EXPECT_TRUE(inlet.verify());
    EXPECT_FALSE(inlet.contains("group/values"));
    EXPECT_TRUE(group.containsFunctionValueAlternative("values"));
    const auto result = group.getFunctionValueAlternative("values").call<FunctionType::Vector>();
    EXPECT_DOUBLE_EQ(result[0], 1.0);
    EXPECT_DOUBLE_EQ(result[1], 2.0);
    EXPECT_DOUBLE_EQ(result[2], 3.0);
  }
}

TEST(inlet_function, function_value_alternative_applies_to_every_collection_kind)
{
  // Arrays and dictionaries use separate read paths.
  // When a function is supplied function, the former are treated as absent.
  const std::string input = R"(
    arr = function() return {1.0, 2.0} end
    dict = function() return {3.0} end
    varr = function() return {4.0} end
    vdict = function() return {5.0} end
    plain = {first = 1.0, second = 2.0}
  )";
  auto inlet = createBasicInlet(input);
  inlet.getGlobalContainer().strict();

  for(const std::string name : {"arr", "dict", "varr", "vdict", "plain"})
  {
    inlet.addFunctionAsValueAlternative(name, FunctionTag::Vector, {});
  }
  inlet.addDoubleArray("arr");
  inlet.addDoubleDictionary("dict");
  inlet.addVariantArray("varr");
  inlet.addVariantDictionary("vdict");
  inlet.addDoubleDictionary("plain");

  EXPECT_TRUE(inlet.verify());
  EXPECT_TRUE(inlet.unexpectedNames().empty());
  const auto& root = inlet.getGlobalContainer();
  EXPECT_EQ(root.getFunctionValueAlternativeNames(),
            (std::vector<std::string> {"arr", "dict", "varr", "vdict"}));
  for(const std::string name : {"arr", "dict", "varr", "vdict"})
  {
    EXPECT_FALSE(inlet.contains(name)) << name;
    EXPECT_TRUE(inlet.isUserProvided(name)) << name;
    EXPECT_TRUE(root.containsFunctionValueAlternative(name)) << name;
  }
  EXPECT_DOUBLE_EQ(root.getFunctionValueAlternative("vdict").call<FunctionType::Vector>()[0], 5.0);

  // A concrete collection is unaffected by its declared alternative
  EXPECT_TRUE(inlet.contains("plain"));
  EXPECT_FALSE(root.containsFunctionValueAlternative("plain"));
  const auto plain = inlet["plain"].get<std::unordered_map<std::string, double>>();
  EXPECT_DOUBLE_EQ(plain.at("first"), 1.0);
  EXPECT_DOUBLE_EQ(plain.at("second"), 2.0);
}

TEST(inlet_function, function_value_alternative_collection_rejects_unrelated_wrong_type)
{
  // A scalar is neither a collection nor a function
  auto inlet = createBasicInlet("values = 3.0");
  inlet.addFunctionAsValueAlternative("values", FunctionTag::Vector, {});
  inlet.addDoubleArray("values");

  EXPECT_FALSE(inlet.verify());
  EXPECT_FALSE(inlet.contains("values"));
  EXPECT_FALSE(inlet.getGlobalContainer().containsFunctionValueAlternative("values"));
}

TEST(inlet_function, function_value_alternative_declaration_is_idempotent)
{
  // Matches addFunction's behavior: redeclaring returns the existing entry
  auto inlet = createBasicInlet("group = { scale = function() return 2.0 end }");
  auto& group = inlet.addStruct("group");

  auto& first = inlet.addFunctionAsValueAlternative("group/scale", FunctionTag::Double, {});
  auto& second = group.addFunctionAsValueAlternative("scale", FunctionTag::Double, {});
  EXPECT_EQ(&first, &second);

  EXPECT_EQ(1u, group.getChildFunctions().size());
  EXPECT_TRUE(inlet.getGlobalContainer().getChildFunctions().empty());
  EXPECT_DOUBLE_EQ(group.getFunctionValueAlternative("scale").call<double>(), 2.0);
}

TEST(inlet_function, function_value_alternative_is_not_documented_under_its_schema_name)
{
  // The alternative is stored under an internal schema name
  // which should not leak into generated documentation
  axom::utilities::filesystem::TempFile doc_file("inlet_function_value_alternative_docs", ".rst");
  {
    auto inlet = createBasicInlet(R"(
      scale = function() return 2.0 end
      ordinary_inlet_function_alternative = function() return 9.0 end
    )");
    inlet.addFunctionAsValueAlternative("scale", FunctionTag::Double, {}, "a scale callback");
    inlet.addDouble("scale", "a scale");
    inlet.addFunction("ordinary_inlet_function_alternative", FunctionTag::Double, {});
    EXPECT_TRUE(inlet.verify());
    // The alternative is a real schema entry, just under an internal name
    EXPECT_EQ(2u, inlet.getGlobalContainer().getChildFunctions().size());
    inlet.write(axom::inlet::SphinxWriter(doc_file.getPath()));
  }

  const std::string contents = doc_file.getFileContents();
  EXPECT_EQ(std::string::npos, contents.find("scale_inlet_function_alternative"));
  EXPECT_NE(std::string::npos, contents.find("ordinary_inlet_function_alternative"));
  EXPECT_NE(std::string::npos, contents.find("scale"));
}

TEST(inlet_function, function_value_alternative_is_not_in_json_schema)
{
  // JSONSchemaWriter documents Fields only; check that the internal name stays out
  axom::utilities::filesystem::TempFile schema_file("inlet_function_value_alternative_schema",
                                                    ".json");
  {
    auto inlet = createBasicInlet("scale = function() return 2.0 end");
    inlet.addFunctionAsValueAlternative("scale", FunctionTag::Double, {}, "a scale callback");
    inlet.addDouble("scale", "a scale");
    EXPECT_TRUE(inlet.verify());
    inlet.write(axom::inlet::JSONSchemaWriter(schema_file.getPath()));
  }

  const std::string contents = schema_file.getFileContents();
  EXPECT_EQ(std::string::npos, contents.find("_inlet_function_alternative"));
  EXPECT_NE(std::string::npos, contents.find("scale"));
}

TEST(inlet_function, function_value_alternative_can_be_required)
{
  // The alternative is an ordinary Function, so schema constraints apply to it
  auto missing = createBasicInlet("scale = 2.0");
  missing.addFunctionAsValueAlternative("scale", FunctionTag::Double, {}).required();
  missing.addDouble("scale");
  EXPECT_FALSE(missing.verify());

  auto supplied = createBasicInlet("scale = function() return 2.0 end");
  supplied.addFunctionAsValueAlternative("scale", FunctionTag::Double, {}).required();
  supplied.addDouble("scale");
  EXPECT_TRUE(supplied.verify());
}

TEST(inlet_function, function_value_alternative_rejects_invalid_declaration)
{
  auto inlet = createBasicInlet("");
  axom::slic::ScopedAbortToThrow abortGuard;

  EXPECT_THROW(inlet.addFunctionAsValueAlternative("", FunctionTag::Double, {}),
               axom::slic::SlicAbortException);
  EXPECT_THROW(inlet.addFunctionAsValueAlternative("value", FunctionTag::Void, {}),
               axom::slic::SlicAbortException);
  EXPECT_TRUE(inlet.getGlobalContainer().getFunctionValueAlternativeNames().empty());
}

TEST(inlet_function, function_value_alternative_rejects_declaration_after_the_value)
{
  // Declaring the concrete entry first cannot work: it has already been read and
  // marked as being of the wrong type. Without this check the schema author sees
  // a verification failure blaming the input instead of the schema.
  axom::slic::ScopedAbortToThrow abortGuard;

  auto scalar = createBasicInlet("scale = function() return 2.0 end");
  scalar.addDouble("scale");
  EXPECT_THROW(scalar.addFunctionAsValueAlternative("scale", FunctionTag::Double, {}),
               axom::slic::SlicAbortException);

  auto collection = createBasicInlet("values = function() return {1.0, 2.0} end");
  collection.addDoubleArray("values");
  EXPECT_THROW(collection.addFunctionAsValueAlternative("values", FunctionTag::Vector, {}),
               axom::slic::SlicAbortException);

  // Declaring through a parent Container is rejected the same way
  auto nested = createBasicInlet("group = { value = function() return 2.0 end }");
  nested.addStruct("group").addDouble("value");
  EXPECT_THROW(nested.addFunctionAsValueAlternative("group/value", FunctionTag::Double, {}),
               axom::slic::SlicAbortException);

  // Collection declarations must check the concrete entries in their elements.
  auto structArray = createBasicInlet("items = {{value = function() return 2.0 end}}");
  auto& arrayEntries = structArray.addStructArray("items");
  arrayEntries.addDouble("value");
  EXPECT_THROW(arrayEntries.addFunctionAsValueAlternative("value", FunctionTag::Double, {}),
               axom::slic::SlicAbortException);

  auto structDictionary =
    createBasicInlet("items = {first = {values = function() return {1.0, 2.0} end}}");
  auto& dictionaryEntries = structDictionary.addStructDictionary("items");
  dictionaryEntries.addDoubleArray("values");
  EXPECT_THROW(dictionaryEntries.addFunctionAsValueAlternative("values", FunctionTag::Vector, {}),
               axom::slic::SlicAbortException);
}

TEST(inlet_function, function_value_alternative_applies_to_every_primitive_type)
{
  // Each primitive type reads its concrete entry separately, so check that a
  // supplied function is treated as an absent value rather than a wrong type
  auto inlet = createBasicInlet(R"(
    b = function() return 'true' end
    i = function() return 2.0 end
    d = function() return 3.0 end
    s = function() return 'four' end
  )");
  inlet.getGlobalContainer().strict();

  inlet.addFunctionAsValueAlternative("b", FunctionTag::String, {});
  inlet.addBool("b");
  inlet.addFunctionAsValueAlternative("i", FunctionTag::Double, {});
  inlet.addInt("i");
  inlet.addFunctionAsValueAlternative("d", FunctionTag::Double, {});
  inlet.addDouble("d");
  inlet.addFunctionAsValueAlternative("s", FunctionTag::String, {});
  inlet.addString("s");

  EXPECT_TRUE(inlet.verify());
  EXPECT_TRUE(inlet.unexpectedNames().empty());
  const auto& root = inlet.getGlobalContainer();
  EXPECT_EQ(root.getFunctionValueAlternativeNames(),
            (std::vector<std::string> {"b", "d", "i", "s"}));
  for(const std::string name : {"b", "i", "d", "s"})
  {
    EXPECT_FALSE(inlet.contains(name)) << name;
    EXPECT_TRUE(inlet.isUserProvided(name)) << name;
  }
  EXPECT_EQ(root.getFunctionValueAlternative("b").call<std::string>(), "true");
  EXPECT_DOUBLE_EQ(root.getFunctionValueAlternative("i").call<double>(), 2.0);
  EXPECT_DOUBLE_EQ(root.getFunctionValueAlternative("d").call<double>(), 3.0);
  EXPECT_EQ(root.getFunctionValueAlternative("s").call<std::string>(), "four");

  // Without an alternative, the same input is the wrong type for each primitive
  for(const std::string name : {"b", "i", "d", "s"})
  {
    auto plain = createBasicInlet(R"(
      b = function() return 'true' end
      i = function() return 2.0 end
      d = function() return 3.0 end
      s = function() return 'four' end
    )");
    if(name == "b")
    {
      plain.addBool(name);
    }
    else if(name == "i")
    {
      plain.addInt(name);
    }
    else if(name == "d")
    {
      plain.addDouble(name);
    }
    else
    {
      plain.addString(name);
    }
    EXPECT_FALSE(plain.verify()) << name;
  }
}

TEST(inlet_function, function_value_alternative_coexists_with_a_default)
{
  // A default populates the concrete entry even when the input supplies a function,
  // so check for the function first
  auto callback = createBasicInlet("scale = function() return 3.0 end");
  callback.addFunctionAsValueAlternative("scale", FunctionTag::Double, {});
  callback.addDouble("scale").defaultValue(1.0);
  EXPECT_TRUE(callback.verify());
  EXPECT_TRUE(callback.contains("scale"));
  EXPECT_TRUE(callback.isUserProvided("scale"));
  EXPECT_DOUBLE_EQ(callback["scale"].get<double>(), 1.0);
  ASSERT_TRUE(callback.getGlobalContainer().containsFunctionValueAlternative("scale"));
  EXPECT_DOUBLE_EQ(callback.getGlobalContainer().getFunctionValueAlternative("scale").call<double>(),
                   3.0);

  auto concrete = createBasicInlet("scale = 2.0");
  concrete.addFunctionAsValueAlternative("scale", FunctionTag::Double, {});
  concrete.addDouble("scale").defaultValue(1.0);
  EXPECT_TRUE(concrete.verify());
  EXPECT_FALSE(concrete.getGlobalContainer().containsFunctionValueAlternative("scale"));
  EXPECT_DOUBLE_EQ(concrete["scale"].get<double>(), 2.0);

  auto absent = createBasicInlet("");
  absent.addFunctionAsValueAlternative("scale", FunctionTag::Double, {});
  absent.addDouble("scale").defaultValue(1.0);
  EXPECT_TRUE(absent.verify());
  EXPECT_FALSE(absent.isUserProvided("scale"));
  EXPECT_FALSE(absent.getGlobalContainer().containsFunctionValueAlternative("scale"));
  EXPECT_DOUBLE_EQ(absent["scale"].get<double>(), 1.0);
}

TEST(inlet_function, function_value_alternative_constraints_apply_to_each_form_separately)
{
  // required() on the concrete entry is not satisfied by a supplied function
  auto required = createBasicInlet("scale = function() return 2.0 end");
  required.addFunctionAsValueAlternative("scale", FunctionTag::Double, {});
  required.addDouble("scale").required();
  EXPECT_FALSE(required.verify());

  // A range on the concrete entry does not apply to the function's result
  auto ranged = createBasicInlet("scale = function() return 5.0 end");
  ranged.addFunctionAsValueAlternative("scale", FunctionTag::Double, {});
  ranged.addDouble("scale").range(0.0, 1.0);
  EXPECT_TRUE(ranged.verify());

  // ...but still applies to a concrete value
  auto concrete = createBasicInlet("scale = 5.0");
  concrete.addFunctionAsValueAlternative("scale", FunctionTag::Double, {});
  concrete.addDouble("scale").range(0.0, 1.0);
  EXPECT_FALSE(concrete.verify());
}

TEST(inlet_function, function_value_alternative_lookup_of_absent_function)
{
  auto inlet = createBasicInlet("scale = 2.0");
  inlet.addFunctionAsValueAlternative("scale", FunctionTag::Double, {});
  inlet.addDouble("scale");
  ASSERT_TRUE(inlet.verify());
  const auto& root = inlet.getGlobalContainer();

  // Declared but not supplied: the Function exists but holds no callable
  EXPECT_FALSE(root.containsFunctionValueAlternative("scale"));
  EXPECT_FALSE(static_cast<bool>(root.getFunctionValueAlternative("scale")));

  // Never declared
  axom::slic::ScopedAbortToThrow abortGuard;
  EXPECT_FALSE(root.containsFunctionValueAlternative("offset"));
  EXPECT_THROW(root.getFunctionValueAlternative("offset"), axom::slic::SlicAbortException);
}

TEST(inlet_function, or_function_selects_the_supplied_form)
{
  auto inlet = createBasicInlet(R"(
    label = function() return 'computed' end
    scale = 4.0
  )");
  inlet.addString("label").orFunction();
  inlet.addDouble("scale").orFunction();
  EXPECT_TRUE(inlet.verify());

  // type() tells the caller which form the input supplied
  EXPECT_EQ(inlet["label"].type(), InletType::Function);
  EXPECT_EQ(inlet["scale"].type(), InletType::Double);

  // evaluate() handles either form
  EXPECT_EQ(inlet["label"].evaluate<std::string>(), "computed");
  EXPECT_DOUBLE_EQ(inlet["scale"].evaluate<double>(), 4.0);

  // The concrete form is still an ordinary value
  EXPECT_DOUBLE_EQ(inlet["scale"].get<double>(), 4.0);
}

TEST(inlet_function, or_function_evaluate_passes_arguments_to_the_function_only)
{
  auto inlet = createBasicInlet(R"(
    doubled = function(x) return 2 * x end
    constant = 4.0
    greeting = function(name) return 'hello ' .. name end
  )");
  inlet.addDouble("doubled").orFunction({FunctionTag::Double});
  inlet.addDouble("constant").orFunction({FunctionTag::Double});
  inlet.addString("greeting").orFunction({FunctionTag::String});
  EXPECT_TRUE(inlet.verify());

  // Numbers are passed as double, so an int literal or an lvalue also matches
  const double x = 3.0;
  EXPECT_DOUBLE_EQ(inlet["doubled"].evaluate<double>(3.0), 6.0);
  EXPECT_DOUBLE_EQ(inlet["doubled"].evaluate<double>(3), 6.0);
  EXPECT_DOUBLE_EQ(inlet["doubled"].evaluate<double>(x), 6.0);

  // A concrete value acts as a constant function of the arguments
  EXPECT_DOUBLE_EQ(inlet["constant"].evaluate<double>(3.0), 4.0);

  // Character strings are passed as std::string
  EXPECT_EQ(inlet["greeting"].evaluate<std::string>("inlet"), "hello inlet");
}

TEST(inlet_function, or_function_entry_is_present_in_either_form)
{
  auto inlet = createBasicInlet(R"(
    from_function = function() return 2.0 end
    from_value = 3.0
  )");
  inlet.addDouble("from_function").orFunction();
  inlet.addDouble("from_value").orFunction();
  inlet.addDouble("absent").orFunction();
  EXPECT_TRUE(inlet.verify());

  for(const std::string name : {"from_function", "from_value"})
  {
    EXPECT_TRUE(inlet.contains(name)) << name;
    EXPECT_TRUE(inlet.isUserProvided(name)) << name;
  }
  EXPECT_FALSE(inlet.contains("absent"));
  EXPECT_FALSE(inlet.isUserProvided("absent"));

  const auto& fields = inlet.getGlobalContainer().getChildFields();
  EXPECT_TRUE(fields.at("from_function")->isFunction());
  EXPECT_FALSE(fields.at("from_value")->isFunction());
  EXPECT_FALSE(fields.at("absent")->isFunction());
}

TEST(inlet_function, or_function_satisfies_required_in_either_form)
{
  for(const std::string input : {"scale = function() return 2.0 end", "scale = 2.0"})
  {
    auto inlet = createBasicInlet(input);
    inlet.addDouble("scale").orFunction().required();
    EXPECT_TRUE(inlet.verify()) << input;
  }

  auto missing = createBasicInlet("");
  missing.addDouble("scale").orFunction().required();
  EXPECT_FALSE(missing.verify());
}

TEST(inlet_function, or_function_default_applies_only_when_neither_form_is_supplied)
{
  // The default may be declared before or after orFunction()
  for(const bool defaultFirst : {true, false})
  {
    auto declare = [defaultFirst](Inlet& inlet) {
      auto& scale = inlet.addDouble("scale");
      if(defaultFirst)
      {
        scale.defaultValue(1.0).orFunction();
      }
      else
      {
        scale.orFunction().defaultValue(1.0);
      }
    };

    auto callback = createBasicInlet("scale = function() return 3.0 end");
    declare(callback);
    EXPECT_TRUE(callback.verify());
    EXPECT_EQ(callback["scale"].type(), InletType::Function) << defaultFirst;
    EXPECT_DOUBLE_EQ(callback["scale"].evaluate<double>(), 3.0) << defaultFirst;
    {
      // The default does not linger behind the function
      axom::slic::ScopedAbortToThrow abortGuard;
      EXPECT_THROW(callback.getGlobalContainer().getChildFields().at("scale")->get<double>(),
                   axom::slic::SlicAbortException)
        << defaultFirst;
      EXPECT_TRUE(slicErrorContains("was supplied as a function"));
    }

    auto concrete = createBasicInlet("scale = 2.0");
    declare(concrete);
    EXPECT_TRUE(concrete.verify());
    EXPECT_DOUBLE_EQ(concrete["scale"].evaluate<double>(), 2.0) << defaultFirst;

    auto absent = createBasicInlet("");
    declare(absent);
    EXPECT_TRUE(absent.verify());
    EXPECT_FALSE(absent.isUserProvided("scale"));
    EXPECT_DOUBLE_EQ(absent["scale"].evaluate<double>(), 1.0) << defaultFirst;
  }
}

TEST(inlet_function, or_function_reading_a_supplied_function_as_a_value_is_an_error)
{
  auto inlet = createBasicInlet("scale = function() return 2.0 end");
  inlet.addDouble("scale").orFunction();
  EXPECT_TRUE(inlet.verify());

  axom::slic::ScopedAbortToThrow abortGuard;
  EXPECT_THROW(inlet["scale"].get<double>(), axom::slic::SlicAbortException);
  EXPECT_TRUE(slicErrorContains("was supplied as a function"));
  EXPECT_THROW(inlet.getGlobalContainer().getChildFields().at("scale")->get<double>(),
               axom::slic::SlicAbortException);
  EXPECT_TRUE(slicErrorContains("was supplied as a function"));
}

TEST(inlet_function, or_function_rejects_an_unrelated_wrong_type)
{
  auto inlet = createBasicInlet("scale = 'neither a number nor a function'");
  inlet.addDouble("scale").orFunction();
  EXPECT_FALSE(inlet.verify());
  EXPECT_FALSE(inlet.contains("scale"));
  EXPECT_FALSE(inlet.getGlobalContainer().getChildFields().at("scale")->isFunction());
}

TEST(inlet_function, or_function_rejects_unsupported_field_types)
{
  // FunctionTag has no Bool or Int return type
  auto inlet = createBasicInlet("");
  axom::slic::ScopedAbortToThrow abortGuard;
  EXPECT_THROW(inlet.addBool("flag").orFunction(), axom::slic::SlicAbortException);
  EXPECT_TRUE(slicErrorContains("supports double and string Fields"));
  EXPECT_THROW(inlet.addInt("count").orFunction(), axom::slic::SlicAbortException);
  EXPECT_TRUE(slicErrorContains("supports double and string Fields"));
}

TEST(inlet_function, or_function_redeclaration_must_match_the_signature)
{
  auto inlet = createBasicInlet("group = { scale = function(x) return 2 * x end }");
  auto& group = inlet.addStruct("group");
  group.addDouble("scale").orFunction({FunctionTag::Double});

  // Same signature: no-op. Redeclaring the Field returns the existing Field
  group.addDouble("scale").orFunction({FunctionTag::Double});
  EXPECT_TRUE(inlet.verify());
  EXPECT_DOUBLE_EQ(inlet["group/scale"].evaluate<double>(2.0), 4.0);

  axom::slic::ScopedAbortToThrow abortGuard;
  EXPECT_THROW(group.addDouble("scale").orFunction(), axom::slic::SlicAbortException);
  EXPECT_TRUE(slicErrorContains("different signature"));
}

TEST(inlet_function, or_function_is_recognized_by_strict_containers)
{
  auto inlet = createBasicInlet(R"(
    group = { scale = function() return 2.0 end, offset = 1.0 }
  )");
  auto& group = inlet.addStruct("group").strict();
  group.addDouble("scale").orFunction();
  group.addDouble("offset").orFunction();
  EXPECT_TRUE(inlet.verify());
  EXPECT_TRUE(inlet.unexpectedNames().empty());
  EXPECT_TRUE(group.exists());
  EXPECT_TRUE(group.isUserProvided());
}

TEST(inlet_function, or_function_constraints_apply_to_the_value_only)
{
  // Inlet does not call the function during verification
  auto callback = createBasicInlet("scale = function() return 5.0 end");
  callback.addDouble("scale").range(0.0, 1.0).orFunction();
  EXPECT_TRUE(callback.verify());

  auto concrete = createBasicInlet("scale = 5.0");
  concrete.addDouble("scale").range(0.0, 1.0).orFunction();
  EXPECT_FALSE(concrete.verify());
}

TEST(inlet_function, or_function_with_a_reader_that_has_no_functions)
{
  // A reader without function support only ever supplies the value
  auto reader = std::make_unique<axom::inlet::YAMLReader>();
  reader->parseString("scale: 2.5");
  Inlet inlet(std::move(reader));
  inlet.addDouble("scale").orFunction();
  EXPECT_TRUE(inlet.verify());
  EXPECT_EQ(inlet["scale"].type(), InletType::Double);
  EXPECT_DOUBLE_EQ(inlet["scale"].evaluate<double>(), 2.5);
}

TEST(inlet_function, simple_void_to_double_through_container)
{
  std::string testString = "function foo () return 9.64 end";
  auto inlet = createBasicInlet(testString);

  inlet.addFunction("foo", FunctionTag::Double, {}, "foo's description");

  auto callable = inlet["foo"].get<std::function<FunctionType::Double()>>();
  double result = callable();
  EXPECT_FLOAT_EQ(result, 9.64);
}

TEST(inlet_function, simple_double_to_void_through_container)
{
  // Test a function that returns nothing by using it to modify a global
  std::string testString = "bar = 19.9; function foo (a) bar = a end";
  auto inlet = createBasicInlet(testString);

  inlet.addFunction("foo", FunctionTag::Void, {FunctionTag::Double}, "foo's description");

  auto callable = inlet["foo"].get<std::function<FunctionType::Void(FunctionType::Double)>>();
  double arg = -6.37;
  callable(arg);

  inlet.addDouble("bar", "bar's description");
  double result = inlet["bar"];
  EXPECT_FLOAT_EQ(result, arg);
}

TEST(inlet_function, simple_string_to_double_through_container)
{
  std::string testString =
    "function foo(s) "
    "  if s == 'a' then return 9.1 "
    "  elseif s == 'b' then return -6.3 "
    "  else return 66.5 end "
    "end";
  auto inlet = createBasicInlet(testString);

  inlet.addFunction("foo", FunctionTag::Double, {FunctionTag::String}, "foo's description");

  auto callable = inlet["foo"].get<std::function<FunctionType::Double(FunctionType::String)>>();

  EXPECT_FLOAT_EQ(callable("a"), 9.1);
  EXPECT_FLOAT_EQ(callable("b"), -6.3);
  EXPECT_FLOAT_EQ(callable("c"), 66.5);
}

TEST(inlet_function, simple_double_to_string_through_container)
{
  std::string testString =
    "function foo(d) "
    "  if d == 1 then return 'a' "
    "  elseif d == 2 then return 'b' "
    "  else return 'c' end "
    "end";
  auto inlet = createBasicInlet(testString);

  inlet.addFunction("foo", FunctionTag::String, {FunctionTag::Double}, "foo's description");

  auto callable = inlet["foo"].get<std::function<FunctionType::String(FunctionType::Double)>>();
  EXPECT_EQ(callable(1), "a");
  EXPECT_EQ(callable(2), "b");
  EXPECT_EQ(callable(3), "c");
}

TEST(inlet_function, simple_vec3_to_double_through_container_call)
{
  std::string testString = "function foo (v) return v.x + v.y + v.z end";
  auto inlet = createBasicInlet(testString);

  inlet.addFunction("foo", FunctionTag::Double, {FunctionTag::Vector}, "foo's description");

  auto result = inlet["foo"].call<double>(FunctionType::Vector {1, 2, 3});
  EXPECT_FLOAT_EQ(result, 6);
}

TEST(inlet_function, simple_vec3_to_vec3_through_container_call)
{
  std::string testString = "function foo (v) return 2*v end";
  auto inlet = createBasicInlet(testString);

  inlet.addFunction("foo", FunctionTag::Vector, {FunctionTag::Vector}, "foo's description");

  auto result = inlet["foo"].call<FunctionType::Vector>(FunctionType::Vector {1, 2, 3});
  EXPECT_FLOAT_EQ(result[0], 2);
  EXPECT_FLOAT_EQ(result[1], 4);
  EXPECT_FLOAT_EQ(result[2], 6);
}

TEST(inlet_function, simple_vec3_double_to_double_through_container_call)
{
  std::string testString = "function foo (v, t) return t * (v.x + v.y + v.z) end";
  auto inlet = createBasicInlet(testString);

  inlet.addFunction("foo",
                    FunctionTag::Double,
                    {FunctionTag::Vector, FunctionTag::Double},
                    "foo's description");

  auto result = inlet["foo"].call<double>(FunctionType::Vector {1, 2, 3}, 2.0);
  EXPECT_FLOAT_EQ(result, 12);
}

TEST(inlet_function, simple_vec3_double_to_vec3_through_container_call)
{
  std::string testString = "function foo (v, t) return t*v end";
  auto inlet = createBasicInlet(testString);

  inlet.addFunction("foo",
                    FunctionTag::Vector,
                    {FunctionTag::Vector, FunctionTag::Double},
                    "foo's description");

  auto result = inlet["foo"].call<FunctionType::Vector>(FunctionType::Vector {1, 2, 3}, 2.0);
  EXPECT_FLOAT_EQ(result[0], 2);
  EXPECT_FLOAT_EQ(result[1], 4);
  EXPECT_FLOAT_EQ(result[2], 6);
}

TEST(inlet_function, simple_vec3_to_vec3_verify_lambda_pass)
{
  std::string testString = "function foo (v) return 2*v end";
  auto inlet = createBasicInlet(testString);

  auto& func =
    inlet.addFunction("foo", FunctionTag::Vector, {FunctionTag::Vector}, "foo's description").required();
  func.registerVerifier([](const axom::inlet::Function& func) {
    auto result = func.call<FunctionType::Vector>(FunctionType::Vector {1, 0, 0});
    return std::abs(result[0] - 2) < 1e-5;
  });

  EXPECT_TRUE(inlet.verify());
}

TEST(inlet_function, simple_vec3_to_vec3_verify_lambda_fail)
{
  std::string testString = "function foo (v) return 2*v end";
  auto inlet = createBasicInlet(testString);

  auto& func =
    inlet.addFunction("foo", FunctionTag::Vector, {FunctionTag::Vector}, "foo's description").required();
  func.registerVerifier([](const axom::inlet::Function& func) {
    auto result = func.call<FunctionType::Vector>(FunctionType::Vector {2, 0, 0});
    return std::abs(result[0] - 2) < 1e-5;
  });

  EXPECT_FALSE(inlet.verify());
}

TEST(inlet_function, simple_vec3_to_vec3_verify_lambda_with_errors_fail)
{
  std::string testString = "function foo (v) return 2*v end";
  auto inlet = createBasicInlet(testString);

  auto& func =
    inlet.addFunction("foo", FunctionTag::Vector, {FunctionTag::Vector}, "foo's description").required();
  func.registerVerifier([](const axom::inlet::Function& func, std::vector<VerificationError>* errors) {
    INLET_VERIFICATION_WARNING("foo", "Something bad happened", errors);
    auto result = func.call<FunctionType::Vector>(FunctionType::Vector {2, 0, 0});
    return std::abs(result[0] - 2) < 1e-5;
  });

  std::vector<VerificationError> errors;
  EXPECT_FALSE(inlet.verify(&errors));
  ASSERT_FALSE(errors.empty());
  ASSERT_EQ(axom::Path("foo"), errors[0].path);
  ASSERT_EQ("Something bad happened", errors[0].message);
}

struct Foo
{
  bool bar;
  std::function<FunctionType::Vector(FunctionType::Vector)> baz;
};

template <>
struct FromInlet<Foo>
{
  Foo operator()(const axom::inlet::Container& base)
  {
    Foo f {base["bar"], base["baz"]};
    return f;
  }
};

struct FooWithValueAlternative
{
  double bar;
};

template <>
struct FromInlet<FooWithValueAlternative>
{
  FooWithValueAlternative operator()(const axom::inlet::Container& base)
  {
    if(base.containsFunctionValueAlternative("bar"))
    {
      return {base.getFunctionValueAlternative("bar").call<double>()};
    }
    return {base["bar"].get<double>()};
  }
};

struct FooWithValueAlternativeDictionary
{
  std::unordered_map<std::string, FooWithValueAlternative> values;
};

template <>
struct FromInlet<FooWithValueAlternativeDictionary>
{
  FooWithValueAlternativeDictionary operator()(const axom::inlet::Container& base)
  {
    return {base["foo"].get<std::unordered_map<std::string, FooWithValueAlternative>>()};
  }
};

TEST(inlet_function, simple_vec3_to_vec3_struct)
{
  std::string testString = "foo = { bar = true; baz = function (v) return 2*v end }";
  auto inlet = createBasicInlet(testString);

  // Define schema
  inlet.addBool("foo/bar", "bar's description");
  inlet.addFunction("foo/baz", FunctionTag::Vector, {FunctionTag::Vector}, "baz's description")
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
  auto inlet = createBasicInlet(testString);

  auto& arr_container = inlet.addStructArray("foo");

  // Define schema
  arr_container.addBool("bar", "bar's description");
  arr_container.addFunction("baz", FunctionTag::Vector, {FunctionTag::Vector}, "baz's description")
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

TEST(inlet_function, function_value_alternative_in_nested_dictionary_of_struct)
{
  // Both representations appear in the same collection, so the alternative must
  // be expanded across the collection's elements exactly as the concrete entry is
  auto inlet = createBasicInlet(R"(
    groups = {
      [0] = {
        foo = {
          first = {bar = 2},
          second = {bar = function () return 3 end}
        }
      }
    }
  )");

  auto& groups = inlet.addStructArray("groups");
  auto& foos = groups.addStructDictionary("foo");
  foos.addFunctionAsValueAlternative("bar", FunctionTag::Double, {});
  foos.addDouble("bar");

  EXPECT_TRUE(inlet.verify());
  const auto values =
    inlet["groups"].get<std::unordered_map<int, FooWithValueAlternativeDictionary>>();
  EXPECT_DOUBLE_EQ(values.at(0).values.at("first").bar, 2.0);
  EXPECT_DOUBLE_EQ(values.at(0).values.at("second").bar, 3.0);
}

struct OperatorWithOrFunction
{
  double rotate;
};

template <>
struct FromInlet<OperatorWithOrFunction>
{
  OperatorWithOrFunction operator()(const axom::inlet::Container& base)
  {
    return {base["rotate"].evaluate<double>()};
  }
};

struct OperatorsWithOrFunctionDictionary
{
  std::unordered_map<std::string, OperatorWithOrFunction> operators;
};

template <>
struct FromInlet<OperatorsWithOrFunctionDictionary>
{
  OperatorsWithOrFunctionDictionary operator()(const axom::inlet::Container& base)
  {
    return {base["ops"].get<std::unordered_map<std::string, OperatorWithOrFunction>>()};
  }
};

TEST(inlet_function, or_function_applies_to_every_element_of_a_struct_array)
{
  // Elements may use different forms
  auto inlet = createBasicInlet(R"(
    operators = {
      { rotate = 30.0 },
      { rotate = function() return 45.0 end },
      { }
    }
  )");
  auto& ops = inlet.addStructArray("operators").strict();
  ops.addDouble("rotate").orFunction().defaultValue(0.0);
  EXPECT_TRUE(inlet.verify());
  EXPECT_TRUE(inlet.unexpectedNames().empty());

  const auto values = inlet["operators"].get<std::vector<OperatorWithOrFunction>>();
  ASSERT_EQ(values.size(), 3u);
  EXPECT_DOUBLE_EQ(values[0].rotate, 30.0);
  EXPECT_DOUBLE_EQ(values[1].rotate, 45.0);
  EXPECT_DOUBLE_EQ(values[2].rotate, 0.0);
}

TEST(inlet_function, or_function_in_nested_dictionary_of_struct)
{
  auto inlet = createBasicInlet(R"(
    groups = {
      [0] = {
        ops = {
          first = {rotate = 2},
          second = {rotate = function () return 3 end}
        }
      }
    }
  )");

  auto& groups = inlet.addStructArray("groups");
  auto& ops = groups.addStructDictionary("ops");
  ops.addDouble("rotate").orFunction().required();

  EXPECT_TRUE(inlet.verify());
  const auto values =
    inlet["groups"].get<std::unordered_map<int, OperatorsWithOrFunctionDictionary>>();
  EXPECT_DOUBLE_EQ(values.at(0).operators.at("first").rotate, 2.0);
  EXPECT_DOUBLE_EQ(values.at(0).operators.at("second").rotate, 3.0);
}

TEST(inlet_function, dimension_dependent_result)
{
  std::string testString =
    "foo = function (v) "
    "first = 2 * v.x "
    "last = 2 * v.y "
    "if v.dim == 2 then "
    "return Vector.new(first, last) "
    "else "
    "return Vector.new(first, 0, last) "
    "end "
    "end";
  auto inlet = createBasicInlet(testString);

  inlet.addFunction("foo", FunctionTag::Vector, {FunctionTag::Vector}, "foo's description");

  auto callable = inlet["foo"].get<std::function<FunctionType::Vector(FunctionType::Vector)>>();

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

struct FooWithScalarFunc
{
  std::function<double(double)> bar;
};

template <>
struct FromInlet<FooWithScalarFunc>
{
  FooWithScalarFunc operator()(const axom::inlet::Container& base) { return {base["foo/bar"]}; }
};

TEST(inlet_function, nested_function_in_struct)
{
  std::string testString =
    "quux = { [0] = { foo = { bar = function (x) return x + 1 end } }, "
    "         [1] = { foo = { bar = function (x) return x + 3 end } } }";
  auto inlet = createBasicInlet(testString);

  auto& quux_schema = inlet.addStructArray("quux");
  auto& foo_schema = quux_schema.addStruct("foo");

  foo_schema.addFunction("bar", FunctionTag::Double, {FunctionTag::Double}, "bar's description");

  auto foos = inlet["quux"].get<std::vector<FooWithScalarFunc>>();
  EXPECT_EQ(foos.size(), 2);

  auto& first_func = foos[0].bar;
  // Check that the function object contains a valid target
  EXPECT_TRUE(static_cast<bool>(first_func));
  EXPECT_DOUBLE_EQ(first_func(4.0), 5.0);

  auto& second_func = foos[1].bar;
  // Check that the function object contains a valid target
  EXPECT_TRUE(static_cast<bool>(second_func));
  EXPECT_DOUBLE_EQ(second_func(4.0), 7.0);
}

template <typename Ret, typename... Args>
Ret checkedCall(const axom::sol::protected_function& func, Args&&... args)
{
  auto tentative_result = func(std::forward<Args>(args)...);
  EXPECT_TRUE(tentative_result.valid());
  return tentative_result;
}

/**
 * The inlet_function_usertype suite is intended to verify the correctness of the
 * definition of the correspondence between the FunctionType::Vector type and its
 * lua usertype equivalent.  Instead of using the Inlet interface to define and
 * access functions, the axom::sol::state member is interrogated directly via the
 * derived class SolStateReader to avoid mixing concerns in these tests.
 * 
 * Each entry in the Lua table/metatable for this usertype has a corresponding
 * test, i.e., one for each operator overload/constructor/member variable.
 */
TEST(inlet_function_usertype, lua_usertype_basic)
{
  std::string testString = "function func(vec) return 7 end";
  SolStateReader reader;
  reader.parseString(testString);
  axom::sol::protected_function func = (*reader.solState())["func"];
  axom::inlet::FunctionType::Vector vec {1, 2, 3};
  int result = checkedCall<int>(func, vec);
  EXPECT_EQ(result, 7);
}

TEST(inlet_function_usertype, lua_usertype_basic_ret)
{
  std::string testString = "function func(x, y, z) return Vector.new(x, y, z) end";
  SolStateReader reader;
  reader.parseString(testString);
  axom::sol::protected_function func = (*reader.solState())["func"];
  axom::inlet::FunctionType::Vector vec {1, 2, 3};
  auto result = checkedCall<axom::inlet::FunctionType::Vector>(func, 1, 2, 3);
  EXPECT_EQ(vec, result);
}

TEST(inlet_function_usertype, lua_usertype_basic_ret_2d)
{
  std::string testString = "function func(x, y, z) return Vector.new(x, y) end";
  SolStateReader reader;
  reader.parseString(testString);
  axom::sol::protected_function func = (*reader.solState())["func"];
  axom::inlet::FunctionType::Vector vec {1, 2};
  auto result = checkedCall<axom::inlet::FunctionType::Vector>(func, 1, 2, 3);
  EXPECT_EQ(vec, result);
}

TEST(inlet_function_usertype, lua_usertype_basic_ret_default)
{
  std::string testString = "function func(x, y, z) return Vector.new() end";
  SolStateReader reader;
  reader.parseString(testString);
  axom::sol::protected_function func = (*reader.solState())["func"];
  axom::inlet::FunctionType::Vector vec {0, 0, 0};
  auto result = checkedCall<axom::inlet::FunctionType::Vector>(func, 1, 2, 3);
  EXPECT_EQ(vec, result);
}

TEST(inlet_function_usertype, lua_usertype_basic_add)
{
  std::string testString = "function func(vec1, vec2) return vec1 + vec2 end";
  SolStateReader reader;
  reader.parseString(testString);
  axom::sol::protected_function func = (*reader.solState())["func"];
  axom::inlet::FunctionType::Vector vec1 {1, 2, 3};
  axom::inlet::FunctionType::Vector vec2 {4, 5, 6};
  const axom::inlet::FunctionType::Vector sum {5, 7, 9};
  auto result = checkedCall<axom::inlet::FunctionType::Vector>(func, vec1, vec2);
  EXPECT_EQ(result, sum);
}

TEST(inlet_function_usertype, lua_usertype_basic_sub)
{
  std::string testString = "function func(vec1, vec2) return vec1 - vec2 end";
  SolStateReader reader;
  reader.parseString(testString);
  axom::sol::protected_function func = (*reader.solState())["func"];
  axom::inlet::FunctionType::Vector vec1 {1, 2, 3};
  axom::inlet::FunctionType::Vector vec2 {4, 5, 6};
  const axom::inlet::FunctionType::Vector difference {-3, -3, -3};
  auto result = checkedCall<axom::inlet::FunctionType::Vector>(func, vec1, vec2);
  EXPECT_EQ(result, difference);
}

TEST(inlet_function_usertype, lua_usertype_basic_negate)
{
  std::string testString = "function func(vec) return -vec end";
  SolStateReader reader;
  reader.parseString(testString);
  axom::sol::protected_function func = (*reader.solState())["func"];
  axom::inlet::FunctionType::Vector vec {1, 2, 3};
  const axom::inlet::FunctionType::Vector negated {-1, -2, -3};
  auto result = checkedCall<axom::inlet::FunctionType::Vector>(func, vec);
  EXPECT_EQ(result, negated);
}

TEST(inlet_function_usertype, lua_usertype_basic_scalar_mult)
{
  std::string testString =
    "function func1(vec, x) return vec * x end; function func2(vec, x) return "
    "x * vec end";
  SolStateReader reader;
  reader.parseString(testString);
  axom::sol::protected_function func1 = (*reader.solState())["func1"];
  axom::sol::protected_function func2 = (*reader.solState())["func2"];
  axom::inlet::FunctionType::Vector vec {1, 2, 3};
  const axom::inlet::FunctionType::Vector doubled {2, 4, 6};
  auto result = checkedCall<axom::inlet::FunctionType::Vector>(func1, vec, 2.0);
  EXPECT_EQ(result, doubled);
  const axom::inlet::FunctionType::Vector tripled {3, 6, 9};
  result = checkedCall<axom::inlet::FunctionType::Vector>(func2, vec, 3.0);
  EXPECT_EQ(result, tripled);
}

TEST(inlet_function_usertype, lua_usertype_basic_index_get)
{
  std::string testString = "function func(vec, idx) return vec[idx] end";
  SolStateReader reader;
  reader.parseString(testString);
  axom::sol::protected_function func = (*reader.solState())["func"];
  axom::inlet::FunctionType::Vector vec {1, 2, 3};
  // Use 1-based indexing in these tests as lua is 1-indexed
  auto result = checkedCall<double>(func, vec, 1);
  EXPECT_FLOAT_EQ(1, result);
  result = checkedCall<double>(func, vec, 2);
  EXPECT_FLOAT_EQ(2, result);
  result = checkedCall<double>(func, vec, 3);
  EXPECT_FLOAT_EQ(3, result);
}

TEST(inlet_function_usertype, lua_usertype_basic_index_set)
{
  std::string testString =
    "function func(idx) vec = Vector.new(1,1,1); vec[idx] = -1; return vec end";
  SolStateReader reader;
  reader.parseString(testString);
  axom::sol::protected_function func = (*reader.solState())["func"];
  auto result = checkedCall<axom::inlet::FunctionType::Vector>(func, 1);
  EXPECT_FLOAT_EQ(-1, result[0]);
  result = checkedCall<axom::inlet::FunctionType::Vector>(func, 2);
  EXPECT_FLOAT_EQ(-1, result[1]);
  result = checkedCall<axom::inlet::FunctionType::Vector>(func, 3);
  EXPECT_FLOAT_EQ(-1, result[2]);
}

TEST(inlet_function_usertype, lua_usertype_basic_norm)
{
  std::string testString = "function func(vec) return vec:norm() end";
  SolStateReader reader;
  reader.parseString(testString);
  axom::sol::protected_function func = (*reader.solState())["func"];
  axom::inlet::FunctionType::Vector vec {1, 2, 3};
  const double l2_norm = std::sqrt((1 * 1) + (2 * 2) + (3 * 3));
  auto result = checkedCall<double>(func, vec);
  EXPECT_FLOAT_EQ(l2_norm, result);
}

TEST(inlet_function_usertype, lua_usertype_basic_squared_norm)
{
  std::string testString = "function func(vec) return vec:squared_norm() end";
  SolStateReader reader;
  reader.parseString(testString);
  axom::sol::protected_function func = (*reader.solState())["func"];
  axom::inlet::FunctionType::Vector vec {1, 2, 3};
  const double squared_l2_norm = (1 * 1) + (2 * 2) + (3 * 3);
  auto result = checkedCall<double>(func, vec);
  EXPECT_FLOAT_EQ(squared_l2_norm, result);
}

TEST(inlet_function_usertype, lua_usertype_basic_unit_vec)
{
  std::string testString = "function func(vec) return vec:unitVector() end";
  SolStateReader reader;
  reader.parseString(testString);
  axom::sol::protected_function func = (*reader.solState())["func"];
  axom::inlet::FunctionType::Vector vec {1, 2, 3};
  const double l2_norm = std::sqrt((1 * 1) + (2 * 2) + (3 * 3));
  const axom::inlet::FunctionType::Vector unit {1 / l2_norm, 2 / l2_norm, 3 / l2_norm};
  auto result = checkedCall<axom::inlet::FunctionType::Vector>(func, vec);
  EXPECT_EQ(unit, result);
}

TEST(inlet_function_usertype, lua_usertype_basic_dot)
{
  std::string testString = "function func(vec1, vec2) return vec1:dot(vec2) end";
  SolStateReader reader;
  reader.parseString(testString);
  axom::sol::protected_function func = (*reader.solState())["func"];
  axom::inlet::FunctionType::Vector vec1 {1, 2, 3};
  axom::inlet::FunctionType::Vector vec2 {4, 5, 6};
  const double dot = (1 * 4) + (2 * 5) + (3 * 6);
  auto result = checkedCall<double>(func, vec1, vec2);
  EXPECT_EQ(dot, result);
}

TEST(inlet_function_usertype, lua_usertype_basic_cross)
{
  std::string testString = "function func(vec1, vec2) return vec1:cross(vec2) end";
  SolStateReader reader;
  reader.parseString(testString);
  axom::sol::protected_function func = (*reader.solState())["func"];
  axom::inlet::FunctionType::Vector vec1 {1, 2, 3};
  axom::inlet::FunctionType::Vector vec2 {4, 5, 6};
  const double i = (2 * 6) - (3 * 5);
  const double j = (3 * 4) - (1 * 6);
  const double k = (1 * 5) - (2 * 4);
  const axom::inlet::FunctionType::Vector cross {i, j, k};
  auto result = checkedCall<axom::inlet::FunctionType::Vector>(func, vec1, vec2);
  EXPECT_EQ(cross, result);
}

TEST(inlet_function_usertype, lua_usertype_check_dim)
{
  std::string testString = "function func(vec) return vec.dim end";
  SolStateReader reader;
  reader.parseString(testString);
  axom::sol::protected_function func = (*reader.solState())["func"];
  axom::inlet::FunctionType::Vector vec1 {1, 2, 3};
  axom::inlet::FunctionType::Vector vec2 {4, 5};
  auto result = checkedCall<double>(func, vec1);
  EXPECT_EQ(result, 3);

  result = checkedCall<double>(func, vec2);
  EXPECT_EQ(result, 2);
}

TEST(inlet_function_usertype, lua_usertype_named_access)
{
  std::string testString =
    "function func(vec, comp) if comp == 1 then return vec.x elseif comp == 2 "
    "then return vec.y else return vec.z end end";
  SolStateReader reader;
  reader.parseString(testString);
  axom::sol::protected_function func = (*reader.solState())["func"];
  axom::inlet::FunctionType::Vector vec1 {4, 5, 6};
  auto result = checkedCall<double>(func, vec1, 1);
  EXPECT_EQ(result, 4);

  result = checkedCall<double>(func, vec1, 2);
  EXPECT_EQ(result, 5);

  result = checkedCall<double>(func, vec1, 3);
  EXPECT_EQ(result, 6);
}

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;
  axom::slic::addStreamToMsgLevel(new axom::slic::GenericOutputStream(&slicErrors(), "<MESSAGE>\n"),
                                  axom::slic::message::Error);

  result = RUN_ALL_TESTS();

  return result;
}

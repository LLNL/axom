// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/inlet/Inlet.hpp"
#include "axom/inlet/tests/inlet_test_utils.hpp"

using axom::inlet::Inlet;
using axom::inlet::VerificationError;
using axom::sidre::DataStore;

template <typename InletReader>
Inlet createBasicInlet(const std::string& luaString)
{
  std::unique_ptr<InletReader> reader(new InletReader());
  reader->parseString(axom::inlet::detail::fromLuaTo<InletReader>(luaString));
  const bool enableDocs = true;
  return Inlet(std::move(reader), enableDocs);
}

template <typename InletReader>
class inlet_errors : public ::testing::Test
{ };

TYPED_TEST_SUITE(inlet_errors, axom::inlet::detail::ReaderTypes);

TYPED_TEST(inlet_errors, required_field)
{
  std::string testString = "";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  inlet.addInt("foo", "foo's description").required();

  std::vector<VerificationError> errors;
  EXPECT_FALSE(inlet.verify(&errors));
  // Need something about "foo" and "Required"
  EXPECT_TRUE(
    std::any_of(errors.begin(), errors.end(), [](const VerificationError& err) {
      return (err.path == axom::Path {"foo"}) &&
        err.messageContains("Required");
    }));
}

TYPED_TEST(inlet_errors, wrong_type)
{
  std::string testString = "foo = 'hello'";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  inlet.addInt("foo", "foo's description");

  std::vector<VerificationError> errors;
  EXPECT_FALSE(inlet.verify(&errors));
  // Need something about "foo" and "wrong type"
  EXPECT_TRUE(
    std::any_of(errors.begin(), errors.end(), [](const VerificationError& err) {
      return (err.path == axom::Path {"foo"}) &&
        err.messageContains("wrong type");
    }));
}

TYPED_TEST(inlet_errors, wrong_type_nested)
{
  std::string testString = "foo = { bar = { baz = 'hello' } }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& foo = inlet.addStruct("foo", "foo's description");
  auto& bar = foo.addStruct("bar", "bar's description");
  bar.addInt("baz", "baz's description");

  std::vector<VerificationError> errors;
  EXPECT_FALSE(inlet.verify(&errors));
  // Need something about "foo/bar/baz" and "wrong type"
  EXPECT_TRUE(
    std::any_of(errors.begin(), errors.end(), [](const VerificationError& err) {
      return (err.path == axom::Path {"foo/bar/baz"}) &&
        err.messageContains("wrong type");
    }));
}

TYPED_TEST(inlet_errors, heterogeneous_array)
{
  std::string testString = "foo = { [0] = 1, [1] = 2, [2] = 'hello' }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  inlet.addIntArray("foo", "foo's description");

  std::vector<VerificationError> errors;
  EXPECT_FALSE(inlet.verify(&errors));
  // Need something about "foo" and "not homogeneous"
  EXPECT_TRUE(
    std::any_of(errors.begin(), errors.end(), [](const VerificationError& err) {
      // FIXME: Do we want to strip out the _inlet_collection from the error messages?
      return (err.path ==
              axom::inlet::appendPrefix(
                "foo",
                axom::inlet::detail::COLLECTION_GROUP_NAME)) &&
        err.messageContains("not homogeneous");
    }));
}

TYPED_TEST(inlet_errors, heterogeneous_array_nested)
{
  std::string testString =
    "foo = { bar = { baz = { [0] = 1, [1] = 2, [2] = 'hello' } } }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& foo = inlet.addStruct("foo", "foo's description");
  auto& bar = foo.addStruct("bar", "bar's description");
  bar.addIntArray("baz", "baz's description");

  std::vector<VerificationError> errors;
  EXPECT_FALSE(inlet.verify(&errors));
  // Need something about "foo/bar/baz" and "not homogeneous"
  EXPECT_TRUE(
    std::any_of(errors.begin(), errors.end(), [](const VerificationError& err) {
      return (err.path ==
              axom::inlet::appendPrefix(
                "foo/bar/baz",
                axom::inlet::detail::COLLECTION_GROUP_NAME)) &&
        err.messageContains("not homogeneous");
    }));
}

TYPED_TEST(inlet_errors, invalid_range)
{
  std::string testString = "foo = 7";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  inlet.addInt("foo", "foo's description").range(0, 5);

  std::vector<VerificationError> errors;
  EXPECT_FALSE(inlet.verify(&errors));
  // Need something about "foo" and "range"
  EXPECT_TRUE(
    std::any_of(errors.begin(), errors.end(), [](const VerificationError& err) {
      return (err.path == axom::Path {"foo"}) && err.messageContains("range");
    }));
}

TYPED_TEST(inlet_errors, invalid_range_nested)
{
  std::string testString = "foo = { bar = { baz = 7 } }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& foo = inlet.addStruct("foo", "foo's description");
  auto& bar = foo.addStruct("bar", "bar's description");
  bar.addInt("baz", "baz's description").range(0, 5);

  std::vector<VerificationError> errors;
  EXPECT_FALSE(inlet.verify(&errors));
  // Need something about "foo/bar/baz" and "range"
  EXPECT_TRUE(
    std::any_of(errors.begin(), errors.end(), [](const VerificationError& err) {
      return (err.path == axom::Path {"foo/bar/baz"}) &&
        err.messageContains("range");
    }));
}

TYPED_TEST(inlet_errors, invalid_value)
{
  std::string testString = "foo = 7";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  inlet.addInt("foo", "foo's description").validValues({2, 4, 6, 8});

  std::vector<VerificationError> errors;
  EXPECT_FALSE(inlet.verify(&errors));
  // Need something about "foo" and "valid value"
  EXPECT_TRUE(
    std::any_of(errors.begin(), errors.end(), [](const VerificationError& err) {
      return (err.path == axom::Path {"foo"}) &&
        err.messageContains("valid value");
    }));
}

TYPED_TEST(inlet_errors, invalid_value_nested)
{
  std::string testString = "foo = { bar = { baz = 7 } }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& foo = inlet.addStruct("foo", "foo's description");
  auto& bar = foo.addStruct("bar", "bar's description");
  bar.addInt("baz", "baz's description").validValues({2, 4, 6, 8});

  std::vector<VerificationError> errors;
  EXPECT_FALSE(inlet.verify(&errors));
  // Need something about "foo/bar/baz" and "valid value"
  EXPECT_TRUE(
    std::any_of(errors.begin(), errors.end(), [](const VerificationError& err) {
      return (err.path == axom::Path {"foo/bar/baz"}) &&
        err.messageContains("valid value");
    }));
}

TYPED_TEST(inlet_errors, field_verifier)
{
  std::string testString = "foo = 7";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  inlet.addInt("foo", "foo's description")
    .registerVerifier([](const axom::inlet::Field& field) {
      return (field.get<int>() % 2) == 0;
    });

  std::vector<VerificationError> errors;
  EXPECT_FALSE(inlet.verify(&errors));
  // Need something about "foo" and "lambda"
  EXPECT_TRUE(
    std::any_of(errors.begin(), errors.end(), [](const VerificationError& err) {
      return (err.path == axom::Path {"foo"}) && err.messageContains("lambda");
    }));
}

TYPED_TEST(inlet_errors, field_verifier_nested)
{
  std::string testString = "foo = { bar = { baz = 7 } }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& foo = inlet.addStruct("foo", "foo's description");
  auto& bar = foo.addStruct("bar", "bar's description");
  bar.addInt("baz", "baz's description")
    .registerVerifier([](const axom::inlet::Field& field) {
      return (field.get<int>() % 2) == 0;
    });

  std::vector<VerificationError> errors;
  EXPECT_FALSE(inlet.verify(&errors));
  // Need something about "foo/bar/baz" and "lambda"
  EXPECT_TRUE(
    std::any_of(errors.begin(), errors.end(), [](const VerificationError& err) {
      return (err.path == axom::Path {"foo/bar/baz"}) &&
        err.messageContains("lambda");
    }));
}

TYPED_TEST(inlet_errors, container_required)
{
  std::string testString = "foo = { }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& foo = inlet.addStruct("foo", "foo's description").required();
  foo.addInt("bar", "bar's description");

  std::vector<VerificationError> errors;
  EXPECT_FALSE(inlet.verify(&errors));
  // Need something about "foo" and "Required"
  EXPECT_TRUE(
    std::any_of(errors.begin(), errors.end(), [](const VerificationError& err) {
      return (err.path == axom::Path {"foo"}) &&
        err.messageContains("Required");
    }));
}

TYPED_TEST(inlet_errors, container_verifier)
{
  std::string testString = "foo = { bar = 7 }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& foo = inlet.addStruct("foo", "foo's description");
  foo.addInt("bar", "bar's description");

  // Bogus verifier, will fail because "bar" exists
  foo.registerVerifier([](const axom::inlet::Container& container) {
    return !container.exists();
  });

  std::vector<VerificationError> errors;
  EXPECT_FALSE(inlet.verify(&errors));
  // Need something about "foo" and "verification"
  EXPECT_TRUE(
    std::any_of(errors.begin(), errors.end(), [](const VerificationError& err) {
      return (err.path == axom::Path {"foo"}) &&
        err.messageContains("verification");
    }));
}

TYPED_TEST(inlet_errors, container_strict)
{
  std::string testString = "foo = { bar = 7, baz = 12 }";
  Inlet inlet = createBasicInlet<TypeParam>(testString);

  auto& foo = inlet.addStruct("foo", "foo's description").strict();
  foo.addInt("bar", "bar's description");

  std::vector<VerificationError> errors;
  EXPECT_FALSE(inlet.verify(&errors));
  // Need something about "foo" and "unexpected" and "baz"
  EXPECT_TRUE(
    std::any_of(errors.begin(), errors.end(), [](const VerificationError& err) {
      return (err.path == axom::Path {"foo"}) &&
        err.messageContains("unexpected") && err.messageContains("baz");
    }));
}

#ifdef AXOM_USE_SOL

TEST(inlet_errors_lua, function_verifier)
{
  std::string testString = "function foo (v) return 2*v end";

  Inlet inlet = createBasicInlet<axom::inlet::LuaReader>(testString);

  using axom::inlet::FunctionTag;
  using axom::inlet::FunctionType;
  auto& func = inlet.addFunction("foo",
                                 FunctionTag::Vector,
                                 {FunctionTag::Vector},
                                 "foo's description");

  func.registerVerifier([](const axom::inlet::Function& func) {
    auto result = func.call<FunctionType::Vector>(FunctionType::Vector {1, 1, 1});
    return result[0] == 0;
  });

  std::vector<VerificationError> errors;
  EXPECT_FALSE(inlet.verify(&errors));
  // Need something about "foo" and "verification"
  EXPECT_TRUE(
    std::any_of(errors.begin(), errors.end(), [](const VerificationError& err) {
      return (err.path == axom::Path {"foo"}) &&
        err.messageContains("verification");
    }));
}

#endif

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

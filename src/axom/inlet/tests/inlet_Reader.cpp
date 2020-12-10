// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include <string>
#include <vector>
#include <memory>

#include "axom/config.hpp"

#ifdef AXOM_USE_SOL
  #include "axom/inlet/LuaReader.hpp"
#endif

#include "axom/inlet/YAMLReader.hpp"

template <typename InletReader>
class InletReaderTest : public testing::Test
{ };

TYPED_TEST_SUITE_P(InletReaderTest);

namespace detail
{
class LuaToYAML
{
  static std::vector<std::string> tokenize(const std::string& text)
  {
    std::vector<std::string> result;
    std::size_t pos = 0;
    auto add_to_result = [&result](std::string&& str) {
      std::vector<std::string> parsed_after;
      std::size_t pos;
      while(str.length() > 1 &&
            ((pos = str.find_first_of("{}[];,")) != std::string::npos))
      {
        if(pos == 0)
        {
          result.emplace_back(1, str[pos]);
          str = str.substr(1, str.length() - 1);
        }
        else
        {
          parsed_after.emplace_back(1, str[pos]);
          str = str.substr(0, str.length() - 1);
        }
      }
      if(!str.empty())
      {
        result.push_back(str);
      }
      result.insert(result.end(), parsed_after.rbegin(), parsed_after.rend());
    };

    while(pos < text.length())
    {
      if(text[pos] == '"')
      {
        auto ending_quote = text.find('"', pos + 1);
        add_to_result(text.substr(pos, ending_quote - pos + 1));
        pos = ending_quote + 1;
      }
      else if(text[pos] == '[')
      {
        auto ending_bracket = text.find(']', pos + 1);
        add_to_result(text.substr(pos, ending_bracket - pos + 1));
        pos = ending_bracket + 1;
      }
      else
      {
        auto ending_space = text.find(' ', pos + 1);
        if(ending_space == std::string::npos)
        {
          if(text[pos] == ' ')
          {
            pos++;
          }
          add_to_result(text.substr(pos, text.length() - pos));
          pos = text.length();
        }
        else
        {
          add_to_result(text.substr(pos, ending_space - pos));
          pos = ending_space + 1;
        }
      }
    }
    return result;
  }

public:
  static std::string convert(const std::string& luaString)
  {
    const auto tokens = tokenize(luaString);
    std::size_t i = 0;
    std::string indent = "";
    std::string result;
    while(i < tokens.size())
    {
      const auto& token = tokens[i];
      if((i < tokens.size() - 2) && tokens[i + 1] == "=" && tokens[i + 2] == "{")
      {
        result += indent + token + ":\n";
        indent += "  ";
        i += 3;
      }
      else if(token == "}")
      {
        indent = indent.substr(2);
        i += 1;
      }
      else if(token == "," || token == ";")
      {
        i += 1;
      }
      else if(token == "[")
      {
        const auto& key = tokens[i + 1];

        const auto& value = tokens[i + 4];
        if(isdigit(key.front()))
        {
          // Assume this is an integer-keyed array
          result += indent + "- " + value + "\n";
        }
        else
        {
          result += indent + key + ": " + value + "\n";
        }
        i += 5;
      }
      else
      {
        const auto& value = tokens[i + 2];
        result += indent + token + ": " + value + "\n";
        i += 3;
      }
    }
    return result;
  }
};

template <typename InletReader>
std::string fromLuaTo(const std::string& luaString)
{
  return luaString;
}

template <>
std::string fromLuaTo<axom::inlet::YAMLReader>(const std::string& luaString)
{
  return LuaToYAML::convert(luaString);
}
}  // namespace detail

TYPED_TEST_P(InletReaderTest, getTopLevelBools)
{
  TypeParam reader;
  reader.parseString(detail::fromLuaTo<TypeParam>("foo = true; bar = false"));

  bool value, retValue;

  value = false;
  retValue = reader.getBool("foo", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, true);

  value = true;
  retValue = reader.getBool("bar", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, false);
}

TYPED_TEST_P(InletReaderTest, getInsideBools)
{
  TypeParam reader;
  reader.parseString(
    detail::fromLuaTo<TypeParam>("foo = { bar = false; baz = true }"));

  bool value, retValue;

  value = true;
  retValue = reader.getBool("foo/bar", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, false);

  value = false;
  retValue = reader.getBool("foo/baz", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, true);
}

TYPED_TEST_P(InletReaderTest, getTopLevelStrings)
{
  TypeParam reader;
  reader.parseString(detail::fromLuaTo<TypeParam>(
    "foo = \"this is a test string\"; bar = \"TesT StrInG\""));

  bool retValue;
  std::string value;

  value = "";
  retValue = reader.getString("foo", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, "this is a test string");

  value = "";
  retValue = reader.getString("bar", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, "TesT StrInG");
}

TYPED_TEST_P(InletReaderTest, getInsideStrings)
{
  TypeParam reader;
  reader.parseString(detail::fromLuaTo<TypeParam>(
    "foo = { bar = \"this is a test string\"; baz = \"TesT StrInG\" }"));

  bool retValue;
  std::string value;

  value = "";
  retValue = reader.getString("foo/bar", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, "this is a test string");

  value = "";
  retValue = reader.getString("foo/baz", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, "TesT StrInG");
}

TYPED_TEST_P(InletReaderTest, mixLevelTables)
{
  TypeParam reader;
  reader.parseString(detail::fromLuaTo<TypeParam>(
    "t = { innerT = { foo = 1 }, anotherInnerT = {baz = 3}}"));

  bool retValue;
  int value;

  value = 0;
  retValue = reader.getInt("t/innerT/foo", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, 1);

  value = 0;
  retValue = reader.getInt("t/doesntexist", value);
  EXPECT_EQ(retValue, false);
  EXPECT_EQ(value, 0);

  value = 0;
  retValue = reader.getInt("t/anotherInnerT/baz", value);
  EXPECT_EQ(retValue, true);
  EXPECT_EQ(value, 3);
}

TYPED_TEST_P(InletReaderTest, getMap)
{
  // Keep this contiguous in order to test all supported input languages
  std::string testString =
    "luaArray = { [0] = 4, [1] = 5, [2] = 6 , [3] = true, [4] = false, [5] = "
    "2.4, [6] = 'hello', [7] = 'bye' }";
  TypeParam reader;
  reader.parseString(detail::fromLuaTo<TypeParam>(testString));

  std::unordered_map<int, int> ints;
  bool found = reader.getIntMap("luaArray", ints);
  EXPECT_TRUE(found);
  std::unordered_map<int, int> expectedInts {{0, 4}, {1, 5}, {2, 6}, {5, 2}};
  EXPECT_EQ(expectedInts, ints);

  std::unordered_map<int, double> doubles;
  found = reader.getDoubleMap("luaArray", doubles);
  EXPECT_TRUE(found);
  std::unordered_map<int, double> expectedDoubles {{0, 4},
                                                   {1, 5},
                                                   {2, 6},
                                                   {5, 2.4}};
  EXPECT_EQ(expectedDoubles, doubles);

  std::unordered_map<int, bool> bools;
  found = reader.getBoolMap("luaArray", bools);
  EXPECT_TRUE(found);
  std::unordered_map<int, bool> expectedBools {{3, true}, {4, false}};
  EXPECT_EQ(expectedBools, bools);

  // Conduit's YAML parser doesn't distinguish boolean literals from strings
  // so the YAML version will extract the "true" and "false" here
  std::unordered_map<int, std::string> strs;
  found = reader.getStringMap("luaArray", strs);
  EXPECT_TRUE(found);
  // std::unordered_map<int, std::string> expectedStrs {{6, "hello"},
  //                                                    {7, "bye"}};
  // EXPECT_EQ(expectedStrs, strs);
}

// Checks that LuaReader parses array information as expected
// Discontiguous arrays are lua-specific
TEST(LuaReaderTest, getDiscontiguousMap)
{
  std::string testString =
    "luaArray = { [1] = 4, [2] = 5, [3] = 6 , [4] = true, [8] = false, [12] = "
    "2.4, [33] = 'hello', [200] = 'bye' }";
  axom::inlet::LuaReader reader;
  reader.parseString(testString);

  std::unordered_map<int, int> ints;
  bool found = reader.getIntMap("luaArray", ints);
  EXPECT_TRUE(found);
  std::unordered_map<int, int> expectedInts {{1, 4}, {2, 5}, {3, 6}, {12, 2}};
  EXPECT_EQ(expectedInts, ints);

  std::unordered_map<int, double> doubles;
  found = reader.getDoubleMap("luaArray", doubles);
  EXPECT_TRUE(found);
  std::unordered_map<int, double> expectedDoubles {{1, 4},
                                                   {2, 5},
                                                   {3, 6},
                                                   {12, 2.4}};
  EXPECT_EQ(expectedDoubles, doubles);

  std::unordered_map<int, bool> bools;
  found = reader.getBoolMap("luaArray", bools);
  EXPECT_TRUE(found);
  std::unordered_map<int, bool> expectedBools {{4, true}, {8, false}};
  EXPECT_EQ(expectedBools, bools);

  std::unordered_map<int, std::string> strs;
  found = reader.getStringMap("luaArray", strs);
  EXPECT_TRUE(found);
  std::unordered_map<int, std::string> expectedStrs {{33, "hello"},
                                                     {200, "bye"}};
  EXPECT_EQ(expectedStrs, strs);
}

REGISTER_TYPED_TEST_SUITE_P(InletReaderTest,
                            getTopLevelBools,
                            getInsideBools,
                            getTopLevelStrings,
                            getInsideStrings,
                            mixLevelTables,
                            getMap);

#ifdef AXOM_USE_SOL
using ReaderTypes =
  ::testing::Types<axom::inlet::LuaReader, axom::inlet::YAMLReader>;
#else
using ReaderTypes = ::testing::Types<axom::inlet::YAMLReader>;
#endif

INSTANTIATE_TYPED_TEST_SUITE_P(ReaderTests, InletReaderTest, ReaderTypes);

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

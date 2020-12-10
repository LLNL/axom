// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/inlet/tests/inlet_test_utils.hpp"

namespace axom::inlet::detail
{
std::vector<std::string> LuaToYAML::tokenize(const std::string& text)
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
        // The position could be partway through the string, but
        // we need the last entry
        pos = str.find_last_of("{}[];,");
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
      // Remove the leading space if there is one
      if(text[pos] == ' ')
      {
        pos++;
      }
      if(ending_space == std::string::npos)
      {
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

std::string LuaToYAML::convert(const std::string& luaString)
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
      // Check if this is the start of a new table
      if(tokens[i + 4] == "{")
      {
        result += indent + key + ":\n";
        indent += "  ";
      }
      // Or if it's a terminal entry in the current table
      else
      {
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

}  // namespace axom::inlet::detail

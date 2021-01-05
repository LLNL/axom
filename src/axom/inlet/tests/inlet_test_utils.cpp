// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/inlet/tests/inlet_test_utils.hpp"
#include "axom/slic.hpp"

namespace axom
{
namespace inlet
{
namespace detail
{
void LuaToYAML::add_token(std::string&& token, std::vector<std::string>& tokens)
{
  const static std::string punctuation = "{}[];,=";
  std::vector<std::string> parsed_after;
  std::size_t pos;
  while(token.length() > 1 &&
        ((pos = token.find_first_of(punctuation)) != std::string::npos))
  {
    if(pos == 0)
    {
      tokens.emplace_back(1, token[pos]);
      token = token.substr(1, token.length() - 1);
    }
    else
    {
      // The position could be partway through the string, but
      // we need the last entry
      pos = token.find_last_of(punctuation);
      parsed_after.emplace_back(1, token[pos]);
      token = token.substr(0, token.length() - 1);
    }
  }
  if(!token.empty())
  {
    tokens.push_back(token);
  }
  tokens.insert(tokens.end(), parsed_after.rbegin(), parsed_after.rend());
}

std::vector<std::string> LuaToYAML::tokenize(const std::string& text)
{
  std::vector<std::string> result;
  std::size_t pos = 0;

  // Adds a token starting at the current "pos" up through the next instance
  // of the selected "end_char"
  auto find_and_add_until = [&text, &result, &pos](const char end_char) {
    auto ending_pos = text.find(end_char, pos + 1);
    add_token(text.substr(pos, ending_pos - pos + 1), result);
    pos = ending_pos + 1;
  };

  while(pos < text.length())
  {
    // For quoted strings, find the next quote and use that as the end
    if(text[pos] == '"')
    {
      find_and_add_until('"');
    }
    else if(text[pos] == '\'')
    {
      find_and_add_until('\'');
    }
    // If it's a table index, find the next bracket and use that as the end
    else if(text[pos] == '[')
    {
      find_and_add_until(']');
    }
    // Otherwise just find the next instance of whitespace
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
        add_token(text.substr(pos, text.length() - pos), result);
        pos = text.length();
      }
      else
      {
        add_token(text.substr(pos, ending_space - pos), result);
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
    // The start of a new table
    if((i < tokens.size() - 2) && tokens[i + 1] == "=" && tokens[i + 2] == "{")
    {
      result += indent + token + ":\n";
      indent += "  ";
      i += 3;
    }
    // End of a table - only need to reduce the indent
    else if(token == "}")
    {
      indent = indent.substr(2);
      i += 1;
    }
    // Commas and semicolons can be ignored completely
    else if(token == "," || token == ";")
    {
      i += 1;
    }
    // The start of an implicitly indexed array
    else if(token == "{")
    {
      result += indent + "-\n";
      indent += "  ";
      i += 1;
    }
    else if(token == "[")
    {
      const auto& key = tokens[i + 1];
      // Check if this is the start of a new table
      if(tokens[i + 4] == "{")
      {
        result += indent;
        // Check if it's an integer-keyed table
        result += isdigit(key.front()) ? "-\n" : key + ":\n";
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

}  // namespace detail

}  // namespace inlet

}  // namespace axom

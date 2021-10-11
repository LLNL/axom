// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/utilities/StringUtilities.hpp"

#include <algorithm>
#include <cctype>
#include <iostream>

#include "axom/fmt.hpp"

namespace axom
{
namespace utilities
{
namespace string
{
std::vector<std::string> split(const std::string& str, const char delimiter)
{
  std::vector<std::string> tokens;
  std::string token;
  std::istringstream tokenStream(str);
  while(std::getline(tokenStream, token, delimiter))
  {
    tokens.push_back(token);
  }
  return tokens;
}

void toLower(std::string& str)
{
  std::transform(str.begin(), str.end(), str.begin(), [](const unsigned char c) {
    return std::tolower(c);
  });
}

void toUpper(std::string& str)
{
  std::transform(str.begin(), str.end(), str.begin(), [](const unsigned char c) {
    return std::toupper(c);
  });
}

std::vector<std::string> splitLastNTokens(const std::string& input,
                                          const std::size_t n,
                                          const char delim)
{
  std::vector<std::string> result;

  auto last_pos = std::string::npos;
  auto pos = input.find_last_of(delim, last_pos - 1);

  if(n > 0 && !input.empty())
  {
    while((pos != std::string::npos) && (result.size() < n - 1))
    {
      result.push_back(input.substr(pos + 1, last_pos - pos - 1));
      last_pos = pos;
      pos = input.find_last_of(delim, last_pos - 1);
    }

    // Add the rest of the string (first token)
    result.push_back(input.substr(0, last_pos));
    std::reverse(result.begin(), result.end());
  }

  return result;
}

std::string appendPrefix(const std::string& prefix, const std::string& name)
{
  return (prefix == "") ? name : prefix + "/" + name;
}

std::string removePrefix(const std::string& prefix, const std::string& name)
{
  if(prefix.empty())
  {
    return name;
  }
  else if(axom::utilities::string::startsWith(name, prefix + "/"))
  {
    return name.substr(prefix.size());
  }
  std::cerr << axom::fmt::format(
                 "Provided name {0} does not "
                 "contain prefix {1}",
                 name,
                 prefix)
            << std::endl;
  return name;
}

std::string removeBeforeDelimiter(const std::string& str, const char delim)
{
  auto pos = str.find_last_of(delim);
  // Will return an empty string if the delimiter was not found
  return str.substr(pos + 1);
}

std::string removeAllInstances(const std::string& target,
                               const std::string& substr)
{
  std::string result = target;
  auto pos = result.find(substr);
  while(pos != std::string::npos)
  {
    result.erase(pos, substr.length());
    pos = result.find(substr);
  }
  return result;
}

}  // end namespace string
}  // end namespace utilities
}  // end namespace axom

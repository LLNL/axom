// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/utilities/StringUtilities.hpp"

#include <algorithm>

namespace axom
{
namespace utilities
{
namespace string
{
void split(std::vector<std::string>& tokens,
           const std::string& str,
           const char delimiter)
{
  std::string token;
  std::istringstream tokenStream(str);
  while(std::getline(tokenStream, token, delimiter))
  {
    tokens.push_back(token);
  }
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

  while((pos != std::string::npos) && (result.size() < n - 1))
  {
    result.push_back(input.substr(pos + 1, last_pos - pos - 1));
    last_pos = pos;
    pos = input.find_last_of(delim, last_pos - 1);
  }
  // Add the rest of the string (first token)
  result.push_back(input.substr(0, last_pos));
  std::reverse(result.begin(), result.end());
  return result;
}

}  // end namespace string
}  // end namespace utilities
}  // end namespace axom

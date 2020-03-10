// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/utilities/StringUtilities.hpp"

namespace axom
{
namespace utilities
{
namespace string
{

void split(std::vector<std::string>& tokens, const std::string& str, const char delimiter)
{
  std::string token;
  std::istringstream tokenStream(str);
  while (std::getline(tokenStream, token, delimiter)) {
    tokens.push_back(token);
  }
}

} // end namespace string
} // end namespace utilities
} // end namespace axom

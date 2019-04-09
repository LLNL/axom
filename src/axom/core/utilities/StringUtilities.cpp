// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/utilities/StringUtilities.hpp"

#include "axom/config.hpp"

#ifndef AXOM_USE_CXX11
#include <sstream>
#endif

namespace axom
{
namespace utilities
{
namespace string
{

std::string intToString(int intValue)
{
  std::string stringValue = "";
#ifdef AXOM_USE_CXX11
  stringValue += std::to_string(intValue);
#else
  std::ostringstream ss;
  ss << intValue;
  stringValue += ss.str();
#endif
  return stringValue;
}

int stringToInt(const std::string& stringValue)
{
  int intValue = 0;
#ifdef AXOM_USE_CXX11
  intValue = stoi(stringValue);
#else
  std::istringstream(stringValue) >> intValue;
#endif
  return intValue;
}

} // end namespace string
} // end namespace utilities
} // end namespace axom

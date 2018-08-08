/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

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

/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

#include "common/StringUtilities.hpp"

#include "common/config.hpp"

#ifndef ATK_USE_CXX11
#include <sstream>
#endif

namespace asctoolkit {
namespace utilities {
namespace string {

std::string intToString(int intValue)
{
    std::string stringValue = "";
#ifdef ATK_USE_CXX11
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
#ifdef ATK_USE_CXX11
    intValue = stoi(stringValue);
#else
    std::istringstream(stringValue) >> intValue;
#endif
    return intValue;
}

} // end namespace string
} // end namespace utilities
} // end namespace asctoolkit

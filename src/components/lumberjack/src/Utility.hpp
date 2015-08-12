/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

#ifndef UTILITY_HPP
#define UTILITY_HPP

#include <string>

namespace asctoolkit {
namespace lumberjack {

std::string intToString(int intValue);
int stringToInt(const std::string& stringValue);

}
}

#endif
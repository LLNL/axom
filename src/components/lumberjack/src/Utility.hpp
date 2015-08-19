/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

/*!
 *******************************************************************************
 * \file Utility.hpp
 * \author Chris White (white238@llnl.gov)
 *
 * \brief This file contains the definitions of utility functions.
 *******************************************************************************
 */

#ifndef UTILITY_HPP
#define UTILITY_HPP

#include <string>

namespace asctoolkit {
namespace lumberjack {

/*!
*****************************************************************************
* \brief Returns the converted integer as a string.
*
* \param [in] intValue Integer to be converted to a string.
*****************************************************************************
*/
std::string intToString(int intValue);

/*!
*****************************************************************************
* \brief Returns the converted string as a integer.
*
* \param [in] stringValue String to be converted to a integer.
*****************************************************************************
*/
int stringToInt(const std::string& stringValue);

} // end namespace lumberjack
} // end namespace asctoolkit

#endif

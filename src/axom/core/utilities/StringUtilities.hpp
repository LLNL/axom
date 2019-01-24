/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
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

#ifndef COMMON_STRING_UTILITIES_H_
#define COMMON_STRING_UTILITIES_H_

#include <string>

namespace axom
{
namespace utilities
{
namespace string
{

/*!
 * \brief Returns the converted integer as a string.
 *
 * \param [in] intValue Integer to be converted to a string.
 */
std::string intToString(int intValue);

/*!
 * \brief Returns the converted string as a integer.
 *
 * \param [in] stringValue String to be converted to a integer.
 */
int stringToInt(const std::string& stringValue);


} // end namespace string
} // end namespace utilities
} // end namespace axom

#endif //  COMMON_STRING_UTILITIES_H_

// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef CORE_SYSTEM_UTILITIES_H_
#define CORE_SYSTEM_UTILITIES_H_

#include <string>

namespace axom
{
namespace utilities
{
/**
 * @brief Returns the name of the machine
 *
 * @return The name of the current machine, empty string on failure
 */
std::string getHostName();

/**
 * @brief Returns the name of the current user
 *
 * @return The name of the current user, empty string on failure
 */
std::string getUserName();

}  // end namespace utilities
}  // end namespace axom

#endif  //  CORE_SYSTEM_UTILITIES_H_

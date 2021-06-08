// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_UTILS_ABOUT_H_
#define AXOM_UTILS_ABOUT_H_

#include <ostream>  // for std::ostream
#include <string>   // for std::string

namespace axom
{
/*!
 * \brief Prints info about how Axom was configured and built to stdout
 *
 */
void about();

/*!
 * \brief Prints info about how Axom was configured and built to a stream
 *
 * \param [in,out] oss the target stream where to append the Axom info
 */
void about(std::ostream &oss);

/*!
 * \brief Returns a string consisting of the Axom version.
 *
 * \return str string corresponding to the Axom version
 * \post str != ""
 */
std::string getVersion();

}  // end namespace axom

#endif  //  AXOM_UTILS_ABOUT_H_

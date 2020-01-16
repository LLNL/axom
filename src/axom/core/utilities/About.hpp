// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_UTILS_ABOUT_H_
#define AXOM_UTILS_ABOUT_H_

#include <ostream>

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
 */
void about(std::ostream &oss);

} // end namespace axom

#endif //  AXOM_UTILS_ABOUT_H_

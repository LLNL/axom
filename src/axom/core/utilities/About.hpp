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

// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file Annotations.hpp
 *
 * \brief Defines functions to help with performance annotations
 */

#ifndef AXOM_CORE_ANNOTATIONS_HPP_
#define AXOM_CORE_ANNOTATIONS_HPP_

#include "axom/config.hpp"

#ifdef AXOM_USE_ADIAK
  #include "adiak.hpp"
#endif

#ifdef AXOM_USE_CALIPER
  #include "caliper/cali-manager.h"
  #include "caliper/cali.h"
#endif

namespace axom
{
namespace utilities
{
namespace annotations
{
void initialize();

void finalize();

}  // namespace annotations
}  // namespace utilities
}  // namespace axom

#endif  // AXOM_CORE_ANNOTATIONS_HPP_

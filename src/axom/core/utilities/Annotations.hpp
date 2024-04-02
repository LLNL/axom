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
namespace detail
{
void initialize_adiak();
void initialize_caliper(const std::string& mode, int num_ranks);
}  // namespace detail

void initialize(const std::string& mode, int num_ranks);

void finalize();

/// Declares metadata for this run
template <typename T>
void declare_metadata(const std::string& name,
                      const T& value,
                      std::string category = "")
{
#ifdef AXOM_USE_ADIAK
  detail::initialize_adiak();
  adiak::value(name, value, adiak_general, category);
#endif
}

}  // namespace annotations
}  // namespace utilities
}  // namespace axom

#endif  // AXOM_CORE_ANNOTATIONS_HPP_

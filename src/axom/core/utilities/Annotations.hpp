// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file Annotations.hpp
 *
 * \brief Defines functions to help with performance annotations
 * 
 * The annotations API and macros are always available but they are effectively no-ops
 * unless axom is built with caliper and adiak support
 */

#ifndef AXOM_CORE_ANNOTATIONS_HPP_
#define AXOM_CORE_ANNOTATIONS_HPP_

#include "axom/config.hpp"
#include "axom/core/Macros.hpp"

#ifdef AXOM_USE_ADIAK
  #include "adiak.hpp"
#endif

#ifdef AXOM_USE_CALIPER
  #include "caliper/cali-manager.h"
  #include "caliper/cali.h"
#endif

#include <string>

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

bool check_mode(const std::string& mode);
std::string help_string();
}  // namespace detail

/// Initializes the Annotation API
void initialize(const std::string& mode, int num_ranks);

/// Finalizes and flushes the Annotation API
void finalize();

/// Begins an annotation within a region
void begin(const std::string& name);

/// Ends an annotation within a region
void end(const std::string& name);

/// Declares metadata for this run
template <typename T>
void declare_metadata(const std::string& name,
                      const T& value,
                      std::string category = "")
{
#ifdef AXOM_USE_ADIAK
  detail::initialize_adiak();
  adiak::value(name, value, adiak_general, category);
#else
  AXOM_UNUSED_VAR(name);
  AXOM_UNUSED_VAR(value);
  AXOM_UNUSED_VAR(category);
#endif
}

}  // namespace annotations
}  // namespace utilities
}  // namespace axom

#endif  // AXOM_CORE_ANNOTATIONS_HPP_

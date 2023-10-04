// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *
 * \file NumericLimits.hpp
 *
 * \brief Header file containing portability layer for std::numeric_limits
 *        capabilities
 *
 */

#ifndef AXOM_NUMERICLIMITS_HPP_
#define AXOM_NUMERICLIMITS_HPP_

#include "axom/config.hpp"  // for compile-time definitions

#include <limits>

#if defined(AXOM_USE_CUDA)
#include <cuda/std/limits>
#endif

namespace axom
{

#if defined(AXOM_USE_CUDA) && defined(AXOM_DEVICE_CODE)
  template<typename T>
  using numeric_limits = cuda::std::numeric_limits<T>;
#else
  template<typename T>
  using numeric_limits = std::numeric_limits<T>;
#endif

}  // namespace axom

#endif  // AXOM_NUMERICLIMITS_HPP_

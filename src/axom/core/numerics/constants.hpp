// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_NUMERICS_CONSTANTS_HPP_
#define AXOM_NUMERICS_CONSTANTS_HPP_

#include "axom/core/Macros.hpp"  // for Axom macros

namespace axom
{
namespace numerics
{
namespace constants
{
/*!
 * \brief Return the value of e.
 *
 * \return The value of e.
 */
constexpr inline AXOM_HOST_DEVICE double e() { return 2.7182818284590452354; }

/*!
 * \brief Return the value of pi.
 *
 * \return The value of pi.
 *
 * \note This is a constexpr function here since Windows does not define M_PI
 *       without an additional defines. Let's just define it here.
 */
constexpr inline AXOM_HOST_DEVICE double pi() { return 3.14159265358979323846; }

/*!
 * \brief Return the value of sqrt(2)
 *
 * \return The value of sqrt(2)
 */
constexpr inline AXOM_HOST_DEVICE double sqrt2() { return 1.41421356237309504880; }

}  // namespace constants
}  // namespace numerics
}  // namespace axom

#endif

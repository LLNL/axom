// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file in_polygon.hpp
 *
 * \brief Consists of methods that test whether a given query point is
 * inside a given polygon.
 *
 * Uses a ray casting algorithm
 */

#ifndef AXOM_PRIMAL_IN_POLYGON_HPP_
#define AXOM_PRIMAL_IN_POLYGON_HPP_

// Axom includes
#include "axom/config.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Polygon.hpp"
#include "axom/primal/operators/winding_number.hpp"

// C++ includes
#include <cmath>

namespace axom
{
namespace primal
{
/*!
 * \brief Determines containment for a point in a polygon
 *
 * \param [in] query The query point to test
 * \param [in] poly The Polygon object to test for containment
 * \param [in] includeBoundary If true, points on the boundary are considered interior.
 * \param [in] useNonzeroRule If false, use even/odd protocol for inclusion
 * \param [in] EPS The tolerance level for collinearity
 * 
 * Determines containment using the winding number with respect to the 
 * given polygon. 
 * Different protocols determine containment from the winding number differently.
 *   Nonzero Rule:  If the winding number is nonzero, the point is interior.
 *   Even/Odd rule: If the winding number is odd, it is interior. Exterior otherwise.
 *
 * \return boolean value indicating containment.
 */
template <typename T>
bool in_polygon(const Point<T, 2>& query,
                const Polygon<T, 2>& poly,
                bool includeBoundary = false,
                bool useNonzeroRule = true,
                double EPS = 1e-8)
{
  return useNonzeroRule
    ? winding_number(query, poly, includeBoundary, EPS) != 0
    : (winding_number(query, poly, includeBoundary, EPS) % 2) == 1;
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_IN_CURVED_POLYGON_H_

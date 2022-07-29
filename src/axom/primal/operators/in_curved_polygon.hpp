// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file in_curved_polygon.hpp
 *
 * \brief Consists of methods that test whether a given query point is
 * inside a given curved polygon.
 *
 * Uses an adaptive winding number calculation
 */

#ifndef AXOM_PRIMAL_IN_CURVED_POLYGON_HPP_
#define AXOM_PRIMAL_IN_CURVED_POLYGON_HPP_

// Axom includes
#include "axom/config.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/geometry/CurvedPolygon.hpp"
#include "axom/primal/operators/detail/in_curved_polygon_impl.hpp"

namespace axom
{
namespace primal
{
/*!
 * \brief Robustly determine if query point is interior to a curved polygon
 *
 * \param [in] query The query point to test
 * \param [in] cpoly The CurvedPolygon object to test for containment
 * \param [in] linear_tol The tolerance level at which a Bezier curve is linear
 * \param [in] edge_tol The tolerance level at which the query point is on the curve
 * 
 * Determines contianment using the (rounded) winding number with respect
 * to the given curved polygon. This algorithm is robust, as the winding number is rounded 
 * in the final in/out determination. 
 * Different protocols determine containment from the winding number differently.
 *   Nonzero Rule:  If the winding number is nonzero, the point is interior.
 *   Even/Odd rule: If the winding number is odd, it is interior. Exterior otherwise.
 *
 * \return A boolean value indicating containment.
 */
template <typename T>
inline bool in_curved_polygon(const Point<T, 2>& query,
                              const CurvedPolygon<T, 2>& cpoly,
                              const bool useNonzeroRule = true,
                              const double linear_tol = 1e-8,
                              const double edge_tol = 1e-8)
{
  double winding_num = winding_number(query, cpoly, linear_tol, edge_tol);

  // Else, use EvenOdd rule
  return useNonzeroRule ? std::round(winding_num)
                        : (std::lround(winding_num) % 2) == 1;
}

/*!
 * \brief Computes the generalized winding number for a curved polygon
 *
 * \param [in] query The query point to test
 * \param [in] cpoly The CurvedPolygon object 
 * \param [in] linear_tol The tolerance level at which a Bezier curve is linear
 * \param [in] edge_tol The tolerance level at which the query point is on the curve
 *
 * Computes the winding number using a recursive, bisection algorithm.
 * Iterates over the edges of a curved polygon object, and uses nearly-linear 
 * Bezier curves as a base case.
 * 
 * \return float the generalized winding number.
 */
template <typename T>
double winding_number(const Point<T, 2>& q,
                      const CurvedPolygon<T, 2>& cpoly,
                      const double linear_tol = 1e-8,
                      const double edge_tol = 1e-8)
{
  double ret_val = 0.0;
  for(int i = 0; i < cpoly.numEdges(); i++)
    ret_val +=
      detail::adaptive_winding_number(q, cpoly[i], false, linear_tol, edge_tol);

  return ret_val;
}

/*!
 * \brief Computes the generalized winding number for a single Bezier curve
 *
 * \param [in] query The query point to test
 * \param [in] cpoly The Bezier curve object 
 * \param [in] linear_tol The tolerance level at which a Bezier curve is linear
 * \param [in] edge_tol The tolerance level at which the query point is on the curve
 *
 * Computes the winding number using a recursive, bisection algorithm,
 * using nearly-linear Bezier curves as a base case.
 * 
 * \return float the generalized winding number.
 */
template <typename T>
double winding_number(const Point<T, 2>& q,
                      const BezierCurve<T, 2>& c,
                      const double linear_tol = 1e-8,
                      const double edge_tol = 1e-8)
{
  return detail::adaptive_winding_number(q, c, false, linear_tol, edge_tol);
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_IN_CURVED_POLYGON_H_

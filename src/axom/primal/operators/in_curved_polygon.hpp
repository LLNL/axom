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
 * \brief Robustly tests whether a query point lies inside curved polygon.
 * Basically boolean interface for winding_number
 */
template <typename T>
inline bool in_curved_polygon(const Point<T, 2>& query,
                              const CurvedPolygon<T, 2>& cpoly,
                              const double EPS = 1e-8)
{
  double ret_val = winding_number(query, cpoly, EPS);

  return !(std::round(ret_val) == 0);
}

// Base winding number function.
template <typename T>
double winding_number(const Point<T, 2>& q,
                      const CurvedPolygon<T, 2>& cpoly,
                      const double EPS = 1e-8)
{
  double ret_val = 0.0;
  for(int i = 0; i < cpoly.numEdges(); i++)
  {
    double this_val =
      detail::adaptive_winding_number(q, cpoly[i], total_depth, EPS);
    ret_val += this_val;
  }

  return ret_val;
}

// Overload for single bezier curve.
//  Assumes c has a convex bounding box
template <typename T>
double winding_number(const Point<T, 2>& q,
                      const BezierCurve<T, 2>& c,
                      const double EPS = 1e-8)
{
  return detail::adaptive_winding_number(q, c, EPS);
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_IN_CURVED_POLYGON_H_

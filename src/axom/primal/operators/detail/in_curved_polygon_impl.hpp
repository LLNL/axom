// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef PRIMAL_IN_CURVED_POLYGON_IMPL_HPP_
#define PRIMAL_IN_CURVED_POLYGON_IMPL_HPP_

// Axom includes
#include "axom/config.hpp"  // for compile-time configuration options
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/operators/in_polygon.hpp"

// C++ includes
#include <cmath>

namespace axom
{
namespace primal
{
namespace detail
{

// Close the bezier curve linearly, compute winding number of that closure.
template <typename T>
double closure_winding_number(const Point<T, 2>& q, const BezierCurve<T, 2>& c)
{
  Triangle<double, 2> T(c[c.getOrder()], q, c[0]);
  return 0.5 * M_1_PI * T.angle(1) * ((T.signedArea() < 0) ? 1.0 : -1.0);
}

// Recursive winding number function. Assumes convex c, and bisects c
//  until q is exterior to the bounding box
template <typename T>
double adaptive_winding_number(const Point2D& q,
                               const BezierCurve<T, 2>& c,
                               double linear_tol)
{
  double cl_winding_num = closure_winding_number(q, c);

  // Base case for recursion
  if(c.isLinear(linear_tol)) return -cl_winding_num;
  
  Polygon<T, 2> controlPolygon( c.getControlPoints() );

  // Means point is far enough away to ensure we are outside the closure
  if(!in_polygon(q, controlPolygon, linear_tol)) return -cl_winding_num;
  // Otherwise, our quadrature didn't give us a good enough answer, so we try again
  else
  {
    BezierCurve<T, 2> c1, c2;
    c.split(0.5, c1, c2);

    return adaptive_winding_number(q, c1, linear_tol) +
           adaptive_winding_number(q, c2, linear_tol);
  }
}

}  // end namespace detail
}  // end namespace primal
}  // end namespace axom

#endif
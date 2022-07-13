// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef PRIMAL_IN_CURVED_POLYGON_IMPL_HPP_
#define PRIMAL_IN_CURVED_POLYGON_IMPL_HPP_

// Axom includes
#include "axom/config.hpp"  // for compile-time configuration options
#include "axom/primal.hpp"

// MFEM includes
#ifdef AXOM_USE_MFEM
  #include "mfem.hpp"
#else
  #error "Primal's in/out functions for CurvedPolygon require mfem library."
#endif

// C++ includes
#include <cmath>

namespace axom
{
namespace primal
{
namespace detail
{

// Close the bezier curve, compute winding number of that closure.
template <typename T>
double closure_winding_number(BezierCurve<T, 2> c, Point2D q)
{
  Triangle<double, 2> T(q, c[0], c[c.getOrder()]);
  return 0.5 * M_1_PI * T.angle(0) * ((T.signedArea() < 0) ? -1.0 : 1.0);
}

// Get the function to be passed into the evaluate integral function
inline std::function<Vector2D(Point2D)> get_winding_func(Point2D p)
{
  return [p](Point2D x) -> Vector2D {
    double denom =
      2 * M_PI * ((x[0] - p[0]) * (x[0] - p[0]) + (x[1] - p[1]) * (x[1] - p[1]));
    return Vector2D({-(x[1] - p[1]) / denom, (x[0] - p[0]) / denom});
  };
}

// Recursive winding number function
template <typename T>
double adaptive_winding_number(const BezierCurve<T, 2>& c,
                               const Point2D& q,
                               int qnodes,
                               double itol,
                               double ltol,
                               int depth = 0)
{
  auto winding_func = get_winding_func(q);
  double winding_num = evaluate_line_integral(c, winding_func, qnodes);
  double cl_winding_num = closure_winding_number(c, q);
  double closed_winding = winding_num + cl_winding_num;

  // If the value is close to an integer, we know our quadrature is accurate.
  //  Return the value to machine precision.
  if(std::abs(std::round(closed_winding) - closed_winding) < itol)
  {
    return std::round(closed_winding) - cl_winding_num;
  }
  // Or, if the curve is linear, then it has no loops and we only
  //  need to know which side of the curve we are on.
  //  Serves as the base case for our recursion
  else if(c.isLinear(ltol))
  {
    return std::round(closed_winding) - cl_winding_num;
  }
  // Otherwise, our quadrature didn't give us a good enough answer, so we try again
  else
  {
    BezierCurve<double, 2> c1, c2;
    c.split(0.5, c1, c2);

    double wn1 = adaptive_winding_number(c1, q, qnodes, itol, ltol, depth + 1);
    double wn2 = adaptive_winding_number(c2, q, qnodes, itol, ltol, depth + 1);
    return wn1 + wn2;
  }
}

}  // end namespace detail
}  // end namespace primal
}  // end namespace axom

#endif
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
    //std::cout << x << std::endl;
    //std::cout << denom << std::endl;
    return Vector2D({-(x[1] - p[1]) / denom, (x[0] - p[0]) / denom});
  };
}

// Recursive winding number function
template <typename T>
double adaptive_winding_number(const BezierCurve<T, 2>& c,
                               const Point2D& q,
                               int qnodes,
                               int& depth,
                               double atol)
{
  auto winding_func = get_winding_func(q);
  double winding_num = evaluate_line_integral(c, winding_func, qnodes);
  double cl_winding_num = closure_winding_number(c, q);
  double closed_winding = winding_num + cl_winding_num;

  // If the value is close to an integer, we solved it
  if(std::abs(std::round(closed_winding) - closed_winding) < atol)
  {
    //std::cout << depth << std::endl;
    return winding_num;
  }
  // if the value is closed to a half-integer, we are likely on the curve.
  //  In this case, return a value as if we are just barely on the interior
  else if(std::abs(std::round(2 * closed_winding) - 2 * closed_winding) < atol)
  {
    //std::cout << winding_num << ", " << cl_winding_num << std::endl;
    return std::round(winding_num) + 0.5;
  }
  // Otherwise, our quadrature didn't give us a good enough answer, so we try again
  else
  {
    BezierCurve<double, 2> c1, c2;
    c.split(0.5, c1, c2);

    int depth1 = depth + 1;
    double wn1 = adaptive_winding_number(c1, q, qnodes, depth1, atol);

    int depth2 = depth + 1;
    double wn2 = adaptive_winding_number(c2, q, qnodes, depth2, atol);

    depth = std::max(depth1, depth2);
    return wn1 + wn2;
  }
}

}  // end namespace detail
}  // end namespace primal
}  // end namespace axom

#endif
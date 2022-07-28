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
#include "axom/primal/operators/is_convex.hpp"
#include "axom/primal/operators/squared_distance.hpp"

// C++ includes
#include <cmath>

namespace axom
{
namespace primal
{
namespace detail
{
/*
 * \brief Compute the "closure winding number" for a Bezier curve
 *
 * \param [in] query The query point to test
 * \param [in] c The BezierCurve object to close
 *
 * A possible "closure" of a Bezier curve is a straight line segment 
 * connecting its two endpoints. The "closure winding number" is the 
 * winding number of the query point with respect to this segment, which
 * has a direct formula.
 * 
 * \return 
 */
template <typename T>
double closure_winding_number(const Point<T, 2>& q, const BezierCurve<T, 2>& c)
{
  const int ord = c.getOrder();
  Vector<T, 2> V1 = Vector<T, 2>(q, c[0]).unitVector();
  Vector<T, 2> V2 = Vector<T, 2>(q, c[ord]).unitVector();

  // clang-format off
  double orient = axom::numerics::determinant(V1[0] - V2[0], V2[0], 
                                              V1[1] - V2[1], V2[1]);
  // clang-format on

  // If the tangent lines are parallel, winding number is ambiguous
  if(axom::utilities::isNearlyEqual(orient, 0.0)) return 0;

  double dotprod = axom::utilities::clampVal(V1.dot(V2), -1.0, 1.0);

  return -0.5 * M_1_PI * acos(dotprod) * ((orient > 0) ? 1 : -1);
}

/*!
 * \brief Directly compute the winding number at an endpoint of a 
 *        Bezier curve with a convex control polygon
 *
 * \param [in] is_init Boolean value indicating endpoint is at t=0
 * \param [in] c The BezierCurve object to compute the winding number along
 *
 * The winding number for a Bezier curve with a convex control polygon is
 * given by the signed angle between the tangent vector at that endpoint and
 * the vector in the direction of the other endpoint. The query can be located
 * on both endpoints if it is closed.
 * 
 * \return 
 */
template <typename T>
double convex_endpoint_winding_number(bool is_init, const BezierCurve<T, 2>& c)
{
  const int ord = c.getOrder();

  if(ord == 1) return 0;

  Vector<T, 2> V1, V2;
  if(is_init)  // Query point is on initial endpoint
  {
    V1 = Vector<T, 2>(c[0], c[1]).unitVector();
    V2 = Vector<T, 2>(c[0], c[ord]).unitVector();
  }
  else  // Query point is on terminal endpoint
  {
    V1 = Vector<T, 2>(c[ord], c[0]).unitVector();
    V2 = Vector<T, 2>(c[ord], c[ord - 1]).unitVector();
  }

  // clang-format off
  double orient = axom::numerics::determinant(V1[0] - V2[0], V2[0], 
                                              V1[1] - V2[1], V2[1]);
  // clang-format on

  // If the tangent lines are parallel, winding number is ambiguous
  if(axom::utilities::isNearlyEqual(orient, 0.0)) return 0;

  double dotprod = axom::utilities::clampVal(V1.dot(V2), -1.0, 1.0);

  return 0.5 * M_1_PI * acos(dotprod) * ((orient > 0) ? 1 : -1);
}

/*!
 * \brief Recursively compute the winding number for a query point with respect
 *        to a single Bezier curve.
 *
 * \param [in] q The query point at which to compute winding number
 * \param [in] c The BezierCurve object along which to compute the winding number
 * \param [in] convex_cp Boolean flag if the input Bezier curve is already convex
 * \param [in] linear_tol The tolerance level at which a BezierCurve object is to
 *             be interpreted as linear.
 * \param [in] edge_tol The tolerance level at which we consider a query point
 *             to be exactly on the Bezier curve
 *
 * Use a recursive algorithm that checks if the query point is exterior to
 * the convex control polygon of a Bezier curve, in which case we have a direct formula
 * for the winding number. If not, we bisect our curve and run the algorithm on 
 * each half. Use the proximity of the query point to endpoints and approximate
 * linearity of the Bezier curve as base cases.
 * 
 * \return double The winding number.
 */
template <typename T>
double adaptive_winding_number(const Point2D& q,
                               const BezierCurve<T, 2>& c,
                               bool convex_cp,
                               double linear_tol,
                               double edge_tol)
{
  const int ord = c.getOrder();

  double cl_winding_num = closure_winding_number(q, c);
  Polygon<T, 2> controlPolygon(c.getControlPoints());

  // Use linearity as base case for recursion
  if(c.isLinear(linear_tol))
  {
    if(squared_distance(q, Segment<T, 2>(c[0], c[ord])) <= edge_tol)
      return 0;
    else
      return -cl_winding_num;
  }

  // If outside control polygon (with nonzero protocol)
  if(!in_polygon(q, controlPolygon, true, false, edge_tol))
    return -cl_winding_num;

  // Check if our new curve is convex, if we have to
  if(!convex_cp) convex_cp = is_convex(controlPolygon);

  // Can't use endpoint formulas if not convex
  if(convex_cp)
  {
    // Check if we are close to an endpoint
    bool at_init_endpoint = (squared_distance(q, c[0]) <= edge_tol);
    bool at_final_endpoint = (squared_distance(q, c[ord]) <= edge_tol);

    // Use direct formula if at either endpoint.
    if(at_init_endpoint && !at_final_endpoint)
      return convex_endpoint_winding_number(true, c);
    if(at_final_endpoint && !at_init_endpoint)
      return convex_endpoint_winding_number(false, c);
    // If at both endpoints, do a split and avoid a headache
  }

  // Otherwise, our quadrature didn't give us a good enough answer, so we try again
  BezierCurve<T, 2> c1, c2;
  c.split(0.5, c1, c2);

  return adaptive_winding_number(q, c1, convex_cp, linear_tol, edge_tol) +
    adaptive_winding_number(q, c2, convex_cp, linear_tol, edge_tol);
}

}  // end namespace detail
}  // end namespace primal
}  // end namespace axom

#endif
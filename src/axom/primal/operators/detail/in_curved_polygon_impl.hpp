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
template <typename T>
int sgn(T val)
{
  return (T(0) < val) - (val < T(0));
}

/*!
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
  int ord = c.getOrder();
  Vector<T, 2> V1 = Vector<T, 2>(q, c[0]).unitVector();
  Vector<T, 2> V2 = Vector<T, 2>(q, c[ord]).unitVector();

  double dotprod = Vector<T, 2>::dot_product(V1, V2);
  double orient = V2[1] * (V1[0] - V2[0]) - V2[0] * (V1[1] - V2[1]);

  if(dotprod >= 1) return 0;
  if(dotprod <= -1) return -0.5 * sgn(orient);

  return -0.5 * M_1_PI * acos(dotprod) * sgn(orient);
}

/*!
 * \brief Directly compute the winding number at an endpoint of a 
 *        Bezier curve with a convex control polygon
 *
 * \param [in] is_init Boolean value indicating which endpoint to evaluate
 * \param [in] c The BezierCurve object to compute the winding number along
 *
 * The winding number for a Bezier curve with a convex control polygon is
 * given by the signed angle between the tangent vector at that endpoint and
 * the vector in the direction of the other endpoint.
 * 
 * \return 
 */
template <typename T>
double convex_endpoint_winding_number(bool is_init, const BezierCurve<T, 2>& c)
{
  int ord = c.getOrder();

  if(ord == 1) return 0;

  Vector<T, 2> V1, V2;
  if(is_init)
  {
    V1 = c.dt(0).unitVector();
    V2 = Vector<T, 2>(c[0], c[ord]).unitVector();
  }
  else
  {
    V1 = Vector<T, 2>(c[ord], c[0]).unitVector();
    V2 = -c.dt(1).unitVector();
  }

  double orient = V2[1] * (V1[0] - V2[0]) - V2[0] * (V1[1] - V2[1]);

  return 0.5 * M_1_PI * acos(Vector<T, 2>::dot_product(V1, V2)) * sgn(orient);
}

/*!
 * \brief Recursively compute the winding number for a query point with respect
 *        to a single Bezier curve.
 *
 * \param [in] q The query point at which to compute winding number
 * \param [in] c The BezierCurve object along which to compute the winding number
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
                               double linear_tol,
                               double edge_tol)
{
  int ord = c.getOrder();
  double cl_winding_num = closure_winding_number(q, c);
  Polygon<T, 2> controlPolygon(c.getControlPoints());

  // Use direct formulas for endpoints on convex shapes
  if((squared_distance(q, c[0]) <= edge_tol) && is_convex(controlPolygon))
    return convex_endpoint_winding_number(true, c);
  else if((squared_distance(q, c[ord]) <= edge_tol) && is_convex(controlPolygon))
    return convex_endpoint_winding_number(false, c);
  // Use linearity as base case for recursion
  else if(c.isLinear(linear_tol))
  {
    if(squared_distance(q, Segment<T, 2>(c[0], c[ord])) <= edge_tol)
      return 0;
    else
      return -cl_winding_num;
  }
  // Indicates point is far enough away to ensure we are outside the closure
  else if(!in_polygon(q, controlPolygon, edge_tol) && is_convex(controlPolygon))
    return -cl_winding_num;
  // Otherwise, our quadrature didn't give us a good enough answer, so we try again
  else
  {
    BezierCurve<T, 2> c1, c2;
    c.split(0.5, c1, c2);

    return adaptive_winding_number(q, c1, linear_tol, edge_tol) +
      adaptive_winding_number(q, c2, linear_tol, edge_tol);
  }
}

}  // end namespace detail
}  // end namespace primal
}  // end namespace axom

#endif
// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
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
#include "axom/primal/operators/evaluate_integral.hpp"

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
/*
 * \brief Compute the winding number with respect to a line segment
 *
 * \param [in] q The query point to test
 * \param [in] c0 The initial point of the line segment
 * \param [in] c1 The terminal point of the line segment
 * \param [in] edge_tol The tolerance at which a point is on the line
 *
 * The winding number for a point with respect to a straight line
 * is the signed angle subtended by the query point to each endpoint.
 *
 * \return double The winding number
 */
template <typename T>
double linear_winding_number(const Point<T, 2>& q,
                             const Point<T, 2>& c0,
                             const Point<T, 2>& c1,
                             double edge_tol)
{
  Vector<T, 2> V1(q, c0);
  Vector<T, 2> V2(q, c1);

  // clang-format off
  // Measures the signed area of the triangle with vertices q, c0, c1
  double tri_area = axom::numerics::determinant(V1[0] - V2[0], V2[0], 
                                                V1[1] - V2[1], V2[1]);
  // clang-format on

  // Compute distance from line connecting endpoints to query
  if(tri_area * tri_area <= edge_tol * edge_tol * (V1 - V2).squared_norm())
    return 0;

  // Compute signed angle between vectors
  double dotprod = axom::utilities::clampVal(
    Vector<T, 2>::dot_product(V1.unitVector(), V2.unitVector()),
    -1.0,
    1.0);

  return 0.5 * M_1_PI * acos(dotprod) * ((tri_area > 0) ? 1 : -1);
}

/*!
 * \brief Directly compute the winding number at either endpoint of a 
 *        Bezier curve with a convex control polygon
 *
 * \param [in] q The query point
 * \param [in] c The BezierCurve object to compute the winding number along
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance for isNearlyZero
 * \pre Control polygon for c must be convex
 * \pre The query point must be on one of the endpoints
 *
 * The winding number for a Bezier curve with a convex control polygon is
 * given by the signed angle between the tangent vector at that endpoint and
 * the vector in the direction of the other endpoint. 
 * 
 * The query can be located on both endpoints if it is closed, in which case
 * the angle is that between the tangent lines at both endpoints
 * 
 * \return double The winding number
 */
template <typename T>
double convex_endpoint_winding_number(const Point<T, 2>& q,
                                      const BezierCurve<T, 2>& c,
                                      double edge_tol,
                                      double EPS)
{
  const int ord = c.getOrder();
  if(ord == 1) return 0;

  double edge_tol_sq = edge_tol * edge_tol;

  // Verify that the shape is convex, and that the query point is at an endpoint
  SLIC_ASSERT(is_convex(Polygon<T, 2>(c.getControlPoints()), EPS));
  SLIC_ASSERT((squared_distance(q, c[0]) <= edge_tol_sq) ||
              (squared_distance(q, c[ord]) <= edge_tol_sq));

  int idx;

  // Need to find vectors that subtend the entire curve.
  //   We must ignore duplicate nodes
  for(idx = 0; idx <= ord; ++idx)
    if(squared_distance(q, c[idx]) > edge_tol_sq) break;
  Vector<T, 2> V1(q, c[idx]);

  for(idx = ord; idx >= 0; --idx)
    if(squared_distance(q, c[idx]) > edge_tol_sq) break;
  Vector<T, 2> V2(q, c[idx]);

  // clang-format off
  // Measures the signed area of the triangle spanned by V1 and V2
  double tri_area = axom::numerics::determinant(V1[0] - V2[0], V2[0], 
                                                V1[1] - V2[1], V2[1]);
  // clang-format on

  // This means the bounding vectors are anti-parallel.
  //  Parallel tangents can't happen with nontrivial convex control polygons
  if((ord > 3) && axom::utilities::isNearlyEqual(tri_area, 0.0, EPS))
  {
    for(int i = 1; i < ord; ++i)
    {
      // Need to find the first non-parallel control node
      V2 = Vector<T, 2>(q, c[i]);

      // clang-format off
      tri_area = axom::numerics::determinant(V1[0] - V2[0], V2[0], 
                                             V1[1] - V2[1], V2[1]);
      // clang-format on

      // Because we are convex, a single non-collinear vertex tells us the orientation
      if(!axom::utilities::isNearlyEqual(tri_area, 0.0, EPS))
        return (tri_area > 0) ? 0.5 : -0.5;
    }

    // If all vectors are parallel, the curve is linear and return 0
    return 0;
  }

  // Compute signed angle between vectors
  double dotprod = axom::utilities::clampVal(
    Vector<T, 2>::dot_product(V1.unitVector(), V2.unitVector()),
    -1.0,
    1.0);
  return 0.5 * M_1_PI * acos(dotprod) * ((tri_area > 0) ? 1 : -1);
}

/*!
 * \brief Recursively compute the winding number for a query point with respect
 *        to a single Bezier curve.
 *
 * \param [in] q The query point at which to compute winding number
 * \param [in] c The BezierCurve object along which to compute the winding number
 * \param [in] isConvexControlPolygon Boolean flag if the input Bezier curve 
 *                                    is already convex
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance for nonphysical distances, used in
 *   isLinear, isNearlyZero, in_polygon, is_convex
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
double winding_number_recursive(const Point<T, 2>& q,
                                const BezierCurve<T, 2>& c,
                                bool isConvexControlPolygon,
                                double edge_tol = 1e-8,
                                double EPS = 1e-8)
{
  const int ord = c.getOrder();
  if(ord <= 0) return 0.0;  // Catch degenerate cases

  //c.desmos_print(std::cout);

  // Use linearity as base case for recursion
  if(c.isLinear(EPS)) return linear_winding_number(q, c[0], c[ord], edge_tol);

  Polygon<T, 2> controlPolygon(c.getControlPoints());

  // Check if our new curve is convex.
  //  If so, all subcurves will be convex as well
  if(!isConvexControlPolygon)
  {
    isConvexControlPolygon = is_convex(controlPolygon, EPS);
  }
  else  // Formulas for winding number only work if shape is convex
  {
    // If q is outside the control polygon, for an open Bezier curve, the winding
    //  number for the shape connected at the endpoints with straight lines is zero.
    //  We then subtract the contribution of this line segment.
    if(!in_polygon(q, controlPolygon, true, false, EPS))
      return 0.0 - linear_winding_number(q, c[ord], c[0], edge_tol);

    // If the query point is at either endpoint, use direct formula
    if((squared_distance(q, c[0]) <= edge_tol * edge_tol) ||
       (squared_distance(q, c[ord]) <= edge_tol * edge_tol))
      return convex_endpoint_winding_number(q, c, edge_tol, EPS);
  }

  // Recursively split curve until query is outside each control polygon
  BezierCurve<T, 2> c1, c2;
  c.split(0.5, c1, c2);

  return winding_number_recursive(q, c1, isConvexControlPolygon, edge_tol, EPS) +
    winding_number_recursive(q, c2, isConvexControlPolygon, edge_tol, EPS);
}

/*!
 * \brief Compute the orientation of an array of convex Bezier curves
 *
 * \param [in] cvxs The array of convex Bezier curves
 * \param [in] EPS Error tolerance for isNearlyZero and isLinear
 *
 * Bezier curves with Convex control polygons have a constant orientation about
 * around their interior. This method computes that orientation using signed 
 * areas from endpoints to each internal node. 
 * A value of 1 denotes positive orientation, -1 denotes negative orientation,
 * 0 denotes the curve has no interior (i.e. is linear) and thus no orientation
 * 
 * \return Array<int> An above integer for each curve.
 */
template <typename T>
axom::Array<int> convex_orientation(const axom::Array<BezierCurve<T, 2>>& cvxs,
                                    double EPS)
{
  axom::Array<int> orientations;
  for(int i = 0; i < cvxs.size(); i++)
  {
    const BezierCurve<T, 2>& c = cvxs[i];

    // Linear curves have no interior
    if(c.isLinear(EPS))
    {
      orientations.push_back(0);
      break;
    }

    // Get terminal vector for signed area calculation
    Vector<T, 2> V2(c[0], c[c.getOrder()]);

    // Loop over each interior node of the Bezier curve
    int j;
    for(j = 1; j < c.getOrder() - 1; j++)
    {
      Vector<T, 2> V1(c[0], c[j]);

      // clang-format off
      double orient = axom::numerics::determinant(V1[0] - V2[0], V2[0], 
                                                  V1[1] - V2[1], V2[1]);
      // clang-format on

      if(!axom::utilities::isNearlyEqual(orient, 0.0, EPS))
      {
        // Because we are convex, a single non-colinear vertex tells us the orientation
        orientations.push_back((orient > 0) ? 1 : -1);
        break;
      }
      // else, need to try another internal node
    }

    // If we make it through the loop, every node is colinear
    // (although with by a different metric than isLinear above)
    if(j == c.getOrder() - 1) orientations.push_back(0);
  }

  return orientations;
}

/*!
 * \brief Compute the winding number for a query point with respect
 *        to a single Bezier curve directly using quadrature.
 *
 * \param [in] q The query point at which to compute winding number
 * \param [in] cvx The convex BezierCurve object
 * \param [in] cvx_orient Integer value representing orientation of the curve
 * \param [in] quad The MFEM quadrature rule containing 1D nodes and weights
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance for nonphysical distances, used in
 *   isLinear, isNearlyZero, in_polygon, is_convex
 * 
 * Use a direct algorithm evaluates the winding number using quadrature. If the
 * relevant integrand is negative at a quadrature node, then the loop can be 
 * terminated early, as this means the query point is oriented outside a line
 * tangent to the shape at that curve. Otherwise, the quadrature is rounded,
 * and an evaluation of "inside" or "outside" is made accordingly.
 * 
 * \return double The winding number.
 */
template <typename T>
double winding_number_cvx_quad(const Point<T, 2>& q,
                               const BezierCurve<T, 2>& cvx,
                               int cvx_orient,
                               const mfem::IntegrationRule& quad,
                               double edge_tol,
                               double EPS)
{
  const int ord = cvx.getOrder();

  //std::cout << q << ", x_n = [";
  //for(int p = 0; p <= ord; ++p)
  //{
  //  std::cout << cvx[p][0] << (p < ord ? "," : "");
  //}
  //std::cout << "], y_n = [";
  //for(int p = 0; p <= ord; ++p)
  //{
  //  std::cout << cvx[p][1] << (p < ord ? "," : "");
  //}
  //std::cout << "]" << std::endl;

  // This number plus the winding number is either 0 or \pm 1
  double closure_winding_num =
    linear_winding_number(q, cvx[ord], cvx[0], edge_tol);

  // This flag means the curve is flat, and we use the linear winding number
  if(cvx_orient == 0) -closure_winding_num;

  // If the closure winding number and the orientation have opposite signs,
  //  we are external to the closed shape
  if(cvx_orient * closure_winding_num < 0) return 0.0 - closure_winding_num;

  // Evaluate integrand at each quadrature point
  //   Return early if direct answer can be found
  double quad_result = 0;
  for(int k = 0; k < quad.GetNPoints(); k++)
  {
    Point<T, 2> x_quad = cvx.evaluate(quad.IntPoint(k).x);
    Vector<T, 2> dx_quad = cvx.dt(quad.IntPoint(k).x);

    // Compute denominator for quadrature:
    double query_dist = squared_distance(x_quad, q);

    // Early return: Use direct formula if query point is nearly on quadrature node,
    //  Direct formula consistent with `convex_endpoint_winding_number`
    if(query_dist <= edge_tol * edge_tol)
      return cvx_orient *
        (0.5 - linear_winding_number(q, cvx[0], cvx[ord], edge_tol));

    // Compute numerator for quadrature
    // clang-format off
    double query_orient = axom::numerics::determinant(x_quad[0] - q[0], dx_quad[0],
                                                      x_quad[1] - q[1], dx_quad[1]);
    // clang-format on

    // Early return: If numerator is ever negative with positive orientation,
    //  or positive with negative orientation, we are exterior to the closed curve
    if(cvx_orient * query_orient < 0) return 0.0 - closure_winding_num;

    // Add result to quadrature as normal
    quad_result += quad.IntPoint(k).weight * query_orient / query_dist;
  }

  // If rounded total winding number would be +/-1, assume we are inside
  if(std::abs(0.5 * M_1_PI * quad_result + closure_winding_num) >= 0.5)
    return cvx_orient - closure_winding_num;
  else  // else, we are outside the closed shape
    return 0.0 - closure_winding_num;
}

}  // end namespace detail
}  // end namespace primal
}  // end namespace axom

#endif

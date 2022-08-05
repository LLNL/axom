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
#include "axom/primal/operators/evaluate_integral.hpp"

// MFEM includes
#ifdef AXOM_USE_MFEM
  #include "mfem.hpp"
#else
  #error "Primal's in/out functions for CurvedPolygon require mfem library."
#endif

// C++ includes
#include <cmath>
#include <iomanip>

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
 * winding number of the query point with respect to this segment. This
 * is calculated directly by measuring the angle spanned by lines connecting
 * the query point to each endpoint of the Bezier curve.
 * 
 * \return 
 */
template <typename T>
double closure_winding_number(const Point<T, 2>& q,
                              const BezierCurve<T, 2>& c,
                              const double edge_tol)
{
  const int ord = c.getOrder();
  Vector<T, 2> V1 = Vector<T, 2>(q, c[0]);
  Vector<T, 2> V2 = Vector<T, 2>(q, c[ord]);

  // clang-format off
  double orient = axom::numerics::determinant(V1[0] - V2[0], V2[0], 
                                              V1[1] - V2[1], V2[1]);
  // clang-format on

  if(orient * orient / (V1 - V2).squared_norm() <= edge_tol * edge_tol)
    return 0;

  double dotprod = axom::utilities::clampVal(
    Vector<T, 2>::dot_product(V1.unitVector(), V2.unitVector()),
    -1.0,
    1.0);

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
double convex_endpoint_winding_number(bool is_init,
                                      const BezierCurve<T, 2>& c,
                                      const double EPS)
{
  const int ord = c.getOrder();

  if(ord == 1) return 0;

  Vector<T, 2> V1, V2;
  // Using 1e-16 and 1 - 1e-16 is a shortgap measure until
  //  we can get those tangent lines efficiently.
  if(is_init)  // Query point is on initial endpoint
  {
    V1 = c.dt(1e-16);
    V2 = Vector<T, 2>(c[0], c[ord]);
  }
  else  // Query point is on terminal endpoint
  {
    V1 = Vector<T, 2>(c[ord], c[0]);
    V2 = -c.dt(1 - 1e-16);
  }

  // clang-format off
  double orient = axom::numerics::determinant(V1[0] - V2[0], V2[0], 
                                              V1[1] - V2[1], V2[1]);
  // clang-format on

  // If the tangent lines are parallel, winding number is 0
  if(axom::utilities::isNearlyEqual(orient, 0.0, EPS))
  {
    return 0;
  }

  double dotprod = axom::utilities::clampVal(
    Vector<T, 2>::dot_product(V1.unitVector(), V2.unitVector()),
    -1.0,
    1.0);
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
                               double edge_tol = 1e-8,
                               double EPS = 1e-8)
{
  const int ord = c.getOrder();

  //std::cout << std::setprecision(16) << q << ", x_n = [";
  //for(int p = 0; p <= ord; ++p)
  //{
  //  std::cout << c[p][0] << (p < ord ? "," : "");
  //}
  //std::cout << "], y_n = [";
  //for(int p = 0; p <= ord; ++p)
  //{
  //  std::cout << c[p][1] << (p < ord ? "," : "");
  //}
  //std::cout << "]" << std::endl;

  Polygon<T, 2> controlPolygon(c.getControlPoints());

  // Use linearity as base case for recursion
  if(c.isLinear(EPS)) return -closure_winding_number(q, c, edge_tol);

  // If outside control polygon (with nonzero protocol)
  if(!in_polygon(q, controlPolygon, true, false, EPS))
    return -closure_winding_number(q, c, edge_tol);

  // Check if our new curve is convex, if we have to
  if(!convex_cp) convex_cp = is_convex(controlPolygon);

  // Can't use endpoint formulas if not convex
  if(convex_cp)
  {
    // Check if we are close to an endpoint
    bool at_init_endpoint = (squared_distance(q, c[0]) <= edge_tol * edge_tol);
    bool at_final_endpoint = (squared_distance(q, c[ord]) <= edge_tol * edge_tol);

    // Use direct formula if at either endpoint.
    if(at_init_endpoint && !at_final_endpoint)
      return convex_endpoint_winding_number(true, c, EPS);
    if(at_final_endpoint && !at_init_endpoint)
      return convex_endpoint_winding_number(false, c, EPS);
    // If at both endpoints, do a split and avoid a headache
  }

  // Otherwise, our quadrature didn't give us a good enough answer, so we try again
  BezierCurve<T, 2> c1, c2;
  c.split(0.5, c1, c2);

  return adaptive_winding_number(q, c1, convex_cp, edge_tol, EPS) +
    adaptive_winding_number(q, c2, convex_cp, edge_tol, EPS);
}

// Get the function to be passed into the evaluate integral function
inline std::function<Vector2D(Point2D)> get_winding_func(Point2D p)
{
  return [p](Point2D x) -> Vector2D {
    double denom =
      2 * M_PI * ((x[0] - p[0]) * (x[0] - p[0]) + (x[1] - p[1]) * (x[1] - p[1]));
    return Vector2D({-(x[1] - p[1]) / 1.0, (x[0] - p[0]) / denom});
  };
}

template <typename T>
axom::Array<int> convex_orientation(const axom::Array<BezierCurve<T, 2>>& cvxs)
{
  axom::Array<int> orientations;
  for(int i = 0; i < cvxs.size(); i++)
  {
    const BezierCurve<T, 2>& c = cvxs[i];

    if(c.isLinear())
    {
      orientations.push_back(0);
      break;
    }

    // Loop over each interior node of the Bezier curve
    int j;
    for(j = 1; j < c.getOrder() - 1; j++)
    {
      Vector<T, 2> V1(c[0], c[j]);
      Vector<T, 2> V2(c[0], c[c.getOrder()]);

      // clang-format off
      double orient = axom::numerics::determinant(V1[0] - V2[0], V2[0], 
                                                  V1[1] - V2[1], V2[1]);
      // clang-format on

      if(!axom::utilities::isNearlyEqual(orient, 0.0))
      {
        // Because we are convex, a single non-colinear vertex tells us the orientation
        orientations.push_back((orient > 0) ? 1 : -1);
        break;
      }
      // else, need to try another internal node
    }

    // If we make it through the loop, every node is colinear
    if( j == c.getOrder() - 1 )
      orientations.push_back(0);
  }

  return orientations;
}

// Get total winding number with respect to each curve
template <typename T>
double convex_winding_number(const Point<T, 2>& q,
                             const BezierCurve<T, 2>& cvx,
                             const int& cvx_orient,
                             const mfem::IntegrationRule& quad,
                             const double edge_tol = 1e-8)
{
  // Zero orientation indicates the curve is flat
  if(cvx.isLinear(0)) return -closure_winding_number(q, cvx, edge_tol);

  const int ord = cvx.getOrder();

  //std::cout << std::setprecision(16) << q << ", x_n = [";
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

  double quad_result = 0;

  for(int k = 0; k < quad.GetNPoints(); k++)
  {
    Point<T, 2> x_q = cvx.evaluate(quad.IntPoint(k).x);
    Vector<T, 2> dx_q = cvx.dt(quad.IntPoint(k).x);

    // Denominator for quadrature
    double q_dist = squared_distance(x_q, q);
    if(q_dist <= edge_tol * edge_tol)
      return cvx_orient * 0.5 - closure_winding_number(q, cvx, edge_tol);

    // Numerator for quadrature
    // clang-format off
    double q_orient = axom::numerics::determinant(x_q[0] - q[0], dx_q[0],
                                                  x_q[1] - q[1], dx_q[1]);
    // clang-format on
    if(q_orient < 0) return -closure_winding_number(q, cvx, edge_tol);

    quad_result += quad.IntPoint(k).weight * q_orient / q_dist;
  }

  double closure_number = closure_winding_number(q, cvx, edge_tol);

  // Make a call based off of the quadrature
  if((0.5 * M_1_PI * quad_result) + closure_number >= 0.5)
    return 1.0 - closure_number;
  return 0.0 - closure_number;
}

}  // end namespace detail
}  // end namespace primal
}  // end namespace axom

#endif

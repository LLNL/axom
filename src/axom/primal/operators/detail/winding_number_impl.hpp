// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef PRIMAL_WINDING_NUMBER_IMPL_HPP_
#define PRIMAL_WINDING_NUMBER_IMPL_HPP_

// Axom includes
#include "axom/config.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/Polygon.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/operators/in_polygon.hpp"
#include "axom/primal/operators/is_convex.hpp"
#include "axom/primal/operators/squared_distance.hpp"

// C++ includes
#include <math.h>

// MFEM includes
#ifdef AXOM_USE_MFEM
  #include "mfem.hpp"
#endif

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
 * Colinear points return 0 for their winding number.
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
  {
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
  if(ord == 1)
  {
    return 0;
  }

  double edge_tol_sq = edge_tol * edge_tol;

  // Verify that the shape is convex, and that the query point is at an endpoint
  SLIC_ASSERT(is_convex(Polygon<T, 2>(c.getControlPoints()), EPS));
  SLIC_ASSERT((squared_distance(q, c[0]) <= edge_tol_sq) ||
              (squared_distance(q, c[ord]) <= edge_tol_sq));

  int idx;

  // Need to find vectors that subtend the entire curve.
  //   We must ignore duplicate nodes
  for(idx = 0; idx <= ord; ++idx)
  {
    if(squared_distance(q, c[idx]) > edge_tol_sq)
    {
      break;
    }
  }
  Vector<T, 2> V1(q, c[idx]);

  for(idx = ord; idx >= 0; --idx)
  {
    if(squared_distance(q, c[idx]) > edge_tol_sq)
    {
      break;
    }
  }
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
      {
        return (tri_area > 0) ? 0.5 : -0.5;
      }
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
                                      is already convex
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance for nonphysical distances, used in
 *   isLinear, isNearlyZero, in_polygon, is_convex
 * 
 * Use a recursive algorithm that checks if the query point is exterior to
 * some convex shape containing the Bezier curve, in which case we have a direct 
 * formula for the winding number. If not, we bisect our curve and run the algorithm on 
 * each half. Use the proximity of the query point to endpoints and approximate
 * linearity of the Bezier curve as base cases.
 * 
 * \return double The winding number.
 */
template <typename T>
double curve_winding_number_recursive(const Point<T, 2>& q,
                                      const BezierCurve<T, 2>& c,
                                      bool isConvexControlPolygon,
                                      double edge_tol = 1e-8,
                                      double EPS = 1e-8)
{
  const int ord = c.getOrder();
  if(ord <= 0)
  {
    return 0.0;  // Catch degenerate cases
  }

  // If q is outside a convex shape that contains the entire curve, the winding
  //   number for the shape connected at the endpoints with straight lines is zero.
  //   We then subtract the contribution of this line segment.

  // Simplest convex shape containing c is its bounding box
  BoundingBox<T, 2> bBox(c.boundingBox());
  if(!bBox.contains(q))
  {
    return 0.0 - linear_winding_number(q, c[ord], c[0], edge_tol);
  }

  // Use linearity as base case for recursion.
  if(c.isLinear(EPS))
  {
    return linear_winding_number(q, c[0], c[ord], edge_tol);
  }

  // Check if our control polygon is convex.
  //  If so, all subsequent control polygons will be convex as well
  Polygon<T, 2> controlPolygon(c.getControlPoints());
  const bool includeBoundary = true;
  const bool useNonzeroRule = true;

  if(!isConvexControlPolygon)
  {
    isConvexControlPolygon = is_convex(controlPolygon, EPS);
  }
  else  // Formulas for winding number only work if shape is convex
  {
    // Bezier curves are always contained in their convex control polygon
    if(!in_polygon(q, controlPolygon, includeBoundary, useNonzeroRule, EPS))
    {
      return 0.0 - linear_winding_number(q, c[ord], c[0], edge_tol);
    }

    // If the query point is at either endpoint, use direct formula
    if((squared_distance(q, c[0]) <= edge_tol * edge_tol) ||
       (squared_distance(q, c[ord]) <= edge_tol * edge_tol))
    {
      return convex_endpoint_winding_number(q, c, edge_tol, EPS);
    }
  }

  // Recursively split curve until query is outside some known convex region
  BezierCurve<T, 2> c1, c2;
  c.split(0.5, c1, c2);

  return curve_winding_number_recursive(q, c1, isConvexControlPolygon, edge_tol, EPS) +
    curve_winding_number_recursive(q, c2, isConvexControlPolygon, edge_tol, EPS);
}

/// Type to indicate orientation of singularities relative to surface
enum class SingularityAxis
{
  x,
  y,
  z,
  rotated
};

#ifdef AXOM_USE_MFEM
/*!
 * \brief Evaluates an "anti-curl" of the winding number along a curve
 *
 * \param [in] query The query point to test
 * \param [in] curve The BezierCurve object
 * \param [in] ax The axis (relative to query) denoting which anti-curl we use
 * \param [in] npts The number of points used in each Gaussian quadrature
 * \param [in] quad_tol The maximum relative error allowed in each quadrature
 * 
 * Applies a non-adaptive quadrature to a BezierCurve using one of three possible
 * "anti-curl" vector fields, the curl of each of which is equal to <x, y, z>/||x||^3.
 * With the proper "anti-curl" selected, integrating this along a closed curve is equal
 * to evaluating the winding number of a surface with that curve as the boundary.
 * 
 * \note This is only meant to be used for `winding_number<BezierPatch>()`,
 *  and the result does not make sense outside of that context.
 *
 * \return double One component of the winding number
 */
template <typename T>
double stokes_winding_number(const Point<T, 3>& query,
                             const BezierCurve<T, 3>& curve,
                             const SingularityAxis ax,
                             int npts,
                             double quad_tol)
{
  // Generate the quadrature rules in parameter space
  static mfem::IntegrationRules my_IntRules(0, mfem::Quadrature1D::GaussLegendre);
  const mfem::IntegrationRule& quad_rule =
    my_IntRules.Get(mfem::Geometry::SEGMENT, 2 * npts - 1);

  double quadrature = 0.0;
  for(int q = 0; q < quad_rule.GetNPoints(); ++q)
  {
    // Get quadrature points in space (shifted by the query)
    const Vector<T, 3> node(query, curve.evaluate(quad_rule.IntPoint(q).x));
    const Vector<T, 3> node_dt(curve.dt(quad_rule.IntPoint(q).x));
    const double node_norm = node.norm();

    // Compute one of three vector field line integrals depending on
    //  the orientation of the original surface, indicated through ax.
    switch(ax)
    {
    case(SingularityAxis::x):
      quadrature += quad_rule.IntPoint(q).weight *
        (node[2] * node[0] * node_dt[1] - node[1] * node[0] * node_dt[2]) /
        (node[1] * node[1] + node[2] * node[2]) / node_norm;
      break;
    case(SingularityAxis::y):
      quadrature += quad_rule.IntPoint(q).weight *
        (node[0] * node[1] * node_dt[2] - node[2] * node[1] * node_dt[0]) /
        (node[0] * node[0] + node[2] * node[2]) / node_norm;
      break;
    case(SingularityAxis::z):
    case(SingularityAxis::rotated):
      quadrature += quad_rule.IntPoint(q).weight *
        (node[1] * node[2] * node_dt[0] - node[0] * node[2] * node_dt[1]) /
        (node[0] * node[0] + node[1] * node[1]) / node_norm;
      break;
    }
  }

  // Adaptively refine quadrature over curves if query is not far enough away
  //  from the singularity axis. If rotated, assume you need to adapt.
  bool needs_adapt = false;
  BoundingBox<T, 3> cBox(curve.boundingBox());
  Point<T, 3> centroid = cBox.getCentroid();

  switch(ax)
  {
  case(SingularityAxis::x):
    needs_adapt = (query[1] - centroid[1]) * (query[1] - centroid[1]) +
        (query[2] - centroid[2]) * (query[2] - centroid[2]) <=
      cBox.range().squared_norm();
    break;
  case(SingularityAxis::y):
    needs_adapt = (query[0] - centroid[0]) * (query[0] - centroid[0]) +
        (query[2] - centroid[2]) * (query[2] - centroid[2]) <=
      cBox.range().squared_norm();
    break;
  case(SingularityAxis::z):
    needs_adapt = (query[0] - centroid[0]) * (query[0] - centroid[0]) *
        (query[1] - centroid[1]) * (query[1] - centroid[1]) <=
      cBox.range().squared_norm();
    break;
  case(SingularityAxis::rotated):
    needs_adapt = true;
    break;
  }

  if(needs_adapt)
  {
    return stokes_winding_number_adaptive(query,
                                          curve,
                                          ax,
                                          quad_rule,
                                          quadrature,
                                          quad_tol);
  }

  return 0.25 * M_1_PI * quadrature;
}
#endif

#ifdef AXOM_USE_MFEM
/*!
 * \brief Recursively evaluates an "anti-curl" of the winding number on subcurves
 *
 * \param [in] query The query point to test
 * \param [in] curve The BezierCurve object
 * \param [in] ax The axis (relative to query) denoting which anti-curl we use
 * \param [in] quad_rule The mfem quadrature rule object
 * \param [in] quad_coarse The integral evaluated on the original curve
 * \param [in] quad_tol The maximum relative error allowed in each quadrature
 * \param [in] depth The current recursive depth
 * 
 * Recursively apply quadrature for one of three possible integrals along two halfs 
 * of a curve. The sum of this integral along the subcurves should be equal to to
 * quad_coarse. Otherwise, apply this algorithm to each half recursively.
 *
 * \note This is only meant to be used for `winding_number<BezierPatch>()`,
 *  and the result does not make sense outside of that context.
 * 
 * \return double One component of the winding number
 */
template <typename T>
double stokes_winding_number_adaptive(const Point<T, 3>& query,
                                      const BezierCurve<T, 3>& curve,
                                      const SingularityAxis ax,
                                      const mfem::IntegrationRule& quad_rule,
                                      const double quad_coarse,
                                      const double quad_tol,
                                      const int depth = 1)
{
  // Split the curve, do the quadrature over both components
  BezierCurve<T, 3> subcurves[2];
  curve.split(0.5, subcurves[0], subcurves[1]);

  double quad_fine[2] = {0.0, 0.0};
  for(int i = 0; i < 2; ++i)
  {
    for(int q = 0; q < quad_rule.GetNPoints(); ++q)
    {
      // Get quad_rulerature points in space (shifted by the query)
      const Vector<T, 3> node(query,
                              subcurves[i].evaluate(quad_rule.IntPoint(q).x));
      const Vector<T, 3> node_dt(subcurves[i].dt(quad_rule.IntPoint(q).x));
      const double node_norm = node.norm();

      // Compute one of three vector field line integrals depending on
      //  the orientation of the original surface, indicated through ax.
      switch(ax)
      {
      case(SingularityAxis::x):
        quad_fine[i] += quad_rule.IntPoint(q).weight *
          (node[2] * node[0] * node_dt[1] - node[1] * node[0] * node_dt[2]) /
          (node[1] * node[1] + node[2] * node[2]) / node_norm;
        break;
      case(SingularityAxis::y):
        quad_fine[i] += quad_rule.IntPoint(q).weight *
          (node[0] * node[1] * node_dt[2] - node[2] * node[1] * node_dt[0]) /
          (node[0] * node[0] + node[2] * node[2]) / node_norm;
        break;
      case(SingularityAxis::z):
      case(SingularityAxis::rotated):
        quad_fine[i] += quad_rule.IntPoint(q).weight *
          (node[1] * node[2] * node_dt[0] - node[0] * node[2] * node_dt[1]) /
          (node[0] * node[0] + node[1] * node[1]) / node_norm;
        break;
      }
    }
  }

  constexpr int MAX_DEPTH = 12;
  if(depth >= MAX_DEPTH ||
     axom::utilities::isNearlyEqualRelative(quad_fine[0] + quad_fine[1],
                                            quad_coarse,
                                            quad_tol,
                                            1e-10))
  {
    return 0.25 * M_1_PI * (quad_fine[0] + quad_fine[1]);
  }
  else
  {
    return stokes_winding_number_adaptive(query,
                                          subcurves[0],
                                          ax,
                                          quad_rule,
                                          quad_fine[0],
                                          quad_tol,
                                          depth + 1) +
      stokes_winding_number_adaptive(query,
                                     subcurves[1],
                                     ax,
                                     quad_rule,
                                     quad_fine[1],
                                     quad_tol,
                                     depth + 1);
  }
}
#endif

}  // end namespace detail
}  // end namespace primal
}  // end namespace axom

#endif

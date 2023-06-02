// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file winding_number.hpp
 *
 * \brief Consists of methods to compute winding numbers for points 
 *        with respect to various geometric objects.
 */

#ifndef AXOM_PRIMAL_WINDING_NUMBER_HPP_
#define AXOM_PRIMAL_WINDING_NUMBER_HPP_

// Axom includes
#include "axom/config.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/geometry/CurvedPolygon.hpp"
#include "axom/primal/operators/detail/winding_number_impl.hpp"

// C++ includes
#include <cmath>

namespace axom
{
namespace primal
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
double winding_number(const Point<T, 2>& q, const Segment<T, 2>& s, double edge_tol)
{
  return linear_winding_number(q, s[0], s[1], edge_tol);
}

/*!
 * \brief Computes the winding number for a point and a polygon
 *
 * \param [in] R The query point to test
 * \param [in] P The Polygon object to test for containment
 * \param [in] useStrictInclusion If true, points on the boundary are considered exterior.
 * \param [in] EPS The tolerance level for collinearity
 * 
 * Uses an adapted ray-casting approach that counts quarter-rotation
 * of vertices around the query point. Current policy is to return 1 on edges
 * without strict inclusion, 0 on edges with strict inclusion.
 *
 * Directly uses algorithm in 
 * Kai Hormann, Alexander Agathos, "The point in polygon problem for arbitrary polygons"
 * Computational Geometry, Volume 20, Issue 3, 2001,
 * 
 * \return The integer winding number
 */
template <typename T>
int winding_number(const Point<T, 2>& R,
                   const Polygon<T, 2>& P,
                   bool useStrictInclusion = false,
                   double EPS = 1e-8)
{
  const int nverts = P.numVertices();

  // If the query is a vertex, return a value interpreted
  //  as "inside" by evenodd or nonzero protocols
  if(axom::utilities::isNearlyEqual(P[0][0], R[0], EPS) &&
     axom::utilities::isNearlyEqual(P[0][1], R[1], EPS))
    return !useStrictInclusion;

  int winding_num = 0;
  for(int i = 0; i < nverts; i++)
  {
    int j = (i == nverts - 1) ? 0 : i + 1;

    if(axom::utilities::isNearlyEqual(P[j][1], R[1], EPS))
    {
      if(axom::utilities::isNearlyEqual(P[j][0], R[0], EPS))
        return !useStrictInclusion;  // On vertex
      else if(P[i][1] == R[1] && ((P[j][0] > R[0]) == (P[i][0] < R[0])))
        return !useStrictInclusion;  // On horizontal edge
    }

    // Check if edge crosses horizontal line
    if((P[i][1] < R[1]) != (P[j][1] < R[1]))
    {
      double det;
      if(P[i][0] >= R[0])
      {
        if(P[j][0] > R[0])
          winding_num += 2 * (P[j][1] > P[i][1]) - 1;
        else
        {
          // clang-format off
          det = axom::numerics::determinant(P[i][0] - R[0], P[j][0] - R[0],
                                            P[i][1] - R[1], P[j][1] - R[1]);
          // clang-format on

          // On edge
          if(axom::utilities::isNearlyEqual(det, 0.0, EPS))
            return !useStrictInclusion;

          // Check if edge intersects horitonal ray to the right of R
          if((det > 0) == (P[j][1] > P[i][1]))
            winding_num += 2 * (P[j][1] > P[i][1]) - 1;
        }
      }
      else
      {
        if(P[j][0] > R[0])
        {
          // clang-format off
          det = axom::numerics::determinant(P[i][0] - R[0], P[j][0] - R[0],
                                            P[i][1] - R[1], P[j][1] - R[1]);
          // clang-format on

          // On edge
          if(axom::utilities::isNearlyEqual(det, 0.0, EPS))
            return !useStrictInclusion;

          // Check if edge intersects horitonal ray to the right of R
          if((det > 0) == (P[j][1] > P[i][1]))
            winding_num += 2 * (P[j][1] > P[i][1]) - 1;
        }
      }
    }
  }

  return winding_num;
}

/*!
 * \brief Computes the generalized winding number for a single Bezier curve
 *
 * \param [in] query The query point to test
 * \param [in] cpoly The Bezier curve object 
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 *
 * Computes the winding number using a recursive, bisection algorithm,
 * using nearly-linear Bezier curves as a base case.
 * 
 * \return float the generalized winding number.
 */
template <typename T>
double winding_number(const Point<T, 2>& q,
                      const BezierCurve<T, 2>& c,
                      double edge_tol = 1e-8,
                      double EPS = 1e-8)
{
  return detail::curve_winding_number_recursive(q, c, false, edge_tol, EPS);
}

/*!
 * \brief Computes the generalized winding number for a curved polygon
 *
 * \param [in] query The query point to test
 * \param [in] cpoly The CurvedPolygon object
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
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
                      double edge_tol = 1e-8,
                      double EPS = 1e-8)
{
  double ret_val = 0.0;
  for(int i = 0; i < cpoly.numEdges(); i++)
    ret_val +=
      detail::curve_winding_number_recursive(q, cpoly[i], false, edge_tol, EPS);

  return ret_val;
}

/*!
 * \brief Computes the solid angle winding number for a 3D triangle
 *
 * \param [in] query The query point to test
 * \param [in] tri The Triangle object
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 *
 * Computes the winding number using the formula from [Barill 2018]
 * 
 * \return float the generalized winding number.
 */
template <typename T>
double winding_number(const Point<T, 3>& q,
                      const Triangle<T, 3>& tri,
                      const double edge_tol = 1e-8,
                      const double EPS = 1e-8)
{
  if(tri.area() == 0) return 0;

  Vector<T, 3> a(q, tri[0]), b(q, tri[1]), c(q, tri[2]);

  // Compute norms. Possibly return early
  double a_norm = a.norm(), b_norm = b.norm(), c_norm = c.norm();
  if(a_norm < edge_tol || b_norm < edge_tol || c_norm < edge_tol) return 0;

  double num = Vector<T, 3>::scalar_triple_product(a, b, c);
  if(axom::utilities::isNearlyEqual(num, 0.0, EPS)) return 0;

  double denom = a_norm * b_norm * c_norm +
    a_norm * Vector<T, 3>::dot_product(b, c) +
    b_norm * Vector<T, 3>::dot_product(a, c) +
    c_norm * Vector<T, 3>::dot_product(a, b);

  return 0.5 * M_1_PI * atan(num / denom);
}

/*!
 * \brief Computes the solid angle winding number for a 3D planar polygon
 *
 * \param [in] query The query point to test
 * \param [in] poly The Polygon object
 * \param [in] edge_tol The physical distance level at which objects are 
 *                      considered indistinguishable
 * \param [in] EPS Miscellaneous numerical tolerance level for nonphysical distances
 * 
 * \pre Assumes the polygon is truly planar
 * Computes the winding number using the formula from [Paeth 1995]
 * 
 * \return float the generalized winding number.
 */
template <typename T>
double winding_number(const Point<T, 3>& q,
                      const Polygon<T, 3>& poly,
                      const double edge_tol = 1e-8,
                      const double EPS = 1e-8)
{
  const int num_verts = poly.numVertices();
  if(num_verts < 3) return 0;

  // Declare variables
  Vector<T, 3> n1, n2, r1, normal = poly.normal().unitVector();
  double cos_ang, ang, signedVolume, area = 0;
  int idx, max_idx = num_verts - 1, skipped_verts = 0;

  int far_idx = 0;
  double curr_dist, far_dist = 0;

  // Find the last vertex that isn't coincident
  Vector<T, 3> v1(poly[0], poly[max_idx]), v2;
  while(v1.norm() <= edge_tol && max_idx > 0)
  {
    --max_idx;
    v1 = Vector<T, 3>(poly[0], poly[max_idx]);
    ++skipped_verts;
  }

  // Catch the very degenerate case (zero perimeter)
  if(max_idx == 0) return 0;

  // Compute the solid angle at each vertex
  for(int i = 0; i <= max_idx; ++i)
  {
    // Compute vector denoting radius of the sphere
    r1 = Vector<T, 3>(q, poly[i]);

    // Catch case where query is on an endpoint, and find furthest index from query
    curr_dist = r1.norm();
    if(curr_dist <= edge_tol) return 0.0;
    if(curr_dist > far_dist)
    {
      far_idx = i;
      far_dist = curr_dist;
    }

    // Find the next vertex that isn't coincident, advance the loop
    idx = i + 1;
    v2 = Vector<T, 3>(poly[i], poly[idx % num_verts]);
    while(v2.norm() <= edge_tol)
    {
      ++idx;
      v2 = Vector<T, 3>(poly[i], poly[idx % num_verts]);
      ++skipped_verts;
    }
    i = idx - 1;

    // Skip the vertex if v1 and v2 are colinear
    signedVolume = Vector<T, 3>::scalar_triple_product(v2, v1, normal);
    if(axom::utilities::isNearlyEqual(signedVolume, 0.0, EPS))
    {
      ++skipped_verts;
      continue;
    }

    // Compute normal vectors to planes that contain great circles corresponding
    //  to sides of the polygon
    n1 = Vector<T, 3>::cross_product(v1, r1);
    n2 = Vector<T, 3>::cross_product(r1, v2);

    // Find the angle between the two planes
    cos_ang = Vector<T, 3>::dot_product(n1, n2) / n1.norm() / n2.norm();
    ang = M_1_PI * acos(axom::utilities::clampVal(cos_ang, -1.0, 1.0));

    // Adjust for if angle is convex or concave after projection.
    area += signedVolume > 0.0 ? 1.0 - ang : 1.0 + ang;

    v1 = -v2;
  }

  // Check that the point isn't coplanar with the polygon
  if(axom::utilities::isNearlyEqual(
       Vector<T, 3>::dot_product(normal, Vector<T, 3>(q, poly[far_idx])),
       0.0,
       edge_tol))
    return 0;

  // Convert from solid angle to winding number, subtract excess
  area = 0.25 * (area - (num_verts - skipped_verts - 2));

  // Line differs from Paeth, flips sign.
  // Should account for the orientation of the polygon
  return (Vector<T, 3>::dot_product(normal, r1) > 0.0) ? area : -area;
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_WINDING_NUMBER_H_

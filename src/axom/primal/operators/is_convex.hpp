// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_IS_CONVEX_HPP_
#define AXOM_PRIMAL_IS_CONVEX_HPP_

#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/Polygon.hpp"
#include "axom/primal/operators/orientation.hpp"
#include "axom/primal/operators/squared_distance.hpp"

namespace axom
{
namespace primal
{
/*!
 * \brief Determines if a polygon defined by ordered vertices is convex
 * 
 * \param [in] poly The polygon
 * 
 * A 2D polygon is convex when every line that does not contain an edge
 * intersects the shape at most twice.
 * Checks whether for each pair of vertices P[i-1]P[i+1], the point
 * P[i] and (P[0] or P[n]) lie on the same side of the line connecting them.
 * 
 * Algorithm adapted from:
 * Y. L. Ma, W.T. Hewitt. "Point inversion and projection for NURBS curve 
 * and surface: Control polygon approach"
 * Computer Aided Geometric Design 20(2):79-99, May 2003.
 * \note Only defined in 2D
 * 
 * \return A boolean value indicating convexity
 */
template <typename T>
bool is_convex(const Polygon<T, 2>& poly, double EPS = 1e-8)
{
  int n = poly.numVertices() - 1;
  if(n + 1 < 3)
  {
    return true;  // Triangles and lines are convex
  }

  for(int i = 1; i < n; i++)
  {
    // For each non-endpoint, check if that point and one of the endpoints
    //  are on the same side as the segment connecting the adjacent nodes
    Segment<T, 2> seg(poly[i - 1], poly[i + 1]);
    int res1 = orientation(poly[i], seg, EPS);

    // Edge case
    if(res1 == primal::ON_BOUNDARY)
    {
      continue;
    }

    // Ensure other point to check against isn't adjacent
    if(res1 == orientation(poly[(i <= n / 2) ? n : 0], seg, EPS))
    {
      return false;
    }
  }

  return true;
}

/*!
 * \brief Determines if a 3D polygon defined by ordered vertices is convex
 * 
 * \param [in] poly The polygon
 * 
 * Uses dot products to detect whether vertices extend in the "convex" direction.
 * Uses the edge P[0]P[N] as a reference, meaning its adjacent edges are
 * always considered to be oriented correctly.
 * 
 * Algorithm adapted from:
 * Y. L. Ma, W.T. Hewitt. "Point inversion and projection for NURBS curve 
 * and surface: Control polygon approach"
 * Computer Aided Geometric Design 20(2):79-99, May 2003.
 * \note Only defined in 2D
 * 
 * \return A boolean value indicating convexity
 */
template <typename T>
bool is_convex(const Polygon<T, 3>& poly, double EPS = 1e-8)
{
  int n = poly.numVertices() - 1;
  if(n + 1 < 3)
  {
    return true;  // Triangles and lines are convex
  }

  for(int i = 1; i < n; i++)
  {
    // For each non-endpoint, check if that point and one of the endpoints
    //  are on the same side as the segment connecting the adjacent nodes
    const Vector<T, 3> v0(poly[i - 1], poly[i]);
    const Vector<T, 3> v1(poly[i - 1], poly[i + 1]);
    const Vector<T, 3> v2(poly[i - 1], poly[(i <= n / 2) ? n : 0]);

    // Equivalent to (v1 x v0) dot (v1 x v2) > 0
    if(v1.squared_norm() * v0.dot(v2) - v1.dot(v2) * v1.dot(v0) > EPS)
    {
      return false;
    }
  }

  return true;
}

/* 
 * Check if a convex polygon is "shallow," which we define as every line
 * that is normal to an edge of the polygon intersects only that edge
 * and the edge P[n]P[0].
 * 
 * This is true when ang(P[n]P[0]P[1]) + ang(P[n-1]P[n]P[0]) < 90,
 * and every interior angle is greater than 90 degrees
 */
template <typename T, int NDIMS>
bool is_shallow(const Polygon<T, NDIMS>& poly, double EPS = 1e-8)
{
  // max index for vertices in the full polygon
  const int n_full = poly.numVertices() - 1;
  if(n_full + 1 < 2)
  {
    return true;  // lines and points are shallow
  }

  int n = 0;  // max index for vertices in the nondegenerate polygon
  Polygon<T, NDIMS> simple_poly;
  simple_poly.addVertex(poly[0]);

  for(int i = 1; i < n_full; ++i)
  {
    if(squared_distance(simple_poly[n], poly[i]) > 1e-8)
    {
      simple_poly.addVertex(poly[i]);
      n++;
    }
  }

  if((n > 0) && (squared_distance(simple_poly[0], poly[n_full]) > 1e-8) &&
     (squared_distance(simple_poly[n - 1], poly[n_full]) > 1e-8))
  {
    simple_poly.addVertex(poly[poly.numVertices() - 1]);
    n++;
  }

  if(n + 1 < 2)
  {
    return true;  // lines and points are shallow
  }

  // Check two edges connected to P[n]P[0].
  auto compute_angle = [](const Point<T, NDIMS>& a,
                          const Point<T, NDIMS>& b,
                          const Point<T, NDIMS>& c) -> double {
    Vector<T, NDIMS> v1, v2;
    v1 = Vector<T, NDIMS>(b, a).unitVector();
    v2 = Vector<T, NDIMS>(b, c).unitVector();

    return acos(axom::utilities::clampVal(v1.dot(v2), -1.0, 1.0));
  };

  // Check endpoints
  if(compute_angle(simple_poly[n], simple_poly[0], simple_poly[1]) +
       compute_angle(simple_poly[n - 1], simple_poly[n], simple_poly[0]) >=
     0.5 * M_PI)
  {
    return false;
  }

  // Iterate over the middle vertices
  for(int i = 1; i < n; ++i)
  {
    if(compute_angle(simple_poly[i - 1], simple_poly[i], simple_poly[i + 1]) <
       0.5 * M_PI)
    {
      return false;
    }
  }

  return true;
}

}  // namespace primal
}  // namespace axom

#endif

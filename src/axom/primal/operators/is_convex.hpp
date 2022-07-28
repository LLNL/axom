// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_IS_CONVEX_HPP_
#define AXOM_PRIMAL_IS_CONVEX_HPP_

#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/Polygon.hpp"
#include "axom/primal/operators/orientation.hpp"

namespace axom
{
namespace primal
{
/*!
 * \brief Determines if a polygon defined by ordered vertices is convex
 * 
 * \param [in] poly The polygon
 * 
 * Uses dot products to detect whether vertices extend in the "convex" direction. 
 * Algorithm adapted from:
 * Y. L. Ma, W.T. Hewitt. "Point inversion and projection for NURBS curve 
 * and surface: Control polygon approach"
 * Computer Aided Geometric Design 20(2):79-99, May 2003.
 * \note Only defined in 2D
 * 
 * \return A boolean value indicating convexity
 */
template <typename T>
bool is_convex(const Polygon<T, 2>& poly)
{
  int n = poly.numVertices() - 1;
  if(n + 1 < 3) return true;  // Triangles and lines are convex

  for(int i = 1; i < n; i++)
  {
    // For each non-endpoint, check if that point and one of the endpoints
    //  are on the same side as the segment connecting the adjacent nodes
    Segment<T, 2> seg(poly[i - 1], poly[i + 1]);
    int res1 = orientation(poly[i], seg);

    // Edge case
    if(res1 == primal::ON_BOUNDARY) continue;

    // Ensure other point to check against isn't adjacent
    if(res1 == orientation(poly[(i < n / 2) ? n : 0], seg)) return false;
  }

  return true;
}

}  // namespace primal
}  // namespace axom

#endif
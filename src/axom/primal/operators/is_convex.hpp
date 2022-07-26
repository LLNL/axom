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
 * Uses algorithm adapted from [Maa 99], which uses dot products to detect
 * whether vertices extend in the "convex" direction. 
 *
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

    if(i < n / 2)
      if(res1 == orientation(poly[n], seg))
        return false;
      else if(res1 == orientation(poly[0], seg))
        return false;
  }

  return true;
}

}  // namespace primal
}  // namespace axom

#endif
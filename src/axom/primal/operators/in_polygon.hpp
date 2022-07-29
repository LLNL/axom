// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file in_polygon.hpp
 *
 * \brief Consists of methods that test whether a given query point is
 * inside a given polygon.
 *
 * Uses a ray casting algorithm
 */

#ifndef AXOM_PRIMAL_IN_POLYGON_HPP_
#define AXOM_PRIMAL_IN_POLYGON_HPP_

// Axom includes
#include "axom/config.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Polygon.hpp"

// C++ includes
#include <cmath>

namespace axom
{
namespace primal
{
/*!
 * \brief Computes the winding number for a point and a polygon
 *
 * \param [in] R The query point to test
 * \param [in] P The Polygon object to test for containment
 * \param [in[ strict If true, points on the boundary are considered exterior.
 * \param [in] EPS The tolerance level for collinearity
 * 
 * Uses an adapted ray-casting approach that counts quarter-rotation
 * of vertices around the query point. 
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
                   const bool strict = false,
                   const double EPS = 1e-8)
{
  const int nverts = P.numVertices();

  // If the query is a vertex, return a value interpreted
  //  as "inside" by evenodd or nonzero protocols
  if(R == P[0]) return !strict;

  int winding_num = 0;
  for(int i = 0; i < nverts; i++)
  {
    int j = (i == nverts - 1) ? 0 : i + 1;

    if(P[j][1] == R[1])
    {
      if(P[j][0] == R[0])
        return !strict;  // On vertex
      else if(P[i][1] == R[1] && ((P[j][0] > R[0]) == (P[i][0] < R[0])))
        return !strict;  // On horizontal edge
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
          if(axom::utilities::isNearlyEqual(det, 0.0, EPS)) return !strict;

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
          if(axom::utilities::isNearlyEqual(det, 0.0, EPS)) return !strict;

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
 * \brief Determines containment for a point in a polygon
 *
 * \param [in] query The query point to test
 * \param [in] poly The Polygon object to test for containment
 * \param [in] nonzero If false, use even/odd protocol for inclusion
 * \param [in] strict If true, points on the boundary are considered exterior.
 * \param [in] EPS The tolerance level for collinearity
 * 
 * Determines contianment using the winding number with respect to the 
 * given polygon. 
 * Uses different protocols to determine containment from the winding number.
 *   nonzero (true): If the winding number is nonzero, the point is interior.
 *   evenodd (false): If the winding number is odd, it is interior. Exterior otherwise.
 
 * \return boolean value indicating containment.
 */
template <typename T>
bool in_polygon(const Point<T, 2>& query,
                const Polygon<T, 2>& poly,
                const bool nonzero = true,
                const bool strict = false,
                const double EPS = 1e-8)
{
  if(nonzero == true) return winding_number(query, poly, strict, EPS) != 0;
  // else, use evenodd rule
  return (winding_number(query, poly, strict, EPS) % 2) == 1;
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_IN_CURVED_POLYGON_H_

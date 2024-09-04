// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_IS_CONVEX_HPP_
#define AXOM_PRIMAL_IS_CONVEX_HPP_

#include "axom/primal/geometry/Hexahedron.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/Polygon.hpp"
#include "axom/primal/operators/intersect.hpp"
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
 * \brief Determines if a hexahedron is convex
 *
 * \param [in] hex The hexahedron
 * \param [in] eps The tolerance
 *
 * A hexahedron is convex when:
 * - The faces are planar
 * - Rays sent out from an arbitrary point inside the
 *   hexahedron intersect only one face.
 *
 * \return A boolean value indicating convexity
 */
template <typename T>
AXOM_HOST_DEVICE bool is_convex(const Hexahedron<T, 3>& hex, double EPS = 1e-8)
{
  using PlaneType = Plane<T, 3>;
  using PointType = Point<T, 3>;
  using SegmentType = Segment<T, 3>;

  const int NUM_FACES = 6;

  PlaneType faces[NUM_FACES];
  faces[0] = make_plane(hex[0], hex[1], hex[2]);
  faces[1] = make_plane(hex[0], hex[1], hex[4]);
  faces[2] = make_plane(hex[0], hex[3], hex[4]);
  faces[3] = make_plane(hex[1], hex[2], hex[5]);
  faces[4] = make_plane(hex[2], hex[3], hex[6]);
  faces[5] = make_plane(hex[4], hex[5], hex[6]);

  // Check if faces are planar first
  bool isPlanar = (faces[0].getOrientation(hex[3], EPS) == primal::ON_BOUNDARY) &&
    (faces[1].getOrientation(hex[5], EPS) == primal::ON_BOUNDARY) &&
    (faces[2].getOrientation(hex[7], EPS) == primal::ON_BOUNDARY) &&
    (faces[3].getOrientation(hex[6], EPS) == primal::ON_BOUNDARY) &&
    (faces[4].getOrientation(hex[7], EPS) == primal::ON_BOUNDARY) &&
    (faces[5].getOrientation(hex[7], EPS) == primal::ON_BOUNDARY);

  if(!isPlanar)
  {
    return false;
  }

  // Get arbitrary point inside hexahedron (vertexMean in this case) and
  // face means

  // Hex center (hc)
  PointType hc = hex.vertexMean();

  PointType face_means[NUM_FACES];

  face_means[0] = PointType::midpoint(PointType::midpoint(hex[0], hex[1]),
                                      PointType::midpoint(hex[2], hex[3]));

  face_means[1] = PointType::midpoint(PointType::midpoint(hex[0], hex[1]),
                                      PointType::midpoint(hex[4], hex[5]));

  face_means[2] = PointType::midpoint(PointType::midpoint(hex[0], hex[3]),
                                      PointType::midpoint(hex[4], hex[7]));

  face_means[3] = PointType::midpoint(PointType::midpoint(hex[1], hex[2]),
                                      PointType::midpoint(hex[5], hex[6]));

  face_means[4] = PointType::midpoint(PointType::midpoint(hex[2], hex[3]),
                                      PointType::midpoint(hex[6], hex[7]));

  face_means[5] = PointType::midpoint(PointType::midpoint(hex[4], hex[5]),
                                      PointType::midpoint(hex[6], hex[7]));

  // For each segment from the hexahedron center to the face mean, verify
  // segment passes through only one face.

  //Unused lerp value
  T t;
  for(int i = 0; i < NUM_FACES; i++)
  {
    SegmentType seg(hc, face_means[i]);
    for(int j = 0; j < NUM_FACES; j++)
    {
      // Skip plane corresponding to current segment
      if(i == j)
      {
        continue;
      }

      bool intersect_res = axom::primal::intersect(faces[j], seg, t);
      if(intersect_res)
      {
        return false;
      }
    }
  }

  return true;
}

}  // namespace primal
}  // namespace axom

#endif

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
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/Ray.hpp"
#include "axom/primal/geometry/Polygon.hpp"
#include "axom/primal/operators/detail/intersect_ray_impl.hpp"

namespace axom
{
namespace primal
{

// Check if point is interior to polygon.
// Only possible in 2D
template <typename T>
bool in_polygon(const Point<T, 2>& p, const Polygon<T, 2>& poly, const double EPS = 1e-8)
{
  Ray<T, 2> the_ray(p, Vector<T, 2>({1, 0}));
  double ray_param = -1, seg_param = -1;
  int n = poly.numVertices() - 1;

  Segment<T, 2> the_seg(poly[n], poly[0]);

  // To avoid double counting vertices, let the segment be open at origin
  int num_intersects =
    (detail::intersect_ray(the_ray, the_seg, ray_param, seg_param, EPS) &&
     seg_param > EPS);

  // If the ray intersects the segment at its endpoint, then the query point
  //  is on the polygon. Interpret this as "inside." Allow for tolerance
  //  consistent with intersect_ray
  if(num_intersects == 1 && ray_param < EPS) return true;

  for(int i = 0; i < n; i++)
  {
    bool this_intersect =
      detail::intersect_ray(the_ray,
                            Segment<T, 2>(poly[i], poly[i + 1]),
                            ray_param,
                            seg_param,
                            EPS);

    if(this_intersect == true && seg_param > EPS)
    {
      if(ray_param < EPS)
        return true;
      else
        ++num_intersects;
    }
  }

  // Return true if num_intersects is odd
  return (num_intersects % 2 == 1);
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_IN_CURVED_POLYGON_H_

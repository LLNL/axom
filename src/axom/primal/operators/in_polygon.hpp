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
/*!
 * \brief Determine if query point is interior to a polygon
 *
 * \param [in] query The query point to test
 * \param [in] cpoly The CurvedPolygon object to test for containment
 * \param [in] EPS The tolerance level at which a Bezier curve is linear
 * \param [in[ strict If true, check for strict inclusion
 * 
 * Determines contianment using a ray casting appraoch built from primal 
 * primaltives. Checks intersection of ray with origin at query in the 
 * direction [1, 0] with each edge of the polygon. If the number of 
 * intersections is odd, it is in the interior.
 * 
 * \return A boolean value indicating containment.
 */
template <typename T>
bool in_polygon(const Point<T, 2>& p,
                const Polygon<T, 2>& poly,
                const double EPS = 1e-8,
                const bool strict = false)
{
  Ray<T, 2> the_ray(p, Vector<T, 2>({1, 0}));
  double ray_param = -1, seg_param = -1;
  int n = poly.numVertices() - 1;

  Segment<T, 2> the_seg(poly[n], poly[0]);

  // To avoid double counting vertices, interpret the segment as open at origin
  int num_intersects =
    (detail::intersect_ray(the_ray, the_seg, ray_param, seg_param, EPS) &&
     seg_param > EPS);

  // If the ray intersects the segment at its origin, then the query point
  //  is on the boundary of the polygon. Allow for tolerance
  //  consistent with intersect_ray
  if(num_intersects == 1 && ray_param < EPS) return (strict) ? false : true;

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
        return (strict) ? false : true;
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

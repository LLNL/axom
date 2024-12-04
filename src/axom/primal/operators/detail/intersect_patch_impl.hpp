// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file intersect_patch_impl.hpp
 *
 * This file provides helper functions for testing the intersection
 * of rays and Bezier patches
 */

#ifndef AXOM_PRIMAL_INTERSECT_PATCH_IMPL_HPP_
#define AXOM_PRIMAL_INTERSECT_PATCH_IMPL_HPP_

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Polygon.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/BezierPatch.hpp"

#include "axom/primal/operators/intersect.hpp"
#include "axom/primal/operators/in_polygon.hpp"
#include "axom/primal/operators/detail/intersect_impl.hpp"
#include "axom/primal/operators/detail/intersect_ray_impl.hpp"

#include <vector>

namespace axom
{
namespace primal
{
namespace detail
{
//---------------------------- FUNCTION DECLARATIONS ---------------------------

template <typename T>
bool intersect_line_patch(const Line<T, 3> &line,
                          const BezierPatch<T, 3> &patch,
                          axom::Array<T> &tp,
                          axom::Array<T> &up,
                          axom::Array<T> &vp,
                          int order_u,
                          int order_v,
                          double u_offset,
                          double u_scale,
                          double v_offset,
                          double v_scale,
                          double sq_tol,
                          double EPS,
                          bool isRay);

//------------------------------ IMPLEMENTATIONS ------------------------------

template <typename T>
bool intersect_line_patch(const Line<T, 3> &line,
                          const BezierPatch<T, 3> &patch,
                          axom::Array<T> &tp,
                          axom::Array<T> &up,
                          axom::Array<T> &vp,
                          int order_u,
                          int order_v,
                          double u_offset,
                          double u_scale,
                          double v_offset,
                          double v_scale,
                          double sq_tol,
                          double EPS,
                          bool isRay)
{
  using BPatch = BezierPatch<T, 3>;

  // Check bounding box to short-circuit the intersection
  //  Need to expand the box a bit so that intersections near subdivision boundaries
  //  are accurately recorded
  Point<T, 3> ip;
  if(!intersect(line, patch.boundingBox().scale(1.5), ip))
  {
    return false;
  }

  bool foundIntersection = false;
  if(patch.isBilinear(sq_tol))
  {
    // Store candidate intersection points
    axom::Array<T> tc, uc, vc;

    foundIntersection =
      detail::intersect_line_bilinear_patch(line,
                                            patch(0, 0),
                                            patch(order_u, 0),
                                            patch(order_u, order_v),
                                            patch(0, order_v),
                                            tc,
                                            uc,
                                            vc,
                                            EPS,
                                            isRay);

    if(!foundIntersection)
    {
      return false;
    }

    foundIntersection = false;
    for(int i = 0; i < tc.size(); ++i)
    {
      const T t0 = tc[i];
      const T u0 = uc[i];
      const T v0 = vc[i];

      // Use EPS to record points near the boundary of the bilinear approximation
      if(u0 >= -EPS / u_scale && u0 <= 1.0 + EPS / u_scale &&
         v0 >= -EPS / v_scale && v0 <= 1.0 + EPS / v_scale)
      {
        if(t0 >= -EPS || !isRay)
        {
          up.push_back(u_offset + u0 * u_scale);
          vp.push_back(v_offset + v0 * v_scale);
          tp.push_back(t0);

          foundIntersection = true;
        }
      }
    }
  }
  else
  {
    constexpr double splitVal = 0.5;
    constexpr double scaleFac = 0.5;

    BPatch p1(order_u, order_v), p2(order_u, order_v), p3(order_u, order_v),
      p4(order_u, order_v);

    patch.split(splitVal, splitVal, p1, p2, p3, p4);
    u_scale *= scaleFac;
    v_scale *= scaleFac;

    // Note: we want to find all intersections, so don't short-circuit
    if(intersect_line_patch(line,
                            p1,
                            tp,
                            up,
                            vp,
                            order_u,
                            order_v,
                            u_offset,
                            u_scale,
                            v_offset,
                            v_scale,
                            sq_tol,
                            EPS,
                            isRay))
    {
      foundIntersection = true;
    }
    if(intersect_line_patch(line,
                            p2,
                            tp,
                            up,
                            vp,
                            order_u,
                            order_v,
                            u_offset + u_scale,
                            u_scale,
                            v_offset,
                            v_scale,
                            sq_tol,
                            EPS,
                            isRay))
    {
      foundIntersection = true;
    }
    if(intersect_line_patch(line,
                            p3,
                            tp,
                            up,
                            vp,
                            order_u,
                            order_v,
                            u_offset,
                            u_scale,
                            v_offset + v_scale,
                            v_scale,
                            sq_tol,
                            EPS,
                            isRay))
    {
      foundIntersection = true;
    }
    if(intersect_line_patch(line,
                            p4,
                            tp,
                            up,
                            vp,
                            order_u,
                            order_v,
                            u_offset + u_scale,
                            u_scale,
                            v_offset + v_scale,
                            v_scale,
                            sq_tol,
                            EPS,
                            isRay))
    {
      foundIntersection = true;
    }
  }

  return foundIntersection;
}

}  // end namespace detail
}  // end namespace primal
}  // end namespace axom

#endif  // AXOM_PRIMAL_INTERSECT_PATCH_IMPL_HPP_

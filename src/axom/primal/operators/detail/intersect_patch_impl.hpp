// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file intersect_bezier_impl.hpp
 *
 * This file provides helper functions for testing the intersection
 * of Bezier curves
 */

#ifndef AXOM_PRIMAL_INTERSECT_PATCH_IMPL_HPP_
#define AXOM_PRIMAL_INTERSECT_PATCH_IMPL_HPP_

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Polygon.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/BezierPatch.hpp"

#include "axom/primal/operators/intersect.hpp"
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
bool intersect_ray_patch(const BezierPatch<T, 3> &p,
                         const Ray<T, 3> &r,
                         std::vector<T> &up,
                         std::vector<T> &vp,
                         std::vector<T> &rp,
                         double sq_tol,
                         int order_u,
                         int order_v,
                         double u_offset,
                         double u_scale,
                         double v_offset,
                         double v_scale);


//------------------------------ IMPLEMENTATIONS ------------------------------

template <typename T>
bool intersect_ray_patch(const BezierPatch<T, 3> &p,
                         const Ray<T, 3> &r,
                         std::vector<T> &up,
                         std::vector<T> &vp,
                         std::vector<T> &rp,
                         double sq_tol,
                         int order_u,
                         int order_v,
                         double u_offset,
                         double u_scale,
                         double v_offset,
                         double v_scale)
{
  using BPatch = BezierPatch<T, 3>;

  // Check bounding box to short-circuit the intersection
  T r0, s0, c0;
  Point<T, 2> ip;
  if(!intersect(r, p.boundingBox(), ip))
  {
    return false;
  }

  bool foundIntersection = false;
  
  if(p.isPlanar(sq_tol))
  {
    // Add this later
    Polygon<T, 3>(axom::Array<Point<T, 3>>(
        {p(0, 0), p(order_u, 0), p(order_u, order_v), p(0, order_v)}));

    
    if(intersect(r, seg, r0, s0) && s0 <= 1.0 - 1e-8)
    {
      rp.push_back(r0);
      cp.push_back(c_offset + c_scale * s0);
      foundIntersection = true;
    }
  }
  else
  {
    constexpr double splitVal = 0.5;
    constexpr double scaleFac = 0.5;

    BPatch p1(order_u, order_v), p2(order_u, order_v), p3(order_u, order_v),
      p4(order_u, order_v);

    p.split(splitVal, splitVal, p1, p2, p3, p4);
    nevals += 1;
    u_scale *= scaleFac;
    v_scale *= scaleFac;

    // Note: we want to find all intersections, so don't short-circuit
    if(intersect_ray_patch(p1,
                           r,
                           up,
                           vp,
                           rp,
                           sq_tol,
                           order_u,
                           order_v,
                           u_offset,
                           u_scale,
                           v_offset,
                           v_scale,
                           nevals))
    {
      foundIntersection = true;
    }
    if(intersect_ray_patch(p2,
                           r,
                           up,
                           vp,
                           rp,
                           sq_tol,
                           order_u,
                           order_v,
                           u_offset + u_scale,
                           u_scale,
                           v_offset,
                           v_scale,
                           nevals))
    {
      foundIntersection = true;
    }
    if(intersect_ray_patch(p3,
                           r,
                           up,
                           vp,
                           rp,
                           sq_tol,
                           order_u,
                           order_v,
                           u_offset,
                           u_scale,
                           v_offset + v_scale,
                           v_scale,
                           nevals))
    {
      foundIntersection = true;
    }
    if(intersect_ray_patch(p4,
                           r,
                           up,
                           vp,
                           rp,
                           sq_tol,
                           order_u,
                           order_v,
                           u_offset + u_scale,
                           u_scale,
                           v_offset + v_scale,
                           v_scale,
                           nevals))
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

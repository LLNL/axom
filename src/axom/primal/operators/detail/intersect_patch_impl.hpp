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
bool intersect_ray_patch_approximate(const BezierPatch<T, 3> &p,
                                     const Ray<T, 3> &r,
                                     std::vector<T> &up,
                                     std::vector<T> &vp,
                                     double sq_tol,
                                     int order_u,
                                     int order_v,
                                     double u_offset,
                                     double u_scale,
                                     double v_offset,
                                     double v_scale);

//------------------------------ IMPLEMENTATIONS ------------------------------

template <typename T>
bool intersect_ray_patch_approximate(const BezierPatch<T, 3> &p,
                                     const Ray<T, 3> &r,
                                     std::vector<T> &up,
                                     std::vector<T> &vp,
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
  Point<T, 3> ip;
  if(!intersect(r, p.boundingBox().scale(1.1), ip))
  {
    return false;
  }

  bool foundIntersection = false;

  if(p.isPlanar(sq_tol))
  {
    // For the purposes of computing the winding number,
    //  we don't need the exact parameters of the intersection.
    // If we intersect the planar patch,
    //  just take the center of the remaining parameters

    // Project the corners onto the plane normal to the ray
    Vector<T, 3> ray_dir = r.direction().unitVector();
    Vector<T, 3> n1, n2;

    // Define vectors perpendicular to the ray
    if(!axom::utilities::isNearlyEqual(ray_dir[1], ray_dir[0]))
      n1 = Vector<T, 3>({-ray_dir[1], ray_dir[0], 0.0}).unitVector();
    else
      n1 = Vector<T, 3>({-ray_dir[2], 0.0, ray_dir[0]}).unitVector();
    n2 = Vector<T, 3>::cross_product(ray_dir, n1);

    // Find the constants for the projection
    Point<T, 3> ray_origin = r.origin();
    double e1 =
      -(n1[0] * ray_origin[0] + n1[1] * ray_origin[1] + n1[2] * ray_origin[2]);
    double e2 =
      -(n2[0] * ray_origin[0] + n2[1] * ray_origin[1] + n2[2] * ray_origin[2]);

    // Define the 2D polygon
    Polygon<T, 2> poly;
    poly.addVertex(Point<T, 2>(
      {e1 + n1[0] * p(0, 0)[0] + n1[1] * p(0, 0)[1] + n1[2] * p(0, 0)[2],
       e2 + n2[0] * p(0, 0)[0] + n2[1] * p(0, 0)[1] + n2[2] * p(0, 0)[2]}));
    poly.addVertex(
      Point<T, 2>({e1 + n1[0] * p(order_u, 0)[0] + n1[1] * p(order_u, 0)[1] +
                     n1[2] * p(order_u, 0)[2],
                   e2 + n2[0] * p(order_u, 0)[0] + n2[1] * p(order_u, 0)[1] +
                     n2[2] * p(order_u, 0)[2]}));
    poly.addVertex( Point<T, 2>(
      {e1 + n1[0] * p(order_u, order_v)[0] + n1[1] * p(order_u, order_v)[1] +
         n1[2] * p(order_u, order_v)[2],
       e2 + n2[0] * p(order_u, order_v)[0] + n2[1] * p(order_u, order_v)[1] +
         n2[2] * p(order_u, order_v)[2]}));
    poly.addVertex(
      Point<T, 2>({e1 + n1[0] * p(0, order_v)[0] + n1[1] * p(0, order_v)[1] +
                     n1[2] * p(0, order_v)[2],
                   e2 + n2[0] * p(0, order_v)[0] + n2[1] * p(0, order_v)[1] +
                     n2[2] * p(0, order_v)[2]}));

    if( p.isRational() )
    {
        poly[0][0] *= p.getWeight(0, 0);
        poly[0][1] *= p.getWeight(0, 0);

        poly[1][0] *= p.getWeight(order_u, 0);
        poly[1][1] *= p.getWeight(order_u, 0);

        poly[2][0] *= p.getWeight(order_u, order_v);
        poly[2][1] *= p.getWeight(order_u, order_v);

        poly[3][0] *= p.getWeight(0, order_v);
        poly[3][1] *= p.getWeight(0, order_v);
    }

    if(in_polygon(Point<T, 2>({0.0, 0.0}), poly))
    {
      up.push_back(u_offset + 0.5 * u_scale);
      vp.push_back(v_offset + 0.5 * v_scale);
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
    u_scale *= scaleFac;
    v_scale *= scaleFac;

    // Note: we want to find all intersections, so don't short-circuit
    if(intersect_ray_patch_approximate(p1,
                                       r,
                                       up,
                                       vp,
                                       sq_tol,
                                       order_u,
                                       order_v,
                                       u_offset,
                                       u_scale,
                                       v_offset,
                                       v_scale))
    {
      foundIntersection = true;
    }
    if(intersect_ray_patch_approximate(p2,
                                       r,
                                       up,
                                       vp,
                                       sq_tol,
                                       order_u,
                                       order_v,
                                       u_offset + u_scale,
                                       u_scale,
                                       v_offset,
                                       v_scale))
    {
      foundIntersection = true;
    }
    if(intersect_ray_patch_approximate(p3,
                                       r,
                                       up,
                                       vp,
                                       sq_tol,
                                       order_u,
                                       order_v,
                                       u_offset,
                                       u_scale,
                                       v_offset + v_scale,
                                       v_scale))
    {
      foundIntersection = true;
    }
    if(intersect_ray_patch_approximate(p4,
                                       r,
                                       up,
                                       vp,
                                       sq_tol,
                                       order_u,
                                       order_v,
                                       u_offset + u_scale,
                                       u_scale,
                                       v_offset + v_scale,
                                       v_scale))
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

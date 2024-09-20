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
bool intersect_line_patch(const BezierPatch<T, 3> &patch,
                          const Line<T, 3> &line,
                          std::vector<T> &up,
                          std::vector<T> &vp,
                          std::vector<T> &tp,
                          double sq_tol,
                          int order_u,
                          int order_v,
                          double u_offset,
                          double u_scale,
                          double v_offset,
                          double v_scale);

//------------------------------ IMPLEMENTATIONS ------------------------------

template <typename T>
bool intersect_line_patch(const BezierPatch<T, 3> &patch,
                          const Line<T, 3> &line,
                          std::vector<T> &up,
                          std::vector<T> &vp,
                          std::vector<T> &tp,
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
  if(!intersect(line, patch.boundingBox().scale(1.1), ip))
  {
    return false;
  }

  bool foundIntersection = false;
  if(patch.isBilinear(1e-4))
  {
    std::vector<T> u, v, t;
    foundIntersection =
      detail::intersect_bilinear_patch_line(patch(0, 0),
                                            patch(order_u, 0),
                                            patch(order_u, order_v),
                                            patch(0, order_v),
                                            line,
                                            u,
                                            v,
                                            t);

    for(size_t i = 0; i < u.size(); ++i)
    {
      if(axom::utilities::isNearlyEqual(u[i], 1.0, 1e-8) && u_scale != 1.0)
      {
        continue;
      }
      if(axom::utilities::isNearlyEqual(v[i], 1.0, 1e-8) && v_scale != 1.0)
      {
        continue;
      }

      up.push_back(u_offset + u[i] * u_scale);
      vp.push_back(v_offset + v[i] * v_scale);
      tp.push_back(t[i]);
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
    if(intersect_line_patch(p1,
                           line,
                           up,
                           vp,
                           tp,
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
    if(intersect_line_patch(p2,
                           line,
                           up,
                           vp,
                           tp,
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
    if(intersect_line_patch(p3,
                           line,
                           up,
                           vp,
                           tp,
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
    if(intersect_line_patch(p4,
                           line,
                           up,
                           vp,
                           tp,
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

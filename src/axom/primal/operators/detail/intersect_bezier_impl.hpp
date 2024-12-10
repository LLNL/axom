// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file intersect_bezier_impl.hpp
 *
 * This file provides helper functions for testing the intersection
 * of Bezier curves with other Bezier curves and other geometric objects
 */

#ifndef AXOM_PRIMAL_INTERSECT_BEZIER_IMPL_HPP_
#define AXOM_PRIMAL_INTERSECT_BEZIER_IMPL_HPP_

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"

#include "axom/primal/operators/intersect.hpp"
#include "axom/primal/operators/detail/intersect_impl.hpp"

#include <vector>

namespace axom
{
namespace primal
{
namespace detail
{
//---------------------------- FUNCTION DECLARATIONS ---------------------------

/*!
 * \brief Recursive function to find the intersections between two Bezier curves
 *
 * \param [in] c1 The first input curve
 * \param [in] c2 The second input curves
 * \param [out] sp Parametric coordinates of intersections in \a c1
 * \param [out] tp Parametric coordinates of intersections in \a c2
 * \param [in] sq_tol The squared tolerance parameter for determining if a
 * Bezier curve is linear
 * \param [in] order1 The order of \a c1
 * \param [in] order2 The order of \a c2
 * \param s_offset The offset in parameter space for \a c1
 * \param s_scale The scale in parameter space for \a c1
 * \param t_offset The offset in parameter space for \a c2
 * \param t_scale The scale in parameter space for \a c2
 *
 * Bezier curves can only intersect when their bounding boxes intersect.
 * The base case of the recursion is when we can approximate the curves as
 * line segments, where we directly find their intersections. Otherwise,
 * check for intersections recursively after bisecting one of the curves.
 *
 * \note A BezierCurve is parametrized in [0,1). The scale and offset parameters
 * are used to track the local curve parameters during subdivisions
 *
 * \return True if the two curves intersect, False otherwise
 * \sa intersect_bezier
 */
template <typename T>
bool intersect_bezier_curves(const BezierCurve<T, 2> &c1,
                             const BezierCurve<T, 2> &c2,
                             axom::Array<T> &sp,
                             axom::Array<T> &tp,
                             double sq_tol,
                             int order1,
                             int order2,
                             double s_offset,
                             double s_scale,
                             double t_offset,
                             double t_scale);

/*!
 * \brief Tests intersection of two line segments defined by
 * their end points (a,b) and (c,d)
 *
 * \param [in] a,d,c,b the endpoints of the segments
 * \param [out] The parametrized s and t values at which intersection occurs
 * Range of output values for \a s and \a t is [0,1).
 *
 * \return True, if the two line segments intersect, false otherwise.
 * When the two line segments intersect, their intersection parameters
 * are stored in output parameters \a s and \a t, respectively, such that,
 * if \a p is the intersection point,
 * \f$    p = a + s * (b-a) = c + t * (d-c).  \f$
 *
 * \note We consider both line segments to have half-open intervals.
 * As such, the we do not consider the lines to intersect if they do so
 * at the endpoints \a b or \d, respectively.
 *
 * \note This function assumes the all intersections have multiplicity
 * one, i.e. there are no points at which the curves and their derivatives
 * both intersect. Thus, the function does not find tangencies.
 *
 * \note This function assumes two dimensional curves in a plane.
 *
 * \note This function does not properly handle collinear lines
 */
template <typename T>
bool intersect_2d_linear(const Point<T, 2> &a,
                         const Point<T, 2> &b,
                         const Point<T, 2> &c,
                         const Point<T, 2> &d,
                         T &s,
                         T &t);

/*!
 * \brief Recursive function to find intersections between a ray and a Bezier curve
 *
 * \param [in] c The input curve
 * \param [in] r The input ray
 * \param [out] cp Parametric coordinates of intersections in \a c [0, 1)
 * \param [out] rp Parametric coordinates of intersections in \a r [0, inf)
 * \param [in] sq_tol The squared tolerance parameter for distances in physical space
 * \param [in] EPS The tolerance parameter for distances in parameter space
 * \param [in] order The order of \a c
 * \param s_offset The offset in parameter space for \a c
 * \param s_scale The scale in parameter space for \a c
 *
 * A ray can only intersect a Bezier curve if it intersects its bounding box
 * The base case of the recursion is when we can approximate the curves with
 * line segments, where we directly find their intersection with the ray. Otherwise,
 * check for intersections recursively after bisecting the curve.
 *
 * \note A BezierCurve is parametrized in [0,1). The scale and offset parameters
 * are used to track the local curve parameters during subdivisions
 *
 * \note This function assumes the all intersections have multiplicity
 * one, i.e. there are no points at which the curves and their derivatives
 * both intersect. Thus, the function does not find tangencies.
 * 
 * \return True if the two curves intersect, False otherwise
 * \sa intersect_bezier
 */
template <typename T>
bool intersect_ray_bezier(const Ray<T, 2> &r,
                          const BezierCurve<T, 2> &c,
                          axom::Array<T> &rp,
                          axom::Array<T> &cp,
                          double sq_tol,
                          double EPS,
                          int order,
                          double c_offset,
                          double c_scale);
//------------------------------ IMPLEMENTATIONS ------------------------------

template <typename T>
bool intersect_bezier_curves(const BezierCurve<T, 2> &c1,
                             const BezierCurve<T, 2> &c2,
                             axom::Array<T> &sp,
                             axom::Array<T> &tp,
                             double sq_tol,
                             int order1,
                             int order2,
                             double s_offset,
                             double s_scale,
                             double t_offset,
                             double t_scale)
{
  using BCurve = BezierCurve<T, 2>;

  // Check bounding boxes to short-circuit the intersection
  if(!intersect(c1.boundingBox(), c2.boundingBox()))
  {
    return false;
  }

  bool foundIntersection = false;

  if(c1.isLinear(sq_tol) && c2.isLinear(sq_tol))
  {
    T s, t;
    if(intersect_2d_linear(c1[0], c1[order1], c2[0], c2[order2], s, t))
    {
      sp.push_back(s_offset + s_scale * s);
      tp.push_back(t_offset + t_scale * t);
      foundIntersection = true;
    }
  }
  else
  {
    constexpr double splitVal = 0.5;
    constexpr double scaleFac = 0.5;

    BCurve c3(order1);
    BCurve c4(order1);
    c1.split(splitVal, c3, c4);

    s_scale *= scaleFac;

    // Note: we want to find all intersections, so don't short-circuit
    if(intersect_bezier_curves(c2,
                               c3,
                               tp,
                               sp,
                               sq_tol,
                               order2,
                               order1,
                               t_offset,
                               t_scale,
                               s_offset,
                               s_scale))
    {
      foundIntersection = true;
    }
    if(intersect_bezier_curves(c2,
                               c4,
                               tp,
                               sp,
                               sq_tol,
                               order2,
                               order1,
                               t_offset,
                               t_scale,
                               s_offset + s_scale,
                               s_scale))
    {
      foundIntersection = true;
    }
  }

  return foundIntersection;
}

template <typename T>
bool intersect_2d_linear(const Point<T, 2> &a,
                         const Point<T, 2> &b,
                         const Point<T, 2> &c,
                         const Point<T, 2> &d,
                         T &s,
                         T &t)
{
  // Implementation inspired by Section 5.1.9.1 of
  // C. Ericson's Real-Time Collision Detection book

  // Note: Uses exact floating point comparisons since the subdivision algorithm
  // provides both sides of the line segments for interior curve points.

  // compute signed areas of endpoints of segment (c,d) w.r.t. segment (a,b)
  const auto area1 = twoDcross(a, b, c);
  const auto area2 = twoDcross(a, b, d);

  // early return if both have same orientation, or if d is collinear w/ (a,b)
  if(area2 == 0. || (area1 * area2) > 0.)
  {
    return false;
  }

  // compute signed areas of endpoints of segment (a,b) w.r.t. segment (c,d)
  const auto area3 = twoDcross(c, d, a);
  const auto area4 = area3 + area1 - area2;  // equivalent to twoDcross(c,d,b)

  // early return if both have same orientation, or if b is collinear w/ (c,d)
  if(area4 == 0. || (area3 * area4) > 0.)
  {
    return false;
  }

  // Compute intersection parameters using linear interpolation
  // Divisions are safe due to early return conditions
  s = area3 / (area3 - area4);
  t = area1 / (area1 - area2);

  return true;
}

template <typename T>
bool intersect_ray_bezier(const Ray<T, 2> &r,
                          const BezierCurve<T, 2> &c,
                          axom::Array<T> &rp,
                          axom::Array<T> &cp,
                          double sq_tol,
                          double EPS,
                          int order,
                          double c_offset,
                          double c_scale)
{
  using BCurve = BezierCurve<T, 2>;

  // Check bounding box to short-circuit the intersection
  T r0, s0;
  Point<T, 2> ip;
  constexpr T factor = 1e-8;

  // Need to expand the bounding box, since this ray-bb intersection routine
  //  only parameterizes the ray on (0, inf)
  if(!intersect(r, c.boundingBox().expand(factor), ip))
  {
    return false;
  }

  bool foundIntersection = false;

  // For the base case, represent the Bezier curve as a line segment
  if(c.isLinear(sq_tol))
  {
    Segment<T, 2> seg(c[0], c[order]);

    // Need to check intersection with zero tolerance
    //  to handle cases where `intersect` treats the ray as collinear
    if(intersect(r, seg, r0, s0, EPS) && s0 > 0.0 - EPS && s0 < 1.0 - EPS)
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

    BCurve c1(order);
    BCurve c2(order);
    c.split(splitVal, c1, c2);
    c_scale *= scaleFac;

    // Note: we want to find all intersections, so don't short-circuit
    if(intersect_ray_bezier(r, c1, rp, cp, sq_tol, EPS, order, c_offset, c_scale))
    {
      foundIntersection = true;
    }
    if(intersect_ray_bezier(r, c2, rp, cp, sq_tol, EPS, order, c_offset + c_scale, c_scale))
    {
      foundIntersection = true;
    }
  }

  return foundIntersection;
}

}  // end namespace detail
}  // end namespace primal
}  // end namespace axom

#endif  // AXOM_PRIMAL_INTERSECT_BEZIER_IMPL_HPP_

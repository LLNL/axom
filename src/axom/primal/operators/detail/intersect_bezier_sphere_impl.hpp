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

#ifndef AXOM_PRIMAL_INTERSECT_BEZIER_SPHERE_IMPL_HPP_
#define AXOM_PRIMAL_INTERSECT_BEZIER_SPHERE_IMPL_HPP_

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"

#include <vector>
#include <math.h>

namespace axom
{
namespace primal
{
namespace detail
{
//---------------------------- FUNCTION DECLARATIONS ---------------------------

/*!
 * \brief Recursive function to find intersections between a circle and a Bezier curve
 *
 * \param [in] circle The input circle (2D primal::Sphere)
 * \param [in] curve The input curve
 * \param [out] circle_p Parametric coordinates of intersections of the circle [0, 2pi)
 * \param [out] curve_p Parametric coordinates of intersections in \a r [0, 1)
 * \param [in] sq_tol The squared tolerance parameter for distances in physical space
 * \param [in] EPS The tolerance parameter for distances in parameter space
 * \param [in] order The order of \a c
 * \param s_offset The offset in parameter space for \a c
 * \param s_scale The scale in parameter space for \a c
 *
 * A circle can only intersect a Bezier curve if it intersects its bounding box
 * The base case of the recursion is when we can approximate the curves with
 * line segments, where we directly find their intersection with the circle. Otherwise,
 * check for intersections recursively after bisecting the curve.
 *
 * \note A BezierCurve is parametrized in [0,1). The scale and offset parameters
 * are used to track the local curve parameters during subdivisions
 *
 * \note This function assumes the all intersections have multiplicity
 * one, i.e. there are no points at which the curves and their derivatives
 * both intersect. Thus, the function does not find tangencies.
 * 
 * \return True if the circle and the curve intersect, False otherwise
 */
template <typename T>
bool intersect_circle_bezier(const Sphere<T, 2> &circle,
                             const BezierCurve<T, 2> &curve,
                             axom::Array<T> &circle_p,
                             axom::Array<T> &curve_p,
                             double sq_tol,
                             double EPS,
                             int order,
                             double c_offset,
                             double c_scale);

template <typename T>
bool intersect_2d_circle_line(const Sphere<T, 2> &circ,
                              const Point<T, 2> &a,
                              const Point<T, 2> &b,
                              T &c1,
                              T &c2,
                              T &t1,
                              T &t2);

//------------------------------ IMPLEMENTATIONS --------------------------------

template <typename T>
bool intersect_circle_bezier(const Sphere<T, 2> &circle,
                             const BezierCurve<T, 2> &curve,
                             axom::Array<T> &circle_p,
                             axom::Array<T> &curve_p,
                             double sq_tol,
                             double EPS,
                             int order,
                             double c_offset,
                             double c_scale)
{
  using BCurve = BezierCurve<T, 2>;

  // Other function put here to avoid circular dependency
  {
    primal::BoundingBox<T, 2> bb = curve.boundingBox().scale(1.1);
    T dx = axom::utilities::clampVal(circle.getCenter()[0],
                                     bb.getMin()[0],
                                     bb.getMax()[0]);
    T dy = axom::utilities::clampVal(circle.getCenter()[1],
                                     bb.getMin()[1],
                                     bb.getMax()[1]);
    if(circle.computeSignedDistance(primal::Point<T, 2>({dx, dy})) > 0)
    {
      return false;
    }
  }

  bool foundIntersection = false;

  if(curve.isLinear(sq_tol))
  {
    T c1, c2, t1, t2;
    if(intersect_2d_circle_line(circle, curve[0], curve[order], c1, c2, t1, t2))
    {
      if(t1 >= 0.0 && t1 < 1.0)
      {
        circle_p.push_back(c1);
        curve_p.push_back(c_offset + c_scale * t1);
        foundIntersection = true;
      }

      if(t2 >= 0.0 && t2 < 1.0)
      {
        intersect_2d_circle_line(circle, curve[0], curve[order], c1, c2, t1, t2);
        circle_p.push_back(c2);
        curve_p.push_back(c_offset + c_scale * t2);
        foundIntersection = true;
      }
    }
  }
  else
  {
    constexpr double splitVal = 0.5;
    constexpr double scaleFac = 0.5;

    BCurve c1(order);
    BCurve c2(order);
    curve.split(splitVal, c1, c2);
    c_scale *= scaleFac;

    // Note: we want to find all intersections, so don't short-circuit
    if(intersect_circle_bezier(circle,
                               c1,
                               circle_p,
                               curve_p,
                               sq_tol,
                               EPS,
                               order,
                               c_offset,
                               c_scale))
    {
      foundIntersection = true;
    }
    if(intersect_circle_bezier(circle,
                               c2,
                               circle_p,
                               curve_p,
                               sq_tol,
                               EPS,
                               order,
                               c_offset + c_scale,
                               c_scale))
    {
      foundIntersection = true;
    }
  }

  return foundIntersection;
}

// Find intersections between a circle and an *infinite* line passing through a and b,
//  using the algorithm from Wolfram MathWorld. Record intersections
//  for the circle in the domain [0, 2pi).
// Record tangent line as *not* intersecting.
template <typename T>
bool intersect_2d_circle_line(const Sphere<T, 2> &circ,
                              const Point<T, 2> &a,
                              const Point<T, 2> &b,
                              T &c1,
                              T &c2,
                              T &t1,
                              T &t2)
{
  T dx = b[0] - a[0];
  T dy = b[1] - a[1];
  T dr = std::sqrt(dx * dx + dy * dy);

  T D = (a[0] - circ.getCenter()[0]) * (b[1] - circ.getCenter()[1]) -
    (b[0] - circ.getCenter()[0]) * (a[1] - circ.getCenter()[1]);

  T disc = circ.getRadius() * circ.getRadius() * dr * dr - D * D;

  // Treat tangencies as *not* intersecting
  if(disc <= 0.0)
  {
    return false;
  }

  disc = std::sqrt(disc);

  T x1 = (D * dy + (dy < 0 ? -1 : 1) * dx * disc) / (dr * dr);
  T x2 = (D * dy - (dy < 0 ? -1 : 1) * dx * disc) / (dr * dr);

  T y1 = (-D * dx + std::abs(dy) * disc) / (dr * dr);
  T y2 = (-D * dx - std::abs(dy) * disc) / (dr * dr);

  // Find the parameter values for the circle
  c1 = std::atan2(y1, x1);
  c1 = (c1 < 0.0) ? c1 + 2.0 * M_PI : c1;

  c2 = std::atan2(y2, x2);
  c2 = (c2 < 0.0) ? c2 + 2.0 * M_PI : c2;

  // Find the parameter values for the line
  if(std::abs(dx) > std::abs(dy))
  {
    t1 = (x1 - a[0] + circ.getCenter()[0]) / dx;
    t2 = (x2 - a[0] + circ.getCenter()[0]) / dx;
  }
  else
  {
    t1 = (y1 - a[1] + circ.getCenter()[1]) / dy;
    t2 = (y2 - a[1] + circ.getCenter()[1]) / dy;
  }

  return true;
}

}  // end namespace detail
}  // end namespace primal
}  // end namespace axom

#endif  // AXOM_PRIMAL_INTERSECT_BEZIER_SPHERE_IMPL_HPP_

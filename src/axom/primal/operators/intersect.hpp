// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file intersect.hpp
 *
 * \brief Consists of functions to test intersection among geometric primitives.
 */

#ifndef AXOM_PRIMAL_INTERSECT_HPP_
#define AXOM_PRIMAL_INTERSECT_HPP_

#include "axom/config.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/utilities/Utilities.hpp"

#include "axom/primal/constants.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/OrientedBoundingBox.hpp"
#include "axom/primal/geometry/Plane.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Polygon.hpp"
#include "axom/primal/geometry/Ray.hpp"
#include "axom/primal/geometry/Line.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/Sphere.hpp"
#include "axom/primal/geometry/Tetrahedron.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/geometry/BezierPatch.hpp"
#include "axom/primal/geometry/NURBSCurve.hpp"
#include "axom/primal/geometry/NURBSPatch.hpp"

#include "axom/primal/operators/detail/intersect_impl.hpp"
#include "axom/primal/operators/detail/intersect_ray_impl.hpp"
#include "axom/primal/operators/detail/intersect_bounding_box_impl.hpp"
#include "axom/primal/operators/detail/intersect_bezier_impl.hpp"
#include "axom/primal/operators/detail/intersect_patch_impl.hpp"

namespace axom
{
namespace primal
{
/// \name Triangle Intersection Routines
/// \accelerated
/// @{

/*!
 * \brief Tests if 3D Triangles t1 and t2 intersect.
 * 
 * \param [in] t1 The first triangle
 * \param [in] t2 The second triangle
 * \param [in] includeBoundary Indicates if boundaries should be considered
 * when detecting intersections (default: false)
 * \param [in] EPS Tolerance for determining intersections (default: 1E-8)
 * \return status true iff t1 intersects with t2, otherwise, false.
 *
 * If parameter \a includeBoundary is false (default), this function will
 * return true if the interior of t1 intersects the interior of t2.  To include
 * triangle boundaries in intersections, specify \a includeBoundary as true.
 */
template <typename T>
AXOM_HOST_DEVICE bool intersect(const Triangle<T, 3>& t1,
                                const Triangle<T, 3>& t2,
                                bool includeBoundary = false,
                                double EPS = 1E-08)
{
  return detail::intersect_tri3D_tri3D<T>(t1, t2, includeBoundary, EPS);
}

/*!
 * \brief Tests if 2D Triangles t1 and t2 intersect.
 * \param [in] t1 The first triangle
 * \param [in] t2 The second triangle
 * \param [in] includeBoundary Indicates if boundaries should be considered
 * when detecting intersections (default: false)
 * \param [in] EPS Tolerance for determining intersections (default: 1E-8)
 * \return status true iff t1 intersects with t2, otherwise, false.
 *
 * If parameter \a includeBoundary is false (default), this function will
 * return true if the interior of t1 intersects the interior of t2.  To include
 * triangle boundaries in intersections, specify \a includeBoundary as true.
 */
template <typename T>
bool intersect(const Triangle<T, 2>& t1,
               const Triangle<T, 2>& t2,
               bool includeBoundary = false,
               double EPS = 1E-08)
{
  return detail::intersect_tri2D_tri2D<T>(t1, t2, includeBoundary, EPS);
}

/*!
 * \brief Determines if a triangle and a bounding box intersect
 * \param [in] tri user-supplied triangle (with three vertices).
 * \param [in] bb user-supplied axis aligned bounding box.
 * \return true iff tri intersects with bb, otherwise, false.
 */
template <typename T>
bool intersect(const Triangle<T, 3>& tri, const BoundingBox<T, 3>& bb)
{
  return detail::intersect_tri_bbox(tri, bb);
}

/*!
 * \brief Determines if a 3D triangle intersects a 3D ray.
 * \param [in] tri A 3D triangle
 * \param [in] ray A 3D ray
 * \return true iff tri intersects with ray, otherwise, false.
 */
template <typename T>
bool intersect(const Triangle<T, 3>& tri, const Ray<T, 3>& ray)
{
  T t = T();
  Point<double, 3> p;
  return detail::intersect_tri_ray(tri, ray, t, p);
}

/*!
 * \brief Determines if a 3D triangle intersects a 3D ray.
 * \param [in] tri A 3D triangle
 * \param [in] ray A 3D ray
 * \param [out] t Intersection point of tri and R, w.r.t. parametrization of R
 * \note If there is an intersection, the intersection point is:  R.at(t)
 * \return true iff tri intersects with ray, otherwise, false.
 */
template <typename T>
bool intersect(const Triangle<T, 3>& tri, const Ray<T, 3>& ray, T& t)
{
  Point<double, 3> p;
  return detail::intersect_tri_ray(tri, ray, t, p);
}

/*!
 * \brief Determines if a 3D triangle intersects a 3D ray.
 * \param [in] tri A 3D triangle
 * \param [in] ray A 3D ray
 * \param [out] t Intersection point of tri and R, w.r.t. parametrization of R
 * \param [out] p Intersection point of tri and R, in barycentric coordinates
 *   relative to tri.
 * \note If there is an intersection, the intersection point is:  R.at(t)
 * \return true iff tri intersects with ray, otherwise, false.
 * \note \a t and \a p only valid when function returns true
 */
template <typename T>
bool intersect(const Triangle<T, 3>& tri,
               const Ray<T, 3>& ray,
               T& t,
               Point<double, 3>& p)
{
  bool retval = detail::intersect_tri_ray(tri, ray, t, p);

  if(retval)
  {
    // Add a small EPS to avoid dividing by zero
    double normalizer = p[0] + p[1] + p[2] + primal::PRIMAL_TINY;
    p.array() *= 1. / normalizer;
  }

  return retval;
}

/*!
 * \brief Determines if a 3D triangle intersects a 3D segment.
 * \param [in] tri A 3D triangle
 * \param [in] seg A 3D line segment
 * \return true iff tri intersects with seg, otherwise, false.
 */
template <typename T>
bool intersect(const Triangle<T, 3>& tri, const Segment<T, 3>& seg)
{
  T t = T();
  Point<double, 3> p;
  return detail::intersect_tri_segment(tri, seg, t, p);
}

/*!
 * \brief Determines if a 3D triangle intersects a 3D segment.
 * \param [in] tri A 3D triangle
 * \param [in] seg A 3D line segment
 * \param [out] t Intersection point of tri and seg, w.r.t. seg's
 *  parametrization
 * \return true iff tri intersects with seg, otherwise, false.
 */
template <typename T>
bool intersect(const Triangle<T, 3>& tri, const Segment<T, 3>& seg, T& t)
{
  Point<double, 3> p;
  return detail::intersect_tri_segment(tri, seg, t, p);
}

/*!
 * \brief Determines if a 3D triangle intersects a 3D segment.
 * \param [in] tri A 3D triangle
 * \param [in] seg A 3D line segment
 * \param [out] t Intersection point of tri and seg, w.r.t. seg's
 *  parametrization
 * \param [out] p Intersection point of tri and R, in barycentric coordinates
 *   relative to tri.
 * \note If there is an intersection, the intersection point pt is:
 *                     pt = seg.source() + t * ( seg.dest() - seg.target() )
 * \return true iff tri intersects with seg, otherwise, false.
 * \note \a t and \a p only valid when function returns true
 */
template <typename T>
bool intersect(const Triangle<T, 3>& tri,
               const Segment<T, 3>& seg,
               T& t,
               Point<double, 3>& p)
{
  bool retval = detail::intersect_tri_segment(tri, seg, t, p);

  if(retval)
  {
    // Add a small EPS to avoid dividing by zero
    double normalizer = p[0] + p[1] + p[2] + primal::PRIMAL_TINY;
    p.array() *= 1. / normalizer;
  }

  return retval;
}

/// @}

/// \name Ray Intersection Routines
/// @{

/*!
 * \brief Computes the intersection of the given ray, R, with the segment, S.
 *
 * \param [in] R the specified ray
 * \param [in] S the segment to check
 * \param [out] ray_param parametric coordinate of intersection along R, valid only if status=true.
 * \param [out] seg_param parametric coordinate of intersection along S, valid only if status=true.
 * \param [in] EPS tolerance for intersection tests
 *
 * \return status true iff R intersects with S, otherwise, false.
 *
 * \see primal::Ray
 * \see primal::Segment
 */
template <typename T>
bool intersect(const Ray<T, 2>& R,
               const Segment<T, 2>& S,
               T& ray_param,
               T& seg_param,
               const T EPS = 1e-8)
{
  return detail::intersect_ray(R, S, ray_param, seg_param, EPS);
}

/*!
 * \brief Computes the intersection of the given ray, R, with the segment, S.
 *
 * \param [in] R the specified ray
 * \param [in] S the segment to check
 * \param [out] ray_param parametric coordinate of intersection along R, valid only if status=true
 *
 * \note If you need to specify a tolerance for the intersection tests, please use the overload
 * of this function with two [OUT] parameters (\a ray_param and \a seg_param)
 *
 * \return status true iff R intersects with S, otherwise, false.
 *
 * \see primal::Ray
 * \see primal::Segment
 */
template <typename T>
bool intersect(const Ray<T, 2>& R, const Segment<T, 2>& S, T& ray_param)
{
  T seg_param;
  return intersect(R, S, ray_param, seg_param);
}

/*!
 * \brief Computes the intersection of the given ray, R, with the segment, S.
 *
 * \param [in] R the specified ray
 * \param [in] S the segment to check
 * \param [out] ip the intersection point on S, valid only if status=true.
 * \param [in] EPS tolerance for intersection tests
 *
 * \return status true iff R intersects with S, otherwise, false.
 *
 * \see primal::Ray
 * \see primal::Segment
 * \see primal::Point
 */
template <typename T>
bool intersect(const Ray<T, 2>& R,
               const Segment<T, 2>& S,
               Point<T, 2>& ip,
               const T EPS = 1e-8)
{
  T ray_param;
  T seg_param;
  if(detail::intersect_ray(R, S, ray_param, seg_param, EPS))
  {
    ip = R.at(ray_param);
    return true;
  }
  return false;
}

/*!
 * \brief Computes the intersection of the given ray, R, with the Box, bb.
 *
 * \param [in] R the specified ray
 * \param [in] bb the user-supplied axis-aligned bounding box
 *
 * \param [out] ip the intersection point where R intersects bb.
 *
 * \return status true iff bb intersects with R, otherwise, false.
 *
 * \see primal::Ray
 * \see primal::Segment
 * \see primal::BoundingBox
 *
 * \note Computes Ray Box intersection using the slab method from pg 180 of
 *  Real Time Collision Detection by Christer Ericson.
 */
template <typename T, int DIM>
AXOM_HOST_DEVICE bool intersect(const Ray<T, DIM>& R,
                                const BoundingBox<T, DIM>& bb,
                                Point<T, DIM>& ip)
{
  return detail::intersect_ray(R, bb, ip);
}

/*!
 * \brief Computes the intersection of the given line, L, with the Box, bb.
 *
 * \param [in] L the specified line (two-sided ray)
 * \param [in] bb the user-supplied axis-aligned bounding box
 *
 * \param [out] ip the intersection point where L intersects bb.
 *
 * \return status true iff bb intersects with R, otherwise, false.
 *
 * \see primal::Line
 * \see primal::Segment
 * \see primal::BoundingBox
 *
 * \note Computes Ray Box intersection using the slab method from pg 180 of
 *  Real Time Collision Detection by Christer Ericson.
 */
template <typename T, int DIM>
AXOM_HOST_DEVICE bool intersect(const Line<T, DIM>& L,
                                const BoundingBox<T, DIM>& bb,
                                Point<T, DIM>& ip)
{
  return detail::intersect_line(L, bb, ip);
}
/// @}

/// \name Segment-BoundingBox Intersection Routines
/// @{

/*!
 * \brief Computes the intersection of the given segment, \a S, with the Box, \a bb.
 *     If an intersection is found, output parameter \a ip contains an intersection point
 * \return status true iff \a bb intersects with \a S, otherwise, false.
 *
 * \note The intersection between segment \a S and box \a bb intersect, will, in general,
 * be a along a (1D) subset of segment \a S. One variant of this function returns the two
 * parametric coordinates of the intersections along \a S found while determining 
 * if there is a valid intersection. Another variant returns an intersection point along \a S
 * Specifically, it is the point of smallest parametric coordinate that is contained in \a bb
 * (i.e. with parameter \a tmin). These are only valid when the function returns true
 * 
 * Computes Segment-Box intersection using the slab method from pg 180 of
 * Real Time Collision Detection by Christer Ericson.
 */

/// This variant returns the two parametric coordinates of the intersection segment as OUT parameters
template <typename T, int DIM>
bool intersect(const Segment<T, DIM>& S,
               const BoundingBox<T, DIM>& bb,
               T& tmin,
               T& tmax,
               const double& EPS = 1e-8)
{
  const T segLength = S.length();
  tmin = static_cast<T>(0);
  tmax = static_cast<T>(segLength);

  bool intersects = segLength > 0. &&
    detail::intersect_ray(Ray<T, DIM>(S), bb, tmin, tmax, EPS);

  // Scale parametric coordinates with respect to the segment
  if(intersects)
  {
    tmin /= segLength;
    tmax /= segLength;
  }

  return intersects;
}

/// This variant returns a point within the intersection as an OUT parameters
template <typename T, int DIM>
bool intersect(const Segment<T, DIM>& S,
               const BoundingBox<T, DIM>& bb,
               Point<T, DIM>& ip,
               const double& EPS = 1e-8)
{
  T tmin, tmax;
  if(intersect(S, bb, tmin, tmax, EPS))
  {
    ip = S.at(tmin);
    return true;
  }
  return false;
}

template <typename T, int DIM>
bool intersect(const Segment<T, DIM>& S, const BoundingBox<T, DIM>& bb)
{
  T tmin, tmax;
  return intersect(S, bb, tmin, tmax);
}

/// @}

/// \name Axis-Aligned Bounding Box Intersection Routines
/// @{

/*!
 * \brief Determines if two axis aligned bounding boxes intersect
 * \param [in] bb1 user-supplied axis aligned bounding box.
 * \param [in] bb2 user-supplied axis aligned bounding box.
 * \return true iff bb1 intersects with bb2, otherwise, false.
 */
template <typename T, int DIM>
AXOM_HOST_DEVICE bool intersect(const BoundingBox<T, DIM>& bb1,
                                const BoundingBox<T, DIM>& bb2)
{
  return bb1.intersectsWith(bb2);
}

/// @}

/// \name Sphere Intersection Routines
/// @{

/*!
 * \brief Determines if two spheres intersect.
 *
 * \param [in] s1 user-supplied sphere object to check for intersection.
 * \param [in] s2 user-supplied sphere object to check for intersection.
 * \param [in] TOL tolerance used for intersection check (optional)
 *
 * \note If TOL is not supplied, the default is 1.e-9.
 *
 * \return status true iff s1 intersects with s2, otherwise, false.
 */
template <typename T, int DIM>
bool intersect(const Sphere<T, DIM>& s1,
               const Sphere<T, DIM>& s2,
               double TOL = 1.e-9)
{
  return s1.intersectsWith(s2, TOL);
}
/// @}

/// \name Oriented Bounding Box Intersection Routines
/// @{

template <typename T>
bool intersect(const OrientedBoundingBox<T, 1>& b1,
               const OrientedBoundingBox<T, 1>& b2)
{
  return detail::intersect_obb1D_obb1D(b1, b2);
}

/*!
 * \brief Determines if a 2D OBB intersects a 2D OBB.
 * \param [in] b1 A 2D OrientedBoundingBox
 * \param [in] b2 A 2D OrientedBoundingBox
 * \return true iff b1 intersects with b2, otherwise, false.
 */
template <typename T>
bool intersect(const OrientedBoundingBox<T, 2>& b1,
               const OrientedBoundingBox<T, 2>& b2)
{
  return detail::intersect_obb2D_obb2D(b1, b2);
}

/*!
 * \brief Determines if a 3D OBB intersects a 3D OBB.
 * \param [in] b1 A 3D OrientedBoundingBox
 * \param [in] b2 A 3D OrientedBoundingBox
 * \param [in] EPS error tolerance for intersection
 * \return true iff b1 intersects with b2, otherwise, false.
 */
template <typename T>
bool intersect(const OrientedBoundingBox<T, 3>& b1,
               const OrientedBoundingBox<T, 3>& b2,
               double EPS = 1E-4)
{
  return detail::intersect_obb3D_obb3D(b1, b2, EPS);
}

/// @}

/// \name Bezier Curve Intersection Routines
/// @{

/*!
 * \brief Tests if two Bezier Curves \a c1 and \a c2 intersect.
 * \return status true iff \a c1 intersects \a c2, otherwise false.
 *
 * \param [in] c1 the first BezierCurve, parametrized in [0,1)
 * \param [in] c2 the second BezierCurve, parametrized in [0,1)
 * \param [out] sp vector of parameter space intersection points for \a c1
 * \param [out] tp vector of parameter space intersection points for \a c2
 * \param [in] tol tolerance parameter for determining if a curve can
 * be approximated by a line segment.
 * \return True if the curves intersect, false otherwise. Intersection
 * parameters are stored in \a sp and \a tp
 *
 * Finds all intersection points between the two curves.
 *
 * \note This function assumes two dimensional curves in a plane.
 *
 * \note This function assumes that the curves are in general position.
 * Specifically, we assume that all intersections are at points and that
 * the curves don't overlap.
 *
 * \note This function assumes the all intersections have multiplicity
 * one, i.e. there are no points at which the curves and their derivatives
 * both intersect. Thus, the function does not find tangencies.
 *
 * \note This function assumes that the curves are half-open, i.e. they
 * contain their first endpoint, but not their last endpoint. Thus, the
 * curves do not intersect at \f$ s==1 \f$ or at \f$ t==1 \f$.
 */
template <typename T>
bool intersect(const BezierCurve<T, 2>& c1,
               const BezierCurve<T, 2>& c2,
               axom::Array<T>& sp,
               axom::Array<T>& tp,
               double tol = 1E-8)
{
  const double offset = 0.;
  const double scale = 1.;

  // for efficiency, linearity check actually uses a squared tolerance
  const double sq_tol = tol * tol;

  return detail::intersect_bezier_curves(c1,
                                         c2,
                                         sp,
                                         tp,
                                         sq_tol,
                                         c1.getOrder(),
                                         c2.getOrder(),
                                         offset,
                                         scale,
                                         offset,
                                         scale);
}

/*!
 * \brief Function to find intersections between a ray and a Bezier curve
 *
 * \param [in] r The input ray
 * \param [in] c The input curve
 * \param [out] rp Parametric coordinates of intersections in \a r [0, inf)
 * \param [out] cp Parametric coordinates of intersections in \a c [0, 1)
 * Bezier curve is linear
 * \param [in] tol Tolerance parameter for physical distances
 * \param [in] EPS Tolerance parameter for parameter-space distances
 * 
 * \note A BezierCurve is parametrized in [0,1). This function assumes the all
 *  intersections have multiplicity one, i.e. the function does not find tangencies.
 * 
 * \return True if the ray intersects the Bezier curve, False otherwise
 */
template <typename T>
bool intersect(const Ray<T, 2>& r,
               const BezierCurve<T, 2>& c,
               axom::Array<T>& rp,
               axom::Array<T>& cp,
               double tol = 1E-8,
               double EPS = 1E-8)
{
  const double offset = 0.;
  const double scale = 1.;

  // for efficiency, linearity check actually uses a squared tolerance
  const double sq_tol = tol * tol;

  return detail::intersect_ray_bezier(r,
                                      c,
                                      rp,
                                      cp,
                                      sq_tol,
                                      EPS,
                                      c.getOrder(),
                                      offset,
                                      scale);
}

/*!
 * \brief Function to find intersections between a ray and a NURBS curve
 *
 * \param [in] r The input ray
 * \param [in] n The input curve
 * \param [out] rp Parametric coordinates of intersections in \a r [0, inf)
 * \param [out] cp Parametric coordinates of intersections in the knot span of \a n
 * \param [in] tol Tolerance parameter for physical distances
 * \param [in] EPS Tolerance parameter for parameter-space distances
 *
 * \note Assumes the NURBS curve is parameterized on a half-open interval [a, b),
 *  and assumes the all intersections have multiplicity one, i.e. the function does not find tangencies.
 * 
 * \return True if the ray intersects the NURBS curve, False otherwise
 */
template <typename T>
bool intersect(const Ray<T, 2>& r,
               const NURBSCurve<T, 2>& n,
               axom::Array<T>& rp,
               axom::Array<T>& np,
               double tol = 1E-8,
               double EPS = 1E-8)
{
  // Check a bounding box of the entire NURBS first
  Point<T, 2> ip;
  if(!intersect(r, n.boundingBox(), ip))
  {
    return false;
  }

  // Decompose the NURBS curve into Bezier segments
  auto beziers = n.extractBezier();
  axom::Array<T> knot_vals = n.getKnots().getUniqueKnots();

  // Check each Bezier segment, and scale the intersection parameters
  //  back into the span of the original NURBS curve
  for(int i = 0; i < beziers.size(); ++i)
  {
    axom::Array<T> rc, nc;
    intersect(r, beziers[i], rc, nc, tol, EPS);

    // Scale the intersection parameters back into the span of the NURBS curve
    for(int j = 0; j < rc.size(); ++j)
    {
      rp.push_back(rc[j]);
      np.push_back(knot_vals[i] + nc[j] * (knot_vals[i + 1] - knot_vals[i]));
    }
  }

  return !rp.empty();
}

/*!
 * \brief Function to determine *if* a NURBS curve intersects a bounding box
 * \note Does not attempt to compute the actual points of intersection
 */
template <typename T>
bool intersects(const BoundingBox<T, 2>& bb, const NURBSCurve<T, 2>& n)
{
  // If the bounding boxes don't intersect, no intersection is possible
  if(!intersect(bb, n.boundingBox()))
  {
    return false;
  }

  // If the bounding box contains the entire curve, no intersection is possible
  if(bb.contains(n.boundingBox()))
  {
    return false;
  }

  // If the first control point or the last control point *touches* the bounding box,
  //  then that counts as an intersection
  auto first = n[0];
  if(first[0] == bb.getMin()[0] || first[0] == bb.getMax()[0] ||
     first[1] == bb.getMin()[1] || first[1] == bb.getMax()[1])
  {
    return true;
  }

  auto last = n[n.getNumControlPoints() - 1];
  if(last[0] == bb.getMin()[0] || last[0] == bb.getMax()[0] ||
     last[1] == bb.getMin()[1] || last[1] == bb.getMax()[1])
  {
    return true;
  }

  // If the first control point is in the bounding box and the last control point is outside,
  //  or vice versa, then there is guaranteed to be an intersection
  bool first_in = bb.contains(first);
  bool last_in = bb.contains(last);
  if((first_in && !last_in) || (!first_in && last_in))
  {
    return true;
  }

  // Otherwise, need to split the curve in two and check each half
  NURBSCurve<T, 2> n1, n2;
  n.bisect(n1, n2);

  // Since we are only looking if intersections exist, can short circuit
  return intersects(bb, n1) || intersects(bb, n2);
}

/// @}

/// \name Plane Intersection Routines
/// @{

/*!
 * \brief Determines if a 3D plane intersects a 3D bounding box.
 *        By default (checkOverlaps is false), checks if |s| <= r, 
 *        where "s" is the distance of the bounding box center to the plane,
 *        and "r" is the projected radius of the bounding box along the line
 *        parallel to the plane normal and going through the box center.
 *        If checkOverlaps is true, checks if |s| < r,
 *        where the bounding box overlaps both half spaces of the plane.
 * \param [in] p A 3D plane
 * \param [in] bb A 3D bounding box
 * \param [in] checkOverlaps If true, checks if bounding box overlaps both 
 *             halfspaces of the plane.
 *             Otherwise, overlap of both halfspaces is not guaranteed.
 *             Default is false.
 * \param [in] EPS tolerance parameter for determining if "s"
 *             is just within min/max of "r".
 * \return true iff plane intersects with bounding box, otherwise, false.
 *
 * \note Uses method from pg 164 of 
 *       Real Time Collision Detection by Christer Ericson.
 */
template <typename T>
AXOM_HOST_DEVICE bool intersect(const Plane<T, 3>& p,
                                const BoundingBox<T, 3>& bb,
                                bool checkOverlaps = false,
                                double EPS = 1E-08)
{
  return detail::intersect_plane_bbox(p, bb, checkOverlaps, EPS);
}

/*!
 * \brief Determines if a plane intersects a segment.
 * \param [in] plane A plane
 * \param [in] seg A line segment
 * \param [out] t Intersection point of plane and seg, w.r.t. seg's
 *  parametrization
 * \note If there is an intersection, the intersection point pt is:
 *                     pt = seg.at(t)
 * \return true iff plane intersects with seg, otherwise, false.
 * \note \a t is only valid when function returns true
 *
 * \note Uses method from pg 176 of 
 *       Real Time Collision Detection by Christer Ericson.
 */
template <typename T, int DIM>
AXOM_HOST_DEVICE bool intersect(const Plane<T, DIM>& plane,
                                const Segment<T, DIM>& seg,
                                T& t)
{
  return detail::intersect_plane_seg(plane, seg, t);
}

/*!
 * \brief Determines if a 3D plane intersects a tetrahedron.
 *
 * \param [in] p A 3D plane
 * \param [in] tet A 3D tetrahedron
 * \param [out] intersection A polygon containing the intersection.
 *
 * \return true if plane intersects with tetrahedron, otherwise, false.
 *
 * \note If no intersection is found, the output polygon will be empty.
 *       If the plane intersects at a tetrahedron vertex, the polygon
 *       will contain duplicated points.
 */
template <typename T>
AXOM_HOST_DEVICE bool intersect(const Plane<T, 3>& p,
                                const Tetrahedron<T, 3>& tet,
                                Polygon<T, 3>& intersection)
{
  return detail::intersect_plane_tet3d(p, tet, intersection);
}

/*! \brief Determines if a line intersects a bilinear patch.
 * \param [in] patch The bilinear (bezier) patch to intersect with the ray.
 * \param [in] line The line to intersect with the bilinear patch.
 * \param [out] u The u patch parameter(s) of the intersection point.
 * \param [out] v The v patch parameter(s) of the intersection point.
 * \param [out] t The line parameter(s) of the intersection point.
 * \param [in] EPS The tolerance for intersection.
 *
 * Implements GARP algorithm from Chapter 8 of Ray Tracing Gems (2019)
 * 
 * \note Assumes the patch is linear in either axis, but will work for any
 *  order patch using its corners.
 * 
 * \return true iff the ray intersects the bilinear patch, otherwise false.
 */
template <typename T>
AXOM_HOST_DEVICE bool intersect(const Point<T, 3>& p0,
                                const Point<T, 3>& p1,
                                const Point<T, 3>& p2,
                                const Point<T, 3>& p3,
                                const Line<T, 3>& line,
                                std::vector<T>& u,
                                std::vector<T>& v,
                                std::vector<T>& t)
{
  SLIC_ERROR("Not implemented!");
  AXOM_UNUSED_VAR(p0);
  AXOM_UNUSED_VAR(p1);
  AXOM_UNUSED_VAR(p2);
  AXOM_UNUSED_VAR(p3);
  AXOM_UNUSED_VAR(line);
  AXOM_UNUSED_VAR(u);
  AXOM_UNUSED_VAR(v);
  AXOM_UNUSED_VAR(t);
  // const int order_u = patch.getOrder_u();
  // const int order_v = patch.getOrder_v();
  // return detail::intersect_bilinear_patch_ray(patch(0, 0),
  //                                             patch(order_u, 0),
  //                                             patch(order_u, order_v),
  //                                             patch(0, order_v),
  //                                             ray,
  //                                             u,
  //                                             v,
  //                                             t);
}

/// @}

/*! \brief Determines if a line intersects a Bezier patch.
 * \param [in] patch The Bezier patch to intersect with the line.
 * \param [in] line The line to intersect with the patch.
 * \param [out] u The u parameter(s) of intersection point(s).
 * \param [out] v The v parameter(s) of intersection point(s).
 * \param [out] t The t parameter(s) of intersection point(s).
 * \param [in] tol The tolerance for intersection (for physical distances).
 * \param [in] EPS The tolerance for intersection (for parameter distances).
 * \param [in] isHalfOpen True if the patch is parameterized in [0,1)^2.
 * 
 * For bilinear patches, implements GARP algorithm from Chapter 8 of Ray Tracing Gems (2019)
 * For higher order patches, intersections are found through recursive subdivison
 *  until the subpatch is approximated by a bilinear patch.
 * Assumes that the line is not tangent to the patch, and that the intersection
 *  is not at a point of degeneracy for which there are *infinitely* many intersections.
 * For such intersections, the method will hang as it tries records an arbitrarily high
 *  number of intersections with distinct parameter values
 *  
 * \return true iff the line intersects the patch, otherwise false.
 */
template <typename T>
AXOM_HOST_DEVICE bool intersect(const Line<T, 3>& line,
                                const BezierPatch<T, 3>& patch,
                                axom::Array<T>& t,
                                axom::Array<T>& u,
                                axom::Array<T>& v,
                                double tol = 1e-8,
                                double EPS = 1e-8,
                                bool isHalfOpen = false)
{
  const int order_u = patch.getOrder_u();
  const int order_v = patch.getOrder_v();

  // for efficiency, linearity check actually uses a squared tolerance
  const double sq_tol = tol * tol;

  // Store the candidate intersections
  axom::Array<T> tc, uc, vc;

  bool success = false;

  if(order_u < 1 || order_v < 1)
  {
    // Patch has no surface area, ergo no intersections
    return true;
  }
  else if(order_u == 1 && order_v == 1)
  {
    detail::intersect_line_bilinear_patch(line,
                                          patch(0, 0),
                                          patch(order_u, 0),
                                          patch(order_u, order_v),
                                          patch(0, order_v),
                                          t,
                                          u,
                                          v,
                                          EPS,
                                          false);
    return true;
  }
  else
  {
    // primal::Line<T, 3> line(ray.origin(), ray.direction());

    double u_offset = 0., v_offset = 0.;
    double u_scale = 1., v_scale = 1.;

    success = detail::intersect_line_patch(line,
                                           patch,
                                           tc,
                                           uc,
                                           vc,
                                           order_u,
                                           order_v,
                                           u_offset,
                                           u_scale,
                                           v_offset,
                                           v_scale,
                                           sq_tol,
                                           EPS,
                                           false);
  }

  // Remove duplicates from the (u, v) intersection points
  //  (Note it's not possible for (u_1, v_1) == (u_2, v_2) and t_1 != t_2)
  const double sq_EPS = EPS * EPS;

  // The number of reported intersection points will be small,
  //  so we don't need to fully sort the list
  for(int i = 0; i < tc.size(); ++i)
  {
    // Also remove any intersections on the half-interval boundaries
    if(isHalfOpen && (uc[i] >= 1.0 - EPS || vc[i] >= 1.0 - EPS))
    {
      continue;
    }

    Point<T, 2> uv({uc[i], vc[i]});

    bool foundDuplicate = false;
    for(int j = i + 1; !foundDuplicate && j < tc.size(); ++j)
    {
      if(squared_distance(uv, Point<T, 2>({uc[j], vc[j]})) < sq_EPS)
      {
        foundDuplicate = true;
      }
    }

    if(!foundDuplicate)
    {
      t.push_back(tc[i]);
      u.push_back(uc[i]);
      v.push_back(vc[i]);
    }
  }

  // If too many intersections are reported, it can be due to degeneracies in the patch.
  //  To account for this, we filter out duplicates in the ray's parameter space
  if(!success)
  {
    tc.clear();
    uc.clear();
    vc.clear();

    // Remove duplicates from the (t) intersection points
    //  (Note that while it's possible for S(u_2, v_2) == S(u_2, v_2) and for t_1 != t_2,
    //  it's very unlikely for this to occur upwards of 25 times per patch)
    for(int i = 0; i < t.size(); ++i)
    {
      bool foundDuplicate = false;
      for(int j = i + 1; !foundDuplicate && j < t.size(); ++j)
      {
        if(squared_distance(line.at(t[i]), line.at(t[j])) < sq_EPS)
        {
          foundDuplicate = true;
        }
      }

      if(!foundDuplicate)
      {
        tc.push_back(t[i]);
        uc.push_back(u[i]);
        vc.push_back(v[i]);
      }
    }

    t = tc;
    u = uc;
    v = vc;
  }

  return success;
}

template <typename T>
AXOM_HOST_DEVICE bool intersect(const Line<T, 3>& line,
                                const NURBSPatchData<T>& patch_data,
                                axom::Array<T>& t,
                                axom::Array<T>& u,
                                axom::Array<T>& v,
                                double tol = 1e-8,
                                double EPS = 1e-8,
                                bool isHalfOpen = false,
                                bool isTrimmed = true,
                                double buffer = 0.0)
{
  // Check a bounding box of the entire NURBS first
  Point<T, 3> ip;
  if(!intersect(line, BoundingBox<T, 3>(patch_data.bbox).scale(1.5), ip))
  {
    return true;
  }

  // Store candidate intersections
  axom::Array<T> tc, uc, vc;

  // std::cout << "Avoided a bezier extraction! (" << patch_data.beziers.size() << ") " << std::endl;

  // Check each (precomputed) Bezier patch, and scale the intersection parameters
  //  back into the span of the original NURBS patch
  bool success = false;
  for(int ij = 0; ij < patch_data.beziers.size(); ++ij)
  {
    // Store candidate intersections from each Bezier patch
    axom::Array<T> tcc, ucc, vcc;
    // std::cout << patch_data.beziers[ij] << std::endl;
    success = intersect(line, patch_data.beziers[ij], tcc, ucc, vcc, tol, EPS);

    // Scale the intersection parameters back into the span of the NURBS patch
    for(int k = 0; k < tcc.size(); ++k)
    {
      // std::cout << patch_data.u_spans[ij].first << ", " << patch_data.u_spans[ij].second << std::endl;
      tc.push_back(tcc[k]);
      uc.push_back(
        patch_data.u_spans[ij].first +
        ucc[k] * (patch_data.u_spans[ij].second - patch_data.u_spans[ij].first));
      vc.push_back(
        patch_data.v_spans[ij].first +
        vcc[k] * (patch_data.v_spans[ij].second - patch_data.v_spans[ij].first));
    }
  }

  // Do a second pass to remove duplicates from uc, vc
  const double sq_EPS = EPS * EPS;

  // The number of reported intersection points will be small,
  //  so we don't need to fully sort the list

  double max_u_knot =
    patch_data.patch.getKnots_u()[patch_data.patch.getKnots_u().getNumKnots() - 1];
  double max_v_knot =
    patch_data.patch.getKnots_v()[patch_data.patch.getKnots_v().getNumKnots() - 1];

  for(int i = 0; i < tc.size(); ++i)
  {
    // Also remove any intersections on the half-interval boundaries
    if(isHalfOpen &&
       (uc[i] >= max_u_knot - buffer || vc[i] >= max_v_knot - buffer))
    {
      continue;
    }

    // And remove any intersections that are not visible
    if(isTrimmed && !patch_data.patch.isVisible(uc[i], vc[i]))
    {
      continue;
    }

    Point<T, 2> uv({uc[i], vc[i]});

    bool foundDuplicate = false;
    for(int j = i + 1; !foundDuplicate && j < tc.size(); ++j)
    {
      if(squared_distance(uv, Point<T, 2>({uc[j], vc[j]})) < sq_EPS)
      {
        foundDuplicate = true;
      }
    }

    if(!foundDuplicate)
    {
      t.push_back(tc[i]);
      u.push_back(uc[i]);
      v.push_back(vc[i]);
    }
  }

  return success;
}

template <typename T>
AXOM_HOST_DEVICE bool intersect(const Line<T, 3>& line,
                                const NURBSPatch<T, 3>& patch,
                                axom::Array<T>& t,
                                axom::Array<T>& u,
                                axom::Array<T>& v,
                                double tol = 1e-8,
                                double EPS = 1e-8,
                                bool isHalfOpen = false,
                                bool isTrimmed = true,
                                double buffer = 0.0)
{
  // Check a bounding box of the entire NURBS first
  Point<T, 3> ip;
  if(!intersect(line, patch.boundingBox().scale(1.5), ip))
  {
    return true;
  }

  // Decompose the NURBS patch into Bezier patches
  auto beziers = patch.extractBezier();

  axom::Array<T> knot_vals_u = patch.getKnots_u().getUniqueKnots();
  axom::Array<T> knot_vals_v = patch.getKnots_v().getUniqueKnots();

  const auto num_knot_span_u = knot_vals_u.size() - 1;
  const auto num_knot_span_v = knot_vals_v.size() - 1;

  // Store candidate intersections
  axom::Array<T> tc, uc, vc;

  // Check each Bezier patch, and scale the intersection parameters
  //  back into the span of the original NURBS patch
  bool success = false;
  for(int i = 0; i < num_knot_span_u; ++i)
  {
    for(int j = 0; j < num_knot_span_v; ++j)
    {
      auto& bezier = beziers[i * num_knot_span_v + j];

      // Store candidate intersections from each Bezier patch
      axom::Array<T> tcc, ucc, vcc;
      success = intersect(line, bezier, tcc, ucc, vcc, tol, EPS);

      // Scale the intersection parameters back into the span of the NURBS patch
      for(int k = 0; k < tcc.size(); ++k)
      {
        tc.push_back(tcc[k]);
        uc.push_back(knot_vals_u[i] +
                     ucc[k] * (knot_vals_u[i + 1] - knot_vals_u[i]));
        vc.push_back(knot_vals_v[j] +
                     vcc[k] * (knot_vals_v[j + 1] - knot_vals_v[j]));
      }
    }
  }

  // Do a second pass to remove duplicates from uc, vc
  const double sq_EPS = EPS * EPS;

  // The number of reported intersection points will be small,
  //  so we don't need to fully sort the list

  double max_u_knot = patch.getKnots_u()[patch.getKnots_u().getNumKnots() - 1];
  double max_v_knot = patch.getKnots_v()[patch.getKnots_v().getNumKnots() - 1];

  for(int i = 0; i < tc.size(); ++i)
  {
    // Also remove any intersections on the half-interval boundaries
    if(isHalfOpen &&
       (uc[i] >= max_u_knot - buffer || vc[i] >= max_v_knot - buffer))
    {
      continue;
    }

    // And remove any intersections that are not visible
    if(isTrimmed && !patch.isVisible(uc[i], vc[i]))
    {
      continue;
    }

    Point<T, 2> uv({uc[i], vc[i]});

    bool foundDuplicate = false;
    for(int j = i + 1; !foundDuplicate && j < tc.size(); ++j)
    {
      if(squared_distance(uv, Point<T, 2>({uc[j], vc[j]})) < sq_EPS)
      {
        foundDuplicate = true;
      }
    }

    if(!foundDuplicate)
    {
      t.push_back(tc[i]);
      u.push_back(uc[i]);
      v.push_back(vc[i]);
    }
  }

  return success;
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_INTERSECT_HPP_

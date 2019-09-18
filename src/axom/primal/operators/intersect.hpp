// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
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

#include "axom/core/Macros.hpp"
#include "axom/core/utilities/Utilities.hpp"

#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/OrientedBoundingBox.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Ray.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/Sphere.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/geometry/CurvedPolygon.hpp"

#include "axom/primal/operators/detail/intersect_impl.hpp"
#include "axom/primal/operators/detail/intersect_ray_impl.hpp"
#include "axom/primal/operators/detail/intersect_bounding_box_impl.hpp"
#include "axom/primal/operators/detail/intersect_bezier_impl.hpp"
#include "axom/primal/operators/detail/intersect_curved_poly_impl.hpp"

namespace axom
{
namespace primal
{
/// \name Triangle Intersection Routines
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
    const double EPS = 1e-80;
    double normalizer = p[0] + p[1] + p[2] + EPS;
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
    const double EPS = 1e-80;
    double normalizer = p[0] + p[1] + p[2] + EPS;
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
bool intersect(const BoundingBox<T, DIM>& bb1, const BoundingBox<T, DIM>& bb2)
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
 * \brief Tests if two CurvedPolygon \a p1 and \a p2 intersect.
 * \return status true iff \a p1 intersects \a p2, otherwise false.
 *
 * \param [in] p1 the first CurvedPolygon 
 * \param [in] p2 the second CurvedPolygon
 * \param [out] pnew vector of resulting CurvedPolygons
 * \return True if the curvedPolygons intersect, false otherwise. 
 *
 * Finds the set of curved polygons that bound the region of intersection between p1 and p2.
 *
 * \note This function assumes two dimensional curved polygons  in a plane.
 *
 * \note This function assumes that the curved polygons are in general 
 * position. Specifically, we assume that all intersections are at points 
 * and that the component curves don't overlap.
 *
 * \note This function assumes the all intersections have multiplicity 
 * one, i.e. there are no points at which the component curves and their 
 * derivatives both intersect. Thus, the function does not find tangencies.
 *
 * \note This function assumes that the component curves are half-open, 
 * i.e. they contain their first endpoint, but not their last endpoint. 
 * Thus, component curves do not intersect at \f$ s==1 \f$ or at 
 * \f$ t==1 \f$.
 */
template <typename T, int NDIMS>
bool intersect(CurvedPolygon<T, NDIMS>& p1,
               CurvedPolygon<T, NDIMS>& p2,
               std::vector<CurvedPolygon<T, NDIMS>>& pnew,
               double tol = 1e-8)
{
  // for efficiency, linearity check actually uses a squared tolerance
  const double sq_tol = tol * tol;

  return detail::intersect_polygon(p1, p2, pnew, sq_tol);
}

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
template <typename T, int NDIMS>
bool intersect(const BezierCurve<T, NDIMS>& c1,
               const BezierCurve<T, NDIMS>& c2,
               std::vector<T>& sp,
               std::vector<T>& tp,
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

/// @}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_INTERSECT_HPP_

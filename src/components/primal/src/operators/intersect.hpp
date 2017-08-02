/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

/*!
 * \file intersect.hpp
 *
 * \brief Consists of functions to test intersection among geometric primitives.
 */

#ifndef INTERSECTION_HPP_
#define INTERSECTION_HPP_

#include "primal/OrientedBoundingBox.hpp"
#include "primal/BoundingBox.hpp"
#include "primal/Point.hpp"
#include "primal/Ray.hpp"
#include "primal/Segment.hpp"
#include "primal/Triangle.hpp"

#include "primal/intersect_impl.hpp"

namespace axom {
namespace primal {

/*!
 * \brief Tests if 3D Triangles t1 and t2 intersect.
 * \return status true iff t1 intersects with t2, otherwise, false.
 *
 * If parameter includeBoundary is false (default), this function will
 * return true if the interior of t1 intersects the interior of t2.  To include
 * triangle boundaries in intersections, specify includeBoundary as true.
 */
template < typename T >
bool intersect( const Triangle< T, 3 >& t1,
                const Triangle< T, 3 >& t2,
                const bool includeBoundary = false)
{
  return detail::intersect_tri3D_tri3D< T >(t1, t2, includeBoundary);
}

/*!
 * \brief Tests if 2D Triangles t1 and t2 intersect.
 * \return status true iff t1 intersects with t2, otherwise, false.
 *
 * If parameter includeBoundary is false (default), this function will
 * return true if the interior of t1 intersects the interior of t2.  To include
 * triangle boundaries in intersections, specify includeBoundary as true.
 */
template < typename T >
bool intersect( const Triangle< T, 2 >& t1,
                const Triangle< T, 2 >& t2,
                const bool includeBoundary = false)
{
  return detail::intersect_tri2D_tri2D< T >(t1, t2, includeBoundary);
}

/*!
 * \brief Computes the intersection of the given ray, R, with the segment, S.
 *      ip returns the intersection point on S.
 * \return status true iff R intersects with S, otherwise, false.
 */
template < typename T >
bool intersect( const Ray< T,2 >& R, const Segment< T,2 >& S, Point< T,2 >& ip )
{
  return detail::intersect_ray_seg(R, S, ip);
}

/*!
 * \brief Computes the intersection of the given ray, R, with the Box, bb.
 *      ip the point of intersection on R.
 * \return status true iff bb intersects with R, otherwise, false.
 *
 * Computes Ray Box intersection using the slab method from pg 180 of
 * Real Time Collision Detection by Christer Ericson.
 */
template < typename T, int DIM >
bool intersect( const Ray< T,DIM > & R,
                const BoundingBox< T,DIM > & bb,
                Point< T,DIM > & ip)
{
  return detail::intersect_ray_bbox(R, bb, ip);
}

/*!
 * \brief Computes the intersection of the given segment, S, with the Box, bb.
 *     ip the point of intersection on S.
 * \return status true iff bb intersects with S, otherwise, false.
 *
 * Computes Segment Box intersection using the slab method from pg 180 of
 * Real Time Collision Detection by Christer Ericson.
 * WIP: More test cases for this
 */
template < typename T, int DIM >
bool intersect( const Segment< T,DIM > & S,
                const BoundingBox< T,DIM > & bb,
                Point< T,DIM > & ip)
{
  return detail::intersect_seg_bbox(S, bb, ip);
}

/*!
 * \brief Determines if two axis aligned bounding boxes intersect
 * \param [in] bb1 user-supplied axis aligned bounding box.
 * \param [in] bb2 user-supplied axis aligned bounding box.
 * \return true iff bb1 intersects with bb2, otherwise, false.
 */
template < typename T, int DIM >
bool intersect( const BoundingBox< T, DIM >& bb1,
                const BoundingBox< T, DIM >& bb2)
{
  return bb1.intersectsWith(bb2);
}

/*!
 * \brief Determines if a triangle and a bounding box intersect
 * \param [in] tri user-supplied triangle (with three vertices).
 * \param [in] bb user-supplied axis aligned bounding box.
 * \return true iff tri intersects with bb, otherwise, false.
 */
template < typename T >
bool intersect( const Triangle< T, 3 >& tri, const BoundingBox< T, 3 >& bb)
{
  return detail::intersect_tri_bbox(tri, bb);
}

/*!
 * \brief Determines if a 3D triangle intersects a 3D ray.
 * \param [in] tri A 3D triangle
 * \param [in] ray A 3D ray
 * \return true iff tri intersects with ray, otherwise, false.
 */
template < typename T >
bool intersect(const Triangle< T, 3 >& tri, const Ray< T,3 >& ray)
{
  T t = T();
  return detail::intersect_tri_ray(tri, ray, t);
}

/*!
 * \brief Determines if a 3D triangle intersects a 3D ray.
 * \param [in] tri A 3D triangle
 * \param [in] ray A 3D ray
 * \param [out] t Intersection point of tri and R, w.r.t. parametrization of R
 * \note If there is an intersection, the intersection point is:  R.at(t)
 * \return true iff tri intersects with ray, otherwise, false.
 */
template < typename T >
bool intersect(const Triangle< T, 3 >& tri, const Ray< T,3 >& ray, T& t)
{
  return detail::intersect_tri_ray(tri, ray, t);
}

/*!
 * \brief Determines if a 3D triangle intersects a 3D segment.
 * \param [in] tri A 3D triangle
 * \param [in] seg A 3D line segment
 * \return true iff tri intersects with seg, otherwise, false.
 */
template < typename T >
bool intersect(const Triangle< T, 3 >& tri, const Segment< T,3 >& seg)
{
  T t = T();
  return detail::intersect_tri_segment(tri, seg, t);
}

/*!
 * \brief Determines if a 3D triangle intersects a 3D segment.
 * \param [in] tri A 3D triangle
 * \param [in] seg A 3D line segment
 * \param [out] t Intersection point of tri and seg, w.r.t. seg's parametrization
 * \note If there is an intersection, the intersection point pt is:
 *                     pt = seg.source() + t * ( seg.dest() - seg.target() )
 * \return true iff tri intersects with seg, otherwise, false.
 */
template < typename T >
bool intersect(const Triangle< T, 3 >& tri, const Segment< T,3 >& seg, T& t)
{
  return detail::intersect_tri_segment(tri, seg, t);
}

template < typename T >
bool intersect(const OrientedBoundingBox< T, 1 > & b1,
  const OrientedBoundingBox< T, 1 >& b2)
{
  T c1 = b1.getCentroid().array()[0];
  T c2 = b2.getCentroid().array()[0];

  T e1 = b1.getExtents()[0];
  T e2 = b2.getExtents()[0];

  if (c1 + e1 > c2 - e2) return true;
  if (c2 + e2 > c1 - e1) return true;

  return false;
}

// static helper method
template < typename T >
static inline T abs(T x)
{
  if (x < T()) return -x;
  return x;
}

/*!
 *******************************************************************************
 * \brief Determines if a 2D OBB intersects a 2D OBB.
 * \param [in] b1 A 2D OrientedBoundingBox
 * \param [in] b2 A 2D OrientedBoundingBox
 * \return true iff b1 intersects with b2, otherwise, false.
 *******************************************************************************
 */
template < typename T >
bool intersect(const OrientedBoundingBox< T, 2 >& b1,
  const OrientedBoundingBox< T, 2 >& b2)
{
  Vector< T, 2 > c1(b1.getCentroid());
  Vector< T, 2 > c2(b2.getCentroid());

  Vector< T, 2 > e1 = b1.getExtents();
  Vector< T, 2 > e2 = b2.getExtents();

  const Vector< T, 2 > *u1 = b1.getAxes();
  const Vector< T, 2 > *u2 = b2.getAxes();

  Vector< T, 2 > d = c2 - c1;

  for (int i = 0; i < 2; i++) {
    if (abs< T >(d.dot(u1[i])) > e1[i] + abs< T >((e2[0]*u2[0]).dot(u1[i]))
        + abs< T >((e2[1]*u2[1]).dot(u1[i]))) {
      return false;
    }
  }

  for (int i = 0; i < 2; i++) {
    if (abs< T >(d.dot(u2[i])) > e2[i] + abs< T >((e1[0]*u1[0]).dot(u2[i]))
        + abs< T >((e1[1]*u1[1]).dot(u2[i]))) {
      return false;
    }
  }

  return true;
}

/*!
 *******************************************************************************
 * \brief Determines if a 3D OBB intersects a 3D OBB.
 * \param [in] b1 A 3D OrientedBoundingBox
 * \param [in] b2 A 3D OrientedBoundingBox
 * \return true iff b1 intersects with b2, otherwise, false.
 *******************************************************************************
 */
template < typename T >
bool intersect(const OrientedBoundingBox< T, 3 >& b1,
  const OrientedBoundingBox< T, 3 >& b2)
{
  static const double EPS = 1.0e-4;
  Vector< T, 3 > d = Vector< T, 3 >(b1.getCentroid())
    - Vector< T, 3 >(b2.getCentroid());

  Vector< T, 3 > e1 = b1.getExtents();
  Vector< T, 3 > e2 = b2.getExtents();

  const Vector< T, 3 > *u1 = b1.getAxes();
  const Vector< T, 3 > *u2 = b2.getAxes();

  // compute r and r^T here:
  Vector< T, 3 > r[3];
  Vector< T, 3 > rt[3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      r[i][j] = abs< T >(u1[i].dot(u2[j]));
      rt[i][j] = abs< T >(u1[j].dot(u2[i]));
    }
  }

  // check for separating planes parallel to faces
  for (int i = 0; i < 3; i++) {
    if (abs< T >(d.dot(u1[i])) > e1[i] + e2.dot(r[i]) + EPS) {
      return false;
    }
  }

  for (int i = 0; i < 3; i++) {
    if (abs< T >(d.dot(u2[i])) > e2[i] + e1.dot(rt[i]) + EPS) return false;
  }

  // check for separating planes with normals parallel to cross product of edges
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      T left = abs< T >(d.dot(u1[(i + 2) % 3])*r[(i + 1) % 3][j]
        - (d.dot(u1[(i + 1) % 3])*r[(i + 2) % 3][j]));
      T right = (e1[(i + 1) % 3]*r[(i + 2) % 3][j]
        + e1[(i + 2) % 3]*r[(i + 1) % 3][j]);
      right += (e2[(i + 1) % 3]*r[i][(j + 2) % 3]
        + e2[(i + 2) % 3]*r[i][(j + 1) % 3]);
      if (left > right + EPS) return false;
    }
  }

  // didn't find a separating anything
  return true;
}


} /* namespace primal */
} /* namespace axom */

#endif /* INTERSECTION_HPP_ */

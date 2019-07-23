// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file intersect.hpp
 *
 * \brief Consists of functions to test intersection among geometric primitives.
 */

#ifndef INTERSECTION_BEZIER_HPP_
#define INTERSECTION_BEZIER_HPP_

#include "axom/core/numerics/Determinants.hpp"
#include "axom/core/utilities/Utilities.hpp"

#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/OrientedBoundingBox.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Ray.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/Sphere.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"

#include "axom/primal/operators/squared_distance.hpp"
#include "axom/primal/operators/intersect.hpp"

namespace axom
{
namespace primal
{

template < typename T, int NDIMS>
bool intersect_bezier_helper( const BezierCurve< T, NDIMS>& c1,
                              const BezierCurve< T, NDIMS>& c2,
                              std::vector< T >& sp,
                              std::vector< T >& tp,
                              double sq_tol,
                              int order1, int order2,
                              double s_offset, double s_scale,
                              double t_offset, double t_scale)
{
  using BCurve = BezierCurve< T, NDIMS>;

  // Check bounding boxes to short-circuit the intersection
  if(!primal::intersect(c1.boundingBox(), c2.boundingBox()) )
  {
    return false;
  }

  //std::cout<< "(";
  bool foundIntersection = false;

  if ( c1.isLinear(sq_tol) && c2.isLinear(sq_tol))
  {
    T s,t;
    if(intersect_2d_linear(c1[0],c1[order1],c2[0],c2[order2],s,t))
    {
      sp.push_back(s_offset + s_scale * s);
      tp.push_back(t_offset + t_scale * t);
      foundIntersection = true;
      //std::cout << "*" ;
    }
  }
  else
  {
    constexpr double splitVal = 0.5;
    constexpr double scaleFac = 0.5;

    BCurve c3(order1);
    BCurve c4(order1);
    c1.split(splitVal,c3,c4);

    s_scale *= scaleFac;

    if (intersect_bezier_helper(c2,c3,tp,sp, sq_tol, order2, order1,
                                t_offset, t_scale, s_offset, s_scale))
    {
      foundIntersection = true;


    }
    if (intersect_bezier_helper(c2,c4,tp,sp, sq_tol, order2, order1,
                                t_offset, t_scale, s_offset + s_scale, s_scale))
    {
      foundIntersection = true;
    }
  }
  //std::cout<<")";

  return foundIntersection;

}

/*!
 * \brief Tests if Bezier Curves c1 and c2 intersect.
 * \return status true iff c1 intersects with c2, otherwise false.
 *
 * \param c1, c2 BezierCurve objects to intersect
 * \param sp, tp vector of type T parameter space intersection points (t-values
 * and s-values) for c1 and c2, respectively
 */
template < typename T, int NDIMS>
bool intersect_bezier( const BezierCurve< T, NDIMS>& c1,
                       const BezierCurve< T, NDIMS>& c2,
                       std::vector< T >& sp,
                       std::vector< T >& tp,
                       double sq_tol = 1E-16)
{
  const int ord1=c1.getOrder();
  const int ord2=c2.getOrder();
  const double offset = 0.;
  const double scale = 1.;

  return intersect_bezier_helper(c1,c2,sp,tp, sq_tol, ord1, ord2,
                                 offset, scale, offset, scale);
}

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
 * \note This function does not properly handle collinear lines
 */

template < typename T, int NDIMS>
bool intersect_2d_linear( const Point<T,NDIMS> &a, const Point<T,NDIMS> &b,
                          const Point<T,NDIMS> &c, const Point<T,NDIMS> &d,
                          T &s, T &t)
{
  // Implementation inspired by Section 5.1.9.1 of
  // C. Ericson's Real-Time Collision Detection book

  // Note: Uses exact floating point comparisons since the subdivision algorithm
  // provides both sides of the line segments for interior curve points.

  AXOM_STATIC_ASSERT( NDIMS==2);

  // compute signed areas of endpoints of segment (c,d) w.r.t. segment (a,b)
  auto area1 = detail::twoDcross(a,b,c);
  auto area2 = detail::twoDcross(a,b,d);

  // early return if both have same orientation, or if d is collinear w/ (a,b)
  if( area2 == 0. || (area1 * area2) > 0.)
    return false;

  // compute signed areas of endpoints of segment (a,b) w.r.t. segment (c,d)
  auto area3 = detail::twoDcross(c,d,a);
  auto area4 = area3+area1-area2; // equivalent to detail::twoDcross(c,d,b)

  // early return if both have same orientation, or if b is collinear w/ (c,d)
  if( area4 == 0. || (area3 * area4) > 0.)
    return false;

  // Compute intersection parameters using linear interpolation
  // Divisions are safe due to early return conditions
  s = area3 / (area3-area4);
  t = area1 / (area1-area2);
  return true;
}

} // namespace primal
} // namespace axom

#endif // PRIMAL_INTERSECTION_BEZIER_HPP_

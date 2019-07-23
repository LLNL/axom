// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file intersect_bezier.hpp
 *
 * \brief Consists of functions to test intersection of Bezier curves.
 */

#ifndef INTERSECTION_BEZIER_HPP_
#define INTERSECTION_BEZIER_HPP_

#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/operators/detail/intersect_bezier_impl.hpp"

#include <vector>

namespace axom
{
namespace primal
{

/*!
 * \brief Tests if two Bezier Curves \a c1 and \a c2 intersect.
 * \return status true iff \a c1 intersects \a c2, otherwise false.
 *
 * \param [in] c1 The first BezierCurve, parametrized in [0,1)
 * \param [in] c2 The second BezierCurve, parametrized in [0,1)
 * \param [out] sp vector of parameter space intersection points for \a c1
 * \param [out] tp vector of parameter space intersection points for \a c2
 * \param [in] sq_tol Tolerance parameter for determining if a curve can
 * be approximated by a line segment.
 * \return True if the curves intersect, false otherwise. Intersection
 * parameters are stored in \a sp and \a tp
 *
 * Finds all intersection points between the two curves.
 *
 * \note This function assumes that the curves are in general position.
 * Specifically, we assume that all intersections are at points and that
 * the curves don't overlap.
 *
 * \note This function assumes that the curves are half-open, i.e. they
 * contain their first endpoint, but not their last endpoint. Thus, the
 * curves do not intersect at \f$ s==1 \f$ or at \f$ t==1 \f$.
 */
template < typename T, int NDIMS>
bool intersect_bezier( const BezierCurve< T, NDIMS>& c1,
                       const BezierCurve< T, NDIMS>& c2,
                       std::vector< T >& sp,
                       std::vector< T >& tp,
                       double sq_tol = 1E-16)
{
  const double offset = 0.;
  const double scale = 1.;

  return detail::intersect_bezier_helper(c1, c2, sp, tp, sq_tol,
                                         c1.getOrder(), c2.getOrder(),
                                         offset, scale, offset, scale);
}

} // namespace primal
} // namespace axom

#endif // PRIMAL_INTERSECTION_BEZIER_HPP_

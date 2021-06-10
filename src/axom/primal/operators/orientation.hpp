// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file orientation.hpp
 *
 * \brief Consists of a set of templated (overloaded) routines used to calculate
 *  the orientation of a given point to another geometric entity.
 *
 */

#ifndef AXOM_PRIMAL_ORIENTATION_HPP_
#define AXOM_PRIMAL_ORIENTATION_HPP_

#include "axom/core/numerics/Determinants.hpp"
#include "axom/core/utilities/Utilities.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/OrientationResult.hpp"

#include "axom/slic/interface/slic.hpp"

namespace axom
{
namespace primal
{
/*!
 * \brief Computes the orientation of the given point, p, with respect to a
 *  supplied oriented triangle.
 * \param [in] p the query point.
 * \param [in] tri oriented triangle.
 * \return The orientation of the point with respect to the given triangle.
 * \note The triangle lies in a plane that divides space into the positive
 * half-space and the negative half-space.  The triangle's normal vector
 * points into the positive half-space.
 * The return value of this routine can be one of the following:
 * <ul>
 *  <li> ON_BOUNDARY, when the point is coplanar with the triangle </li>
 *  <li> ON_POSITIVE_SIDE, when the point lies in the positive half-space </li>
 *  <li> ON_NEGATIVE_SIDE, when the point lies in the negative half-space </li>
 * </ul>
 */
template <typename T>
inline int orientation(const Point<T, 3>& p, const Triangle<T, 3>& tri)
{
  // clang-format off
  double det = numerics::determinant( tri[0][0], tri[0][1], tri[0][2], 1.0,
                                      tri[1][0], tri[1][1], tri[1][2], 1.0,
                                      tri[2][0], tri[2][1], tri[2][2], 1.0,
                                      p[0],      p[1],      p[2], 1.0  );
  // clang-format on

  int orient = -1;

  if(axom::utilities::isNearlyEqual(det, 0.0, 1.0e-9))
  {
    orient = ON_BOUNDARY;
  }
  else if(det < 0.0f)
  {
    // outside
    orient = ON_POSITIVE_SIDE;
  }
  else
  {
    // inside
    orient = ON_NEGATIVE_SIDE;
  }

  return orient;
}

/*!
 * \brief Computes the orientation of the given point, p, with respect to a
 *  supplied oriented segment.
 * \param [in] p the query point.
 * \param [in] seg the user-supplied segment.
 * \return The orientation of the point with respect to the given segment.
 * \note The return value can be one of the following:
 * <ul>
 *  <li> ON_BOUNDARY, when the point is collinear with the points that define
 *       the segment.
 *  </li>
 *  <li> ON_POSITIVE_SIDE, when the point is clockwise, i.e., to the right of
 *       the directed segment.
 *  </li>
 *  <li> ON_NEGATIVE_SIDE, when the point is counter-clockwise, i.e., to the
 *       left of the directed segment.
 *  </li>
 * </ul>
 */
template <typename T>
inline int orientation(const Point<T, 2>& p, const Segment<T, 2>& seg)
{
  // clang-format off
  double det = numerics::determinant( seg.source()[0], seg.source()[1], 1.0,
                                      seg.target()[0], seg.target()[1], 1.0,
                                      p[0],            p[1], 1.0  );
  // clang-format on

  int orient = -1;

  if(axom::utilities::isNearlyEqual(det, 0.0))
  {
    // collinear
    orient = ON_BOUNDARY;
  }
  else if(det < 0.0f)
  {
    // outside, clockwise, to the right
    orient = ON_POSITIVE_SIDE;
  }
  else
  {
    // inside, counter-clockwise, to the left
    orient = ON_NEGATIVE_SIDE;
  }

  return orient;
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_ORIENTATION_HPP_

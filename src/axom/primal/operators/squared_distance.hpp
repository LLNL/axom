// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file squared_distance.hpp
 *
 * \brief Consists of a set of templated (overloaded) routines used to calculate
 *  the "signed" squared distance between two geometric entities.
 */

#ifndef AXOM_PRIMAL_SQUAREDDISTANCE_HPP_
#define AXOM_PRIMAL_SQUAREDDISTANCE_HPP_

#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/operators/closest_point.hpp"

#include "axom/slic/interface/slic.hpp"

namespace axom
{
namespace primal
{
/*!
 * \brief Computes the squared distance from point A to point B,
 *  represented by arrays of length N.
 * \param [in] A source point
 * \param [in] B end point.
 * \param [in] N length of A and B.
 * \return d the distance from point A to point B.  If N < 1, return 0.
 * \pre A and B are arrays of at least length N.
 */
inline double squared_distance(const double* A, const double* B, int N)
{
  double retval = 0;

  for(int i = 0; i < N; ++i)
  {
    double d = B[i] - A[i];
    retval += d * d;
  }

  return retval;
}

/*!
 * \brief Computes the squared distance from point A to point B.
 * \param [in] A source point
 * \param [in] B end point.
 * \return d the distance from point A to point B.
 */
template <typename T, int NDIMS>
inline double squared_distance(const Point<T, NDIMS>& A, const Point<T, NDIMS>& B)
{
  Vector<T, NDIMS> v(A, B);
  return (v.squared_norm());
}

/*!
 * \brief Computes the minimum squared distance from a query point, P, to a
 *  given axis-aligned bounding box B.
 * \param [in] P the query point.
 * \param [in] B the axis-aligned bounding box.
 * \return d the signed distance from P to the closest point on B.
 */
template <typename T, int NDIMS>
inline double squared_distance(const Point<T, NDIMS>& P,
                               const BoundingBox<T, NDIMS>& B)
{
  if(B.contains(P))
  {
    /* short-circuit */
    return 0.0f;
  }

  // compute closest point to the box
  Point<T, NDIMS> cp;
  for(int i = 0; i < NDIMS; ++i)
  {
    cp[i] = P[i];
    if(cp[i] < B.getMin()[i])
    {
      cp[i] = B.getMin()[i];
    }

    if(cp[i] > B.getMax()[i])
    {
      cp[i] = B.getMax()[i];
    }
  }

  // return squared signed distance to the closest point
  return squared_distance(P, cp);
}

/*!
 * \brief Computes the minimum squared distance from a query point, P, to the
 *  given segment, S.
 * \param [in] P the query point.
 * \param [in] S the input segment.
 * \return d the minimum distance from P on the
 */
template <typename T, int NDIMS>
inline double squared_distance(const Point<T, NDIMS>& P,
                               const Segment<T, NDIMS>& S)
{
  Vector<T, NDIMS> ab(S.source(), S.target());
  Vector<T, NDIMS> ac(S.source(), P);

  const T e = Vector<T, NDIMS>::dot_product(ac, ab);

  // outside segment, on the side of a
  // Testing if closest point is A
  if(e <= 0.0f)
  {
    return ac.squared_norm();
  }

  // outside segment, on the side of b
  // Testing if closest point is B
  const T f = ab.squared_norm();
  if(e >= f)
  {
    Vector<T, NDIMS> bc(S.target(), P);
    return bc.squared_norm();
  }

  // P projects onto the segment
  // Otherwise, we are in between A,B, therefore we project inside A,B.
  const T dist = ac.squared_norm() - (e * e / f);
  return dist;
}

/*!
 * \brief Computes the minimum squared distance from a query point, P, to the
 *  closest point on the given triangle.
 * \param [in] P the query point.
 * \param [in] tri the supplied triangle.
 * \return d the distance from Q to the closest point on the triangle T.
 */
template <typename T, int NDIMS>
inline double squared_distance(const Point<T, NDIMS>& P,
                               const Triangle<T, NDIMS>& tri)
{
  Point<T, NDIMS> cpt = closest_point(P, tri);
  return squared_distance(P, cpt);
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_SQUAREDDISTANCE_HPP_

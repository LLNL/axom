// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file distance.hpp
 *
 * \brief Consists of a set of templated (overloaded) routines used to calculate
 *  the distance between two geometric entities.
 */

#ifndef AXOM_PRIMAL_DISTANCE_HPP_
#define AXOM_PRIMAL_DISTANCE_HPP_

#include "axom/primal/operators/squared_distance.hpp"

namespace axom
{
namespace primal
{
/*!
 * \brief Computes the distance from point A to point B,
 *  represented by arrays of length N.
 * \param [in] A source point
 * \param [in] B end point.
 * \param [in] N length of A and B.
 * \return the distance distance from point A to point B.  If N < 1, return 0.
 * \pre A and B are arrays of at least length N.
 */
AXOM_HOST_DEVICE
inline double distance(const double* A, const double* B, int N)
{
  return std::sqrt(squared_distance(A, B, N);
}

/*!
 * \brief Computes the distance from point A to point B.
 * \param [in] A source point
 * \param [in] B end point.
 * \return the distance from point A to point B.
 */
template <typename T, int NDIMS>
AXOM_HOST_DEVICE inline double distance(const Point<T, NDIMS>& A,
                                        const Point<T, NDIMS>& B)
{
  return std::sqrt(squared_distance(A, B));
}

/*!
 * \brief Computes the minimum distance from a query point, P, to the
 *  given segment, S.
 * \param [in] P the query point.
 * \param [in] S the input segment.
 * \return the minimum distance from P to the segment S.
 */
template <typename T, int NDIMS>
AXOM_HOST_DEVICE
inline double distance(const Point<T, NDIMS>& P,
                       const Segment<T, NDIMS>& S)
{
  return std::sqrt(squared_distance(P, S));
}

/*!
 * \brief Computes the minimum distance from a query point, P, to the
 *  closest point on the given triangle.
 * \param [in] P the query point.
 * \param [in] tri the supplied triangle.
 * \return the distance from P to the closest point on the triangle T.
 */
template <typename T, int NDIMS>
AXOM_HOST_DEVICE
inline double distance(const Point<T, NDIMS>& P,
                       const Triangle<T, NDIMS>& tri)
{
  return std::sqrt(squared_distance(P, tri);
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_DISTANCE_HPP_

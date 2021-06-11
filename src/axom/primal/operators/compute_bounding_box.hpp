// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file compute_bounding_box.hpp
 *
 * \brief Consists of functions to create bounding boxes.
 */

#ifndef AXOM_PRIMAL_COMPUTE_BOUNDING_BOX_HPP_
#define AXOM_PRIMAL_COMPUTE_BOUNDING_BOX_HPP_

#include "axom/primal/geometry/NumericArray.hpp"  // for numeric arrays
#include "axom/core/numerics/Matrix.hpp"          // for Matrix
#include "axom/core/Macros.hpp"                   // for AXOM_HOST__DEVICE
#include "axom/core/numerics/eigen_solve.hpp"     // for eigen_solve
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/OrientedBoundingBox.hpp"

namespace axom
{
namespace primal
{
/*!
 * \brief Creates a bounding box which contains the collection of passed in
 * points.
 *
 * Call srand() to initialize the random number generator before using this
 * function.
 *
 * \param [in] pts array of points
 * \param [in] n number of points
 * \note if n <= 0, invokes default constructor
 */
template <typename T, int NDIMS>
OrientedBoundingBox<T, NDIMS> compute_oriented_bounding_box(
  const Point<T, NDIMS> *pts,
  int n)
{
  return OrientedBoundingBox<T, NDIMS>(pts, n);
}

/*!
 * \brief Creates an oriented bounding box which contains the passed in OBBs.
 *
 * Call srand() to initialize the random number generator before using this
 * function.
 *
 * \param [in] l left obb
 * \param [in] r right obb
 */
template <typename T, int NDIMS>
OrientedBoundingBox<T, NDIMS> merge_boxes(const OrientedBoundingBox<T, NDIMS> &l,
                                          const OrientedBoundingBox<T, NDIMS> &r)
{
  // TODO: See if this initial check can be improved so it's not so costly in
  // cases where it doesn't end up helping
  if(l.contains(r))
  {
    return l;
  }
  if(r.contains(l))
  {
    return r;
  }

  std::vector<Point<T, NDIMS>> lv = l.vertices();
  std::vector<Point<T, NDIMS>> rv = r.vertices();

  const int size = (1 << NDIMS);

  Point<T, NDIMS> pts[2 * size];

  for(int i = 0; i < size; i++)
  {
    pts[i] = lv[i];
    pts[i + size] = rv[i];
  }

  return compute_oriented_bounding_box<T, NDIMS>(pts, 2 * size);
}

/*!
 * \brief Constructor. Creates a bounding box which contains the passed in
 * bounding boxes.
 *
 * \param [in] l left bb
 * \param [in] r right bb
 */

template <typename T, int NDIMS>
BoundingBox<T, NDIMS> merge_boxes(const BoundingBox<T, NDIMS> &l,
                                  const BoundingBox<T, NDIMS> &r)
{
  BoundingBox<T, NDIMS> res(l);
  res.addBox(r);
  return res;
}

/*!
 * \brief Creates a bounding box around a Triangle
 *
 * \param [in] tri The Triangle
 */
template <typename T, int NDIMS>
AXOM_HOST_DEVICE BoundingBox<T, NDIMS> compute_bounding_box(
  const Triangle<T, NDIMS> &tri)
{
  BoundingBox<T, NDIMS> res(tri[0]);
  for(int i = 1; i < 3; i++)
  {
    res.addPoint(tri[i]);
  }
  return res;
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_COMPUTE_BOUNDING_BOX_HPP_

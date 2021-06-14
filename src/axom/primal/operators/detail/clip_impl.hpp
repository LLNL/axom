// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file clip_impl.hpp
 *
 * \brief Helper functions for the primal clipping operators
 */

#ifndef AXOM_PRIMAL_CLIP_IMPL_HPP_
#define AXOM_PRIMAL_CLIP_IMPL_HPP_

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/Polygon.hpp"

#include "axom/primal/operators/orientation.hpp"

namespace axom
{
namespace primal
{
namespace detail
{
/*! Returns true when index is even */
inline bool isEven(int index) { return (index & 1) == 0; }

/*!
 * \brief Specialized point plane classifier for axis aligned planes
 *
 * \param [in] pt The plane to classify
 * \param [in] index The index of the axis aligned plane. See below for mapping
 * \param [in] val The plane's coordinate with respect to the given axis
 * \param [in] eps A parameter for thickening width of the plane (default 1e-8)
 *
 * Mapping of index to axis
 * * 0 -> -x axis
 * * 1 -> +x axis
 * * 2 -> -y axis
 * * 3 -> +y axis
 * * 4 -> -z axis
 * * 5 -> +z axis
 *
 * \return An OrientedSide (ON_POSITIVE_SIDE, ON_NEGATIVE_SIDE, ON_BOUNDARY)
 *         value based on the relative orientations of point pt and the
 *         corresponding plane associated with index.
 * \see OrientedPlane enum
 */
template <typename T, int NDIMS>
int classifyPointAxisPlane(const Point<T, NDIMS>& pt,
                           int index,
                           T val,
                           const double eps = 1e-8)
{
  // Note: we are exploiting the fact that the planes are axis aligned
  // So the dot product is +/- the given coordinate.
  // In general, we would need to call distance(pt, plane) here
  T dist = isEven(index) ? val - pt[index / 2] : pt[index / 2] - val;

  if(dist > eps)
  {
    return ON_POSITIVE_SIDE;
  }
  if(dist < -eps)
  {
    return ON_NEGATIVE_SIDE;
  }

  return ON_BOUNDARY;
}

/*!
 * \brief Finds the clipping intersection point between points a and b.
 *
 * \param [in] a The point behind the plane
 * \param [in] b The point in front of the plane
 * \param [in] index The index of the axis aligned plane.
 * \param [in] val The plane's coordinate with respect to the given axis
 * \return The point between a and b whose corresponding coordinate is val
 *
 * \see classifyPointAxisPlane for description of how index maps to coordinates.
 */
template <typename T, int NDIMS>
Point<T, NDIMS> findIntersectionPoint(const Point<T, NDIMS>& a,
                                      const Point<T, NDIMS>& b,
                                      int index,
                                      T val)
{
  typedef Point<T, NDIMS> PointType;

  // Need to find a parameter t for the point pt, such that,
  // * 0 <= t <= 1
  // * pt = a + t * (b-a)
  // * pt[ index/2]  == val

  T t = (val - a[index / 2]) / (b[index / 2] - a[index / 2]);
  SLIC_ASSERT(0. <= t && t <= 1.);

  PointType ret = PointType(a.array() + t * (b.array() - a.array()));
  SLIC_ASSERT(classifyPointAxisPlane(ret, index, val) == ON_BOUNDARY);

  return ret;
}

/*!
 * \brief Clips the vertices of the polygon to be behind the plane.
 *
 * This is a specialization of the Sutherland-Hodgeman clipping algorithm
 * for axis-aligned planes
 *
 * \param [in] prevPoly  An input polygon with the vertices to clip
 * \param [out] currentPoly An output polygon whose coordinates are clipped
 *                          against this plane.
 * \param [in] index The index of the axis aligned plane.
 * \param [in] val The plane's coordinate with respect to the given axis
 *
 * \note Algorithm for robust clipping against "thick" planes derived from
 *       Section 8.3 of Christer Ericson's "Real-Time Collision Detection"
 *       and is based on the Sutherland-Hodgeman clipping algorithm.
 *       We are only keeping the "back" polygon, w.r.t. that algorithm.
 * \see classifyPointAxisPlane for description of how index maps to coordinates.
 */
template <typename T, int NDIMS>
void clipAxisPlane(const Polygon<T, NDIMS>* prevPoly,
                   Polygon<T, NDIMS>* currentPoly,
                   int index,
                   T val)
{
  typedef Point<T, NDIMS> PointType;

  currentPoly->clear();
  int numVerts = prevPoly->numVertices();

  if(numVerts == 0)
  {
    return;
  }

  // Initialize point a with the last vertex of the polygon
  const PointType* a = &(*prevPoly)[numVerts - 1];
  int aSide = classifyPointAxisPlane(*a, index, val);

  for(int i = 0; i < numVerts; ++i)
  {
    const PointType* b = &(*prevPoly)[i];
    int bSide = classifyPointAxisPlane(*b, index, val);

    switch(bSide)
    {
    case ON_POSITIVE_SIDE:
      if(aSide == ON_NEGATIVE_SIDE)
      {
        currentPoly->addVertex(findIntersectionPoint(*a, *b, index, val));
      }
      break;
    case ON_BOUNDARY:
      if(aSide == ON_NEGATIVE_SIDE)
      {
        currentPoly->addVertex(*b);
      }
      break;
    case ON_NEGATIVE_SIDE:
      switch(aSide)
      {
      case ON_POSITIVE_SIDE:
        currentPoly->addVertex(findIntersectionPoint(*a, *b, index, val));
        currentPoly->addVertex(*b);
        break;
      case ON_BOUNDARY:
        currentPoly->addVertex(*a);
        currentPoly->addVertex(*b);
        break;
      case ON_NEGATIVE_SIDE:
        currentPoly->addVertex(*b);
        break;
      }
      break;
    }

    // swap a and b
    a = b;
    aSide = bSide;
  }
}

}  // namespace detail
}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_CLIP_IMPL_HPP_

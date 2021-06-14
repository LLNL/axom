// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file clip.hpp
 *
 * \brief Consists of a set of methods to clip a primal primitive against
 *        another primal primitive
 */

#ifndef AXOM_PRIMAL_CLIP_HPP_
#define AXOM_PRIMAL_CLIP_HPP_

#include "axom/core/utilities/Utilities.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/Polygon.hpp"

#include "axom/primal/operators/detail/clip_impl.hpp"

namespace axom
{
namespace primal
{
/*!
 * \brief Clips a 3D triangle against an axis-aligned bounding box in 3D
 *
 * \param [in] tri The triangle to clip
 * \param [in] bbox The bounding box to clip against
 * \return A planar polygon of the triangle clipped against the bounding box.
 *         If the triangle is completely outside the bounding box,
 *         the returned polygon is empty (i.e. it has no vertices).
 *
 * \note Using a specialization of the Sutherland-Hodgeman clipping algorithm
 *       for axis aligned planes
 */
template <typename T>
Polygon<T, 3> clip(const Triangle<T, 3>& tri, const BoundingBox<T, 3>& bbox)
{
  using BoundingBoxType = BoundingBox<T, 3>;
  using PolygonType = Polygon<T, 3>;

  // Use two polygons with pointers for 'back-buffer'-like swapping
  const int MAX_VERTS = 6;
  PolygonType poly[2] = {PolygonType(MAX_VERTS), PolygonType(MAX_VERTS)};
  PolygonType* currentPoly = &poly[0];
  PolygonType* prevPoly = &poly[1];

  // First check if the triangle is contained in the bbox, if not we are empty
  BoundingBoxType triBox;
  triBox.addPoint(tri[0]);
  triBox.addPoint(tri[1]);
  triBox.addPoint(tri[2]);

  if(!bbox.intersectsWith(triBox))
  {
    return *currentPoly;  // note: currentPoly is empty
  }

  // Set up the initial polygon
  currentPoly->addVertex(tri[0]);
  currentPoly->addVertex(tri[1]);
  currentPoly->addVertex(tri[2]);

  // If all the vertices are contained, we have no work to do
  if(bbox.contains(triBox))
  {
    return *currentPoly;  // Note current poly has the same verts as tri
  }

  // Loop through the planes of the bbox and clip the vertices
  for(int dim = 0; dim < 3; ++dim)
  {
    // Optimization note: we should be able to save some work based on
    // the clipping plane and the triangle's bounding box

    if(triBox.getMax()[dim] > bbox.getMin()[dim])
    {
      axom::utilities::swap(prevPoly, currentPoly);
      detail::clipAxisPlane(prevPoly, currentPoly, 2 * dim + 0, bbox.getMin()[dim]);
    }

    if(triBox.getMin()[dim] < bbox.getMax()[dim])
    {
      axom::utilities::swap(prevPoly, currentPoly);
      detail::clipAxisPlane(prevPoly, currentPoly, 2 * dim + 1, bbox.getMax()[dim]);
    }
  }

  return *currentPoly;
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_CLIP_HPP_

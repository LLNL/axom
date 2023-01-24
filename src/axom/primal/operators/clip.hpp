// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
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
#include "axom/primal/geometry/Polyhedron.hpp"
#include "axom/primal/geometry/Octahedron.hpp"
#include "axom/primal/geometry/Plane.hpp"
#include "axom/primal/geometry/Tetrahedron.hpp"

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

/*!
 * \brief Clips a 3D octahedron against a tetrahedron in 3D, returning
 *        the geometric intersection of the octahedron and the tetrahedron
 *        as a polyhedron
 *
 *  This function clips the octahedron by the 4 planes obtained from the
 *  tetrahedron's faces (normals point inward). Clipping the
 *  octahedron/polyhedron by each plane gives the polyhedron above that plane.
 *  Clipping the polyhedron by a plane involves
 *  finding new vertices at the intersection of the polyhedron edges and
 *  the plane, removing vertices from the polyhedron that are below the
 *  plane, and redefining the neighbors for each vertex (a vertex is a
 *  neighbor of another vertex if there is an edge between them).
 *
 *
 * \param [in] oct The octahedron to clip
 * \param [in] tet The tetrahedron to clip against
 * \param [in] eps The epsilon value
 * \return A polyhedron of the octahedron clipped against the tetrahedron.
 *
 * \note Function is based off clipPolyhedron() in Mike Owen's PolyClipper.
 *
 */
template <typename T>
AXOM_HOST_DEVICE Polyhedron<T, 3> clip(const Octahedron<T, 3>& oct,
                                       const Tetrahedron<T, 3>& tet,
                                       double eps = 1.e-10)
{
  return detail::clipOctahedron(oct, tet, eps);
}

/*!
 * \brief Clips a 3D tetrahedron against another tetrahedron in 3D, returning
 *        the geometric intersection as a polyhedron
 *
 *  This function clips the first tetrahedron by the 4 planes obtained from the
 *  second tetrahedron's faces (normals point inward). Clipping the
 *  tetrahedron/polyhedron by each plane gives the polyhedron above that plane.
 *  Clipping the polyhedron by a plane involves
 *  finding new vertices at the intersection of the polyhedron edges and
 *  the plane, removing vertices from the polyhedron that are below the
 *  plane, and redefining the neighbors for each vertex (a vertex is a
 *  neighbor of another vertex if there is an edge between them).
 *
 *
 * \param [in] tet1 The tetrahedron to clip
 * \param [in] tet2 The tetrahedron to clip against
 * \param [in] eps The epsilon value
 * \return A polyhedron of the tetrahedron clipped against
 *         the other tetrahedron.
 *
 * \note Function is based off clipPolyhedron() in Mike Owen's PolyClipper.
 *
 */
template <typename T>
AXOM_HOST_DEVICE Polyhedron<T, 3> clip(const Tetrahedron<T, 3>& tet1,
                                       const Tetrahedron<T, 3>& tet2,
                                       double eps = 1.e-10)
{
  return detail::clipTetrahedron(tet1, tet2, eps);
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_CLIP_HPP_

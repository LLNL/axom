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
#include "axom/primal/geometry/Polyhedron.hpp"
#include "axom/primal/geometry/Octahedron.hpp"
#include "axom/primal/geometry/Plane.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/Tetrahedron.hpp"

#include "axom/primal/operators/intersect.hpp"
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
  using PointType = Point<T, NDIMS>;

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
  using PointType = Point<T, NDIMS>;

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

template <typename T, int NDIMS>
AXOM_HOST_DEVICE void poly_clip_vertices(Polyhedron<T, NDIMS>& poly,
                                         const Plane<T, NDIMS>& plane,
                                         const double eps,
                                         unsigned int& out_clipped)
{
  using PointType = Point<T, NDIMS>;
  using SegmentType = Segment<T, NDIMS>;

  // Loop over Polyhedron vertices
  int numVerts = poly.numVertices();
  for(axom::int8 i = 0; i < numVerts; i++)
  {
    int orientation = plane.getOrientation(poly[i], eps);

    // Vertex is under the plane
    if(orientation == ON_NEGATIVE_SIDE)
    {
      // Mark this vertex for removal later
      out_clipped |= 1 << i;

      // Check neighbors for vertex above the plane (edge clipped by plane)
      int numNeighbors = poly.getNumNeighbors(i);
      for(int j = 0; j < numNeighbors; j++)
      {
        axom::int8 neighborIndex = poly.getNeighbors(i)[j];

        int neighborOrientation = plane.getOrientation(poly[neighborIndex], eps);

        // Insert new vertex to polyhedron, where edge intersects plane.
        if(neighborOrientation == ON_POSITIVE_SIDE)
        {
          int expectedVertexIndex = poly.numVertices();

          T lerp_val;
          SegmentType seg(poly[i], poly[neighborIndex]);
          intersect(plane, seg, lerp_val);

          int newVertexIndex = poly.addVertex(seg.at(lerp_val));
          SLIC_ASSERT(newVertexIndex == expectedVertexIndex);

          poly.addNeighbors(newVertexIndex, {i, neighborIndex});

          // Update current vertex's & neighbor's neighbors with the
          // new vertex
          poly.getNeighbors(i)[j] = newVertexIndex;
          for(unsigned int k = 0; k < poly.getNumNeighbors(neighborIndex); k++)
          {
            if(poly.getNeighbors(neighborIndex)[k] == i)
            {
              poly.getNeighbors(neighborIndex)[k] = newVertexIndex;
            }
          }
        }
      }
    }
  }  // end of loop over Polyhedron vertices
}

template <typename T, int NDIMS>
AXOM_HOST_DEVICE void poly_clip_fix_nbrs(Polyhedron<T, NDIMS>& poly,
                                         const Plane<T, NDIMS>& plane,
                                         const int oldVerts,
                                         const double eps,
                                         const unsigned int clipped)
{
  NeighborCollection& poly_nbrs = poly.getNeighbors();
  // Keep copy of old connectivity
  NeighborCollection old_nbrs = poly.getNeighbors();
  for(int i = 0; i < poly.numVertices(); i++)
  {
    // Check clipped created vertices first, then vertices on the plane
    int vIndex = (i + oldVerts) % poly.numVertices();
    int vOrientation = plane.getOrientation(poly[vIndex], eps);

    if(vIndex >= oldVerts || vOrientation == ON_BOUNDARY)
    {
      for(unsigned int j = 0; j < poly_nbrs.getNumNeighbors(vIndex); j++)
      {
        int neighborIndex = poly_nbrs[vIndex][j];
        int neighborOrientation = plane.getOrientation(poly[neighborIndex], eps);

        // This neighbor is below the plane
        if(neighborOrientation == ON_NEGATIVE_SIDE)
        {
          // Look for 1st vertex along this face not below the plane.
          int iprev = vIndex;
          int inext = neighborIndex;
          int itmp = inext;

          int val = 0;

          while((clipped & (1 << inext)) && (val++ < poly.numVertices()))
          {
            itmp = inext;
            unsigned int next_nbrs = poly_nbrs.getNumNeighbors(inext);
            // Find next vertex along face
            for(unsigned int ni = 0; ni < next_nbrs; ni++)
            {
              if(poly_nbrs[inext][ni] == iprev)
              {
                inext = (ni == 0) ? poly_nbrs[inext][next_nbrs - 1]
                                  : poly_nbrs[inext][ni - 1];
                break;
              }
            }

            iprev = itmp;
          }

          // Remove neighbor from list if vertex found was already a neighbor or
          // is the vertex we are currently checking for.
          if(poly_nbrs[vIndex][(j + 1) % poly_nbrs.getNumNeighbors(vIndex)] ==
               inext ||
             inext == vIndex)
          {
            poly_nbrs[vIndex][j] = -1;
          }

          // Otherwise update neighbor lists of vertex found and vertex we are checking for.
          else
          {
            poly_nbrs[vIndex][j] = inext;

            if((clipped & (1 << inext)))
            {
              poly_nbrs.insertNeighborAtPos(inext, vIndex, 0);
              old_nbrs.insertNeighborAtPos(inext, -1, 0);
            }
            else
            {
              int offset;
              for(int oi = 0; oi < old_nbrs.getNumNeighbors(inext); oi++)
              {
                if(old_nbrs[inext][oi] == iprev)
                {
                  offset = oi;
                  break;
                }

                // Max offset
                if(oi == old_nbrs.getNumNeighbors(inext) - 1)
                {
                  offset = old_nbrs.getNumNeighbors(inext);
                }
              }
              poly_nbrs.insertNeighborAtPos(inext, vIndex, offset);
              old_nbrs.insertNeighborAtPos(inext, vIndex, offset);
            }
          }
        }
      }  // end of loop over vertex neighbors
    }
  }  // end of loop over Polyhedron vertices
  poly.getNeighbors().pruneNeighbors();
}

template <typename T, int NDIMS>
AXOM_HOST_DEVICE void poly_clip_reindex(Polyhedron<T, NDIMS>& poly,
                                        const unsigned int clipped)
{
  // Dictionary for old indices to new indices positions
  axom::int8 newIndices[Polyhedron<T, NDIMS>::MAX_VERTS];

  Polyhedron<T, NDIMS> old_poly = poly;

  poly.clear();

  int curIndex = 0;
  for(int i = 0; i < old_poly.numVertices(); i++)
  {
    if(!(clipped & (1 << i)))
    {
      // Non-clipped vertex
      newIndices[i] = curIndex++;
      poly.addVertex(old_poly[i]);
    }
  }

  // Reinsert neighbors into polyhedron
  for(int i = 0; i < old_poly.numVertices(); i++)
  {
    if(!(clipped & (1 << i)))
    {
      for(int j = 0; j < old_poly.getNumNeighbors(i); j++)
      {
        poly.addNeighbors(newIndices[i],
                          {newIndices[old_poly.getNeighbors()[i][j]]});
      }
    }
  }
}

/*!
 * \brief Finds the clipped intersection Polyhedron between Octahedron
 *        oct and Tetrahedron tet.
 *
 * \param [in] oct The octahedron
 * \param [in] tet The tetrahedron
 * \param [in] eps The tolerance for plane point orientation.
 * \return The Polyhedron formed from clipping the octahedron with a tetrahedron.
 *
 */
template <typename T, int NDIMS>
AXOM_HOST_DEVICE Polyhedron<T, NDIMS> clipOctahedron(
  const Octahedron<T, NDIMS>& oct,
  const Tetrahedron<T, NDIMS>& tet,
  double eps = 1.e-24)
{
  using PointType = Point<T, NDIMS>;
  using BoxType = BoundingBox<T, NDIMS>;
  using PlaneType = Plane<T, NDIMS>;

  // Initialize our polyhedron to return
  Polyhedron<T, NDIMS> poly;

  poly.addVertex(oct[0]);
  poly.addVertex(oct[1]);
  poly.addVertex(oct[2]);
  poly.addVertex(oct[3]);
  poly.addVertex(oct[4]);
  poly.addVertex(oct[5]);

  poly.addNeighbors(0, {1, 5, 4, 2});
  poly.addNeighbors(1, {0, 2, 3, 5});
  poly.addNeighbors(2, {0, 4, 3, 1});
  poly.addNeighbors(3, {1, 2, 4, 5});
  poly.addNeighbors(4, {0, 5, 3, 2});
  poly.addNeighbors(5, {0, 1, 3, 4});

  //Bounding Box of Polyhedron
  BoxType polyBox(&oct[0], 6);

  // Initialize planes from tetrahedron vertices
  // (Ordering here matters to get the correct winding)
  PlaneType planes[4] = {make_plane(tet[1], tet[3], tet[2]),
                         make_plane(tet[0], tet[2], tet[3]),
                         make_plane(tet[0], tet[3], tet[1]),
                         make_plane(tet[0], tet[1], tet[2])};

  //Clip octahedron by each plane
  for(int planeIndex = 0; planeIndex < 4; planeIndex++)
  {
    // Each bit value indicates if that Polyhedron vertex is formed from
    // Octahedron clipping with a plane.
    unsigned int clipped = 0;

    PlaneType plane = planes[planeIndex];

    // Check that plane intersects Polyhedron
    if(intersect(plane, polyBox))
    {
      int numVerts = poly.numVertices();

      // Clip polyhedron against current plane, generating extra vertices
      // where edges meet the plane.
      poly_clip_vertices(poly, plane, eps, clipped);

      // Adjust connectivity to link up newly-generated vertices.
      poly_clip_fix_nbrs(poly, plane, numVerts, eps, clipped);

      // Reindex polyhedron connectivity by removing vertices on the negative
      // side of the plane.
      poly_clip_reindex(poly, clipped);

      // Generate new bounding box for polyhedron
      polyBox = BoxType();

      for(int i = 0; i < poly.numVertices(); i++)
      {
        polyBox.addPoint(PointType(poly[i]));
      }
    }
  }

  return poly;
}

}  // namespace detail
}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_CLIP_IMPL_HPP_

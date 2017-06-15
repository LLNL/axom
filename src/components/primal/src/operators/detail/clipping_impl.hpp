/*
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

/*!
 *******************************************************************************
 * \file clipping_impl.hpp
 *
 * \brief Helper functions for the primal clipping operators
 *******************************************************************************
 */

#ifndef PRIMAL_CLIPPING_IMPL_HPP_
#define PRIMAL_CLIPPING_IMPL_HPP_

#include "primal/Point.hpp"
#include "primal/Triangle.hpp"
#include "primal/BoundingBox.hpp"
#include "primal/Polygon.hpp"

namespace axom {
namespace primal {

namespace detail {

/** Enum types for classifying points with respect to a thick plane */
struct PtPlaneClassifier {
  enum {
    POINT_ON_PLANE,
    POINT_IN_FRONT_OF_PLANE,
    POINT_BEHIND_PLANE
   };
};


/**
 * \brief Specialized point plane classifier for axis aligned planes
 *
 * \param pt The plane to classify
 * \param index The index of the axis aligned plane. See below for mapping
 * \param val The plane's coordinate with respect to the given axis
 *
 * Mapping of index to axis
 * * 0 -> -x axis
 * * 1 -> +x axis
 * * 2 -> -y axis
 * * 3 -> +y axis
 * * 4 -> -z axis
 * * 5 -> +z axis
 *
 * \return A PtPlaneClassifier value based on the relative orientations
 */
template<typename T, int NDIMS>
int classifyPointAxisPlane(const Point<T, NDIMS>& pt, int index, T val)
{
  const double EPS = 1e-8;

  // Note: we are exploiting the fact that the planes are axis aligned
  // So the dot product is +/- the given coordinate.
  // In general, we would need to call distance(pt, plane) here
  T dist = (index  % 2 == 0)
        ? -(pt[ index/2 ] - val)
        :  pt[ index/2 ] - val;

  if(dist > EPS)
    return PtPlaneClassifier::POINT_IN_FRONT_OF_PLANE;
  if(dist < -EPS)
    return PtPlaneClassifier::POINT_BEHIND_PLANE;

  return PtPlaneClassifier::POINT_ON_PLANE;
}

/**
 * \brief Finds the clipping intersection point between points a and b.
 *
 * \param a The point behind the plane
 * \param b The point in front of the plane
 * \param index The index of the axis aligned plane.
 * \param val The plane's coordinate with respect to the given axis
 * \return The point between a and b whose corresponding coordinate is val
 *
 * \see classifyPointAxisPlane for a description of how index maps to coordinates.
 */
template<typename T, int NDIMS>
Point<T,NDIMS> findIntersectionPoint(const Point<T, NDIMS>& a, const Point<T, NDIMS>& b, int index, T val)
{
  typedef Point<T,NDIMS> PointType;

  // Need to find a parameter t for the point pt, such that,
  // * 0 <= t <= 1
  // * pt = a + t (b-a)
  // * pt[ index/2]  == val

  T t = (val - a[ index /2 ]) / (b[index/2]- a[index/2]);
  SLIC_ASSERT(0. <= t && t <= 1.);

  PointType ret = PointType(a.array() + t * (b.array()-a.array()) );
  SLIC_ASSERT( classifyPointAxisPlane( ret, index, val) == PtPlaneClassifier::POINT_ON_PLANE);

  return ret;
}


/**
 * \brief Clips the vertices of the polygon to be behind the plane.
 *
 * This is a specialization of the Sutherland-Hodgeman clipping algorithm
 * for axis-aligned planes
 *
 * \param [in] prevPoly  An input polygon with the vertices to clip
 * \param [out] currentPoly An output polygon whose coordinates are clipped
 *                          against this plane.
 * \param index The index of the axis aligned plane.
 * \param val The plane's coordinate with respect to the given axis
 *
 * \note Algorithm for robust clipping against "thick" planes derived from
 *       Section 8.3 of Christer Ericson's "Real-Time Collision Detection"
 *       and is based on the Sutherland-Hodgeman clipping algorithm.
 *       We are only keeping the "back" polygon, w.r.t. that algorithm.
 * \see classifyPointAxisPlane for a description of how index maps to coordinates.
 */
template<typename T, int NDIMS>
void clipAxisPlane(const Polygon<T,NDIMS>* prevPoly, Polygon<T,NDIMS>* currentPoly, int index, T val)
{
  typedef Point<T,NDIMS> PointType;

  currentPoly->clear();
  int numVerts = prevPoly->numVertices();

  if(numVerts == 0)
    return;

  // Initialize point a with the last vertex of the polygon
  const PointType* a = &(*prevPoly)[numVerts-1];
  int aSide = classifyPointAxisPlane(*a, index, val);

  for(int i=0; i< numVerts; ++i)
  {
    const PointType* b = &(*prevPoly)[i];
    int bSide = classifyPointAxisPlane(*b, index, val);

    switch(bSide)
    {
    case PtPlaneClassifier::POINT_IN_FRONT_OF_PLANE:
      if(aSide == PtPlaneClassifier::POINT_BEHIND_PLANE)
      {
        currentPoly->addVertex(findIntersectionPoint(*a, *b, index, val));
      }
      break;
    case PtPlaneClassifier::POINT_ON_PLANE:
      if(aSide == PtPlaneClassifier::POINT_BEHIND_PLANE)
      {
        currentPoly->addVertex(*b);
      }
      break;
    case PtPlaneClassifier::POINT_BEHIND_PLANE:
      switch(aSide)
      {
      case PtPlaneClassifier::POINT_IN_FRONT_OF_PLANE:
        currentPoly->addVertex(findIntersectionPoint(*a, *b, index, val));
        currentPoly->addVertex(*b);
        break;
      case PtPlaneClassifier::POINT_ON_PLANE:
        currentPoly->addVertex(*a);
        currentPoly->addVertex(*b);
        break;
      case PtPlaneClassifier::POINT_BEHIND_PLANE:
        currentPoly->addVertex(*b);
        break;
      default:
        break;
      }
      break;
    default:
        break;
    }

    // swap a and b
    a = b;
    aSide = bSide;
  }

}

} // namespace detail
} // namespace primal
} // namespace axom

#endif // PRIMAL_CLIPPING_IMPL_HPP_

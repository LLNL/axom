// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file MIRUtilities.hpp
 * 
 * \brief Contains a set of functions that provide general utility
 *        within Axom's MIR component.
 * 
 */

#ifndef __MIR_UTILITIES_H__
#define __MIR_UTILITIES_H__

#include "ZooClippingTables.hpp"

//--------------------------------------------------------------------------------

namespace axom
{
namespace mir
{
namespace utilities
{

//--------------------------------------------------------------------------------

  /**
   * \brief Determines the number of vertices of the given shape in the finite element zoo.
   * 
   * \param shape  The shape type from the finite element zoo.
   * 
   * \return  THe number of vertices of the shape.
   */
  inline int  numVerts(mir::Shape shape)
  {
    int numVertices = -1;
    switch (shape)
    {
      case mir::Shape::Triangle:
        numVertices = 3;
        break;
      case mir::Shape::Quad:
        numVertices = 4;
        break;
      case mir::Shape::Tetrahedron:
        numVertices = 4;
        break;
      case mir::Shape::Pyramid:
        numVertices = 5;
        break;
      case mir::Shape::Triangular_Prism:
        numVertices = 6;
        break;
      case mir::Shape::Hexahedron:
        numVertices = 8;
        break;
      default:
        printf("Invalid shape. Cannot determine numVerts().\n");
    }
    return numVertices;
  }

//--------------------------------------------------------------------------------

  /**
   * \brief Determines the maximum number of possible vertices of the given shape in the finite element zoo.
   *        This number includes the midpoint vertices between each of the original shape's vertices.
   * 
   * \param shape  The shape type from the finite element zoo.
   * 
   * \return  THe number of vertices of the shape.
   */
  inline int  maxPossibleNumVerts(mir::Shape shape)
  {
    int numVertices = -1;
    switch (shape)
    {
      case mir::Shape::Triangle:
        numVertices = 6;
        break;
      case mir::Shape::Quad:
        numVertices = 8;
        break;
      case mir::Shape::Tetrahedron:
        numVertices = 10;
        break;
      case mir::Shape::Pyramid:
        numVertices = 13;
        break;
      case mir::Shape::Triangular_Prism:
        numVertices = 15;
        break;
      case mir::Shape::Hexahedron:
        numVertices = 20;
        break;
      default:
        printf("Invalid shape. Cannot determine maxPossibleNumVerts().\n");
    }
    return numVertices;
  }

//--------------------------------------------------------------------------------

  /**
   * \brief Performs linear interpolation between the two given float values.
   * 
   * \param f0  The first float value.
   * \param f1  The second float value.
   * \param t  The percent of the distance from the first float value to the second.
   * 
   * \return  The interpolated value.
   */
  inline axom::float64  lerpFloat(const axom::float64 f0, 
                                  const axom::float64 f1, 
                                  const axom::float64 t)
  {
    return (1 - t) * f0 + t * f1;
  }

//--------------------------------------------------------------------------------

  /**
   * \brief Performs linear interpolation between the two vertex positions.
   * 
   * \param vertexOnePos  The position of the first vertex.
   * \param vertexTwoPos  The position of the second vertex.
   * \param t  The percent of the distance from vertex one to vertex two to interpolate at.
   * 
   * \return  The interpolated position.
   */
  inline mir::Point2  interpolateVertexPosition(const mir::Point2& vertexOnePos, 
                                                const mir::Point2& vertexTwoPos, 
                                                const float t)
  {
    mir::Point2 interpolatedPoint;
    interpolatedPoint.m_x = lerpFloat(vertexOnePos.m_x, vertexTwoPos.m_x, t);
    interpolatedPoint.m_y = lerpFloat(vertexOnePos.m_y, vertexTwoPos.m_y, t);
    return interpolatedPoint;
  }

//--------------------------------------------------------------------------------

  /**
   * \brief  Returns the local vertex ID of the from/to vertex that is one of 
   *         the two endpoints that the edge the given midpoint is on.
   * 
   * \param shapeType  The shape type from the finite element zoo.
   * \param midpointVertexID  The ID of the vertex between the two endpoints.
   * \param isFromVertex  A flag denoting which of the two edge endpoints to return.
   * 
   * \return The vertex ID of one of the endpoints.
   */
  inline int getEdgeEndpoint(const mir::Shape shapeType,
                             const int midpointVertexID, 
                             const bool isFromVertex)
  {
    switch(shapeType)
    {
      case mir::Shape::Triangle:
        if ( midpointVertexID == 3 && isFromVertex  ) { return 0; }
        if ( midpointVertexID == 3 && !isFromVertex ) { return 1; }
        if ( midpointVertexID == 4 && isFromVertex  ) { return 1; }
        if ( midpointVertexID == 4 && !isFromVertex ) { return 2; }
        if ( midpointVertexID == 5 && isFromVertex  ) { return 2; }
        if ( midpointVertexID == 5 && !isFromVertex ) { return 0; }
        break;
      case mir::Shape::Quad:
        if ( midpointVertexID == 4 && isFromVertex  ) { return 0; }
        if ( midpointVertexID == 4 && !isFromVertex ) { return 1; }
        if ( midpointVertexID == 5 && isFromVertex  ) { return 1; }
        if ( midpointVertexID == 5 && !isFromVertex ) { return 2; }
        if ( midpointVertexID == 6 && isFromVertex  ) { return 2; }
        if ( midpointVertexID == 6 && !isFromVertex ) { return 3; }
        if ( midpointVertexID == 7 && isFromVertex  ) { return 3; }
        if ( midpointVertexID == 7 && !isFromVertex ) { return 0; }
        break;
      default:
        printf("Edge endpoint case not implemented.\n");
        return -1;
        break;
    }
    return -1;
  }

//--------------------------------------------------------------------------------

/**
 * \brief Calculate the distance between the two given points.
 * 
 * \param p0  The first point.
 * \param p1  The second point.
 * 
 * \return The distance between the two points.
 */
inline axom::float64 distance(mir::Point2 p0, mir::Point2 p1)
{
  return sqrt( ((p1.m_x - p0.m_x) * (p1.m_x - p0.m_x)) + ((p1.m_y - p0.m_y) * (p1.m_y - p0.m_y)) );
}

//--------------------------------------------------------------------------------

}
}
}

#endif
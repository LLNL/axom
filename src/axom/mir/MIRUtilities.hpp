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
      case mir::Shape::Tetrahedron:
        if ( midpointVertexID == 4 && isFromVertex  ) { return 0; }
        if ( midpointVertexID == 4 && !isFromVertex ) { return 1; }
        if ( midpointVertexID == 5 && isFromVertex  ) { return 0; }
        if ( midpointVertexID == 5 && !isFromVertex ) { return 2; } 
        if ( midpointVertexID == 6 && isFromVertex  ) { return 0; }
        if ( midpointVertexID == 6 && !isFromVertex ) { return 3; }
        if ( midpointVertexID == 7 && isFromVertex  ) { return 1; }
        if ( midpointVertexID == 7 && !isFromVertex ) { return 2; }
        if ( midpointVertexID == 8 && isFromVertex  ) { return 2; }
        if ( midpointVertexID == 8 && !isFromVertex ) { return 3; }
        if ( midpointVertexID == 9 && isFromVertex  ) { return 3; }
        if ( midpointVertexID == 9 && !isFromVertex ) { return 1; }
      case mir::Shape::Pyramid:
        if ( midpointVertexID ==  5 && isFromVertex  ) { return 0; }
        if ( midpointVertexID ==  5 && !isFromVertex ) { return 1; } 
        if ( midpointVertexID ==  6 && isFromVertex  ) { return 0; }
        if ( midpointVertexID ==  6 && !isFromVertex ) { return 2; }
        if ( midpointVertexID ==  7 && isFromVertex  ) { return 0; }
        if ( midpointVertexID ==  7 && !isFromVertex ) { return 3; }
        if ( midpointVertexID ==  8 && isFromVertex  ) { return 0; }
        if ( midpointVertexID ==  8 && !isFromVertex ) { return 4; }
        if ( midpointVertexID ==  9 && isFromVertex  ) { return 1; }
        if ( midpointVertexID ==  9 && !isFromVertex ) { return 2; }
        if ( midpointVertexID == 10 && isFromVertex  ) { return 2; }
        if ( midpointVertexID == 10 && !isFromVertex ) { return 3; }
        if ( midpointVertexID == 11 && isFromVertex  ) { return 3; }
        if ( midpointVertexID == 11 && !isFromVertex ) { return 4; }
        if ( midpointVertexID == 12 && isFromVertex  ) { return 4; }
        if ( midpointVertexID == 12 && !isFromVertex ) { return 1; }
      case mir::Shape::Triangular_Prism:
        if ( midpointVertexID ==  6 && isFromVertex  ) { return 0; }
        if ( midpointVertexID ==  6 && !isFromVertex ) { return 1; }
        if ( midpointVertexID ==  7 && isFromVertex  ) { return 1; }
        if ( midpointVertexID ==  7 && !isFromVertex ) { return 2; }
        if ( midpointVertexID ==  8 && isFromVertex  ) { return 2; }
        if ( midpointVertexID ==  8 && !isFromVertex ) { return 0; }
        if ( midpointVertexID ==  9 && isFromVertex  ) { return 0; }
        if ( midpointVertexID ==  9 && !isFromVertex ) { return 3; }
        if ( midpointVertexID == 10 && isFromVertex  ) { return 1; }
        if ( midpointVertexID == 10 && !isFromVertex ) { return 4; }
        if ( midpointVertexID == 11 && isFromVertex  ) { return 2; }
        if ( midpointVertexID == 11 && !isFromVertex ) { return 5; }
        if ( midpointVertexID == 12 && isFromVertex  ) { return 3; }
        if ( midpointVertexID == 12 && !isFromVertex ) { return 4; }
        if ( midpointVertexID == 13 && isFromVertex  ) { return 4; }
        if ( midpointVertexID == 13 && !isFromVertex ) { return 5; }
        if ( midpointVertexID == 14 && isFromVertex  ) { return 5; }
        if ( midpointVertexID == 14 && !isFromVertex ) { return 3; }
      case mir::Shape::Hexahedron:
        if ( midpointVertexID ==  8 && isFromVertex  ) { return 0; }
        if ( midpointVertexID ==  8 && !isFromVertex ) { return 1; }
        if ( midpointVertexID ==  9 && isFromVertex  ) { return 1; }
        if ( midpointVertexID ==  9 && !isFromVertex ) { return 2; }
        if ( midpointVertexID == 10 && isFromVertex  ) { return 2; }
        if ( midpointVertexID == 10 && !isFromVertex ) { return 3; }
        if ( midpointVertexID == 11 && isFromVertex  ) { return 3; }
        if ( midpointVertexID == 11 && !isFromVertex ) { return 0; }
        if ( midpointVertexID == 12 && isFromVertex  ) { return 0; }
        if ( midpointVertexID == 12 && !isFromVertex ) { return 4; }
        if ( midpointVertexID == 13 && isFromVertex  ) { return 1; }
        if ( midpointVertexID == 13 && !isFromVertex ) { return 5; }
        if ( midpointVertexID == 14 && isFromVertex  ) { return 2; }
        if ( midpointVertexID == 14 && !isFromVertex ) { return 6; }
        if ( midpointVertexID == 15 && isFromVertex  ) { return 3; }
        if ( midpointVertexID == 15 && !isFromVertex ) { return 7; }
        if ( midpointVertexID == 16 && isFromVertex  ) { return 4; }
        if ( midpointVertexID == 16 && !isFromVertex ) { return 5; }
        if ( midpointVertexID == 17 && isFromVertex  ) { return 5; }
        if ( midpointVertexID == 17 && !isFromVertex ) { return 6; }
        if ( midpointVertexID == 18 && isFromVertex  ) { return 6; }
        if ( midpointVertexID == 18 && !isFromVertex ) { return 7; }
        if ( midpointVertexID == 19 && isFromVertex  ) { return 7; }
        if ( midpointVertexID == 19 && !isFromVertex ) { return 4; }
      default:
        printf("Edge endpoint case not implemented.\n");
        return -1;
        break;
    }
    return -1;
  }


//--------------------------------------------------------------------------------

}
}
}

#endif

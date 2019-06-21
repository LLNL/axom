// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file InterfaceReconstructor.hpp
 * 
 * \brief Contains the specification for the InterfaceReconstructor class.
 * 
 */

#ifndef __INTERFACE_RECONSTRUCTOR_H__
#define __INTERFACE_RECONSTRUCTOR_H__

#include "axom/core.hpp"  // for axom macros
#include "axom/slam.hpp"

#include "MIRMesh.hpp"
#include "CellData.hpp"
#include "ZooClippingTables.hpp"

#include <map>

namespace numerics = axom::numerics;
namespace slam = axom::slam;

namespace axom
{
namespace mir
{

  /**
   * \class InterfaceReconstructor
   * 
   * \brief A class that contains the functionality for taking an input mesh
   *        with mixed cells and outputs a mesh with clean cells.
   * 
   * \detail This class requires that the user create a mesh of type MIRMesh 
   *         and pass that into one of the reconstruction methods. There are 
   *         currently two reconstructions methods implemented: one is the 
   *         a zoo-based algorithm described in Meredith 2004, and the other 
   *         is the iterative version described in Meredith and Childs 2010.
   */
  class InterfaceReconstructor
  {
    public:

      /**
       * \brief Default constructor.
       */
      InterfaceReconstructor();


      /**
       * \brief Default destructor.
       */
      ~InterfaceReconstructor();


      /**
       * \brief  Performs material interface reconstruction using the zoo-based algorithm.
       * 
       * \param inputMesh  The mesh composed of mixed cells.
       * 
       * \return  The mesh composed of clean cells.
       */
      mir::MIRMesh  computeReconstructedInterface(mir::MIRMesh& inputMesh);


      /**
       * \brief Performs material interface reconstruction using an iterative version of the zoo-based algorithm.
       * 
       * \param inputMesh  The mesh made up of mixed cells.
       * \param numIterations  The number of iterations for which to run the algorithm.
       * \param percent  The percent of the difference to use when modifying the original element volume fractions.
       * 
       * \return  The mesh made up of clean cells.
       */
      mir::MIRMesh  computeReconstructedInterfaceIterative(mir::MIRMesh& inputMesh, const int numIterations, const axom::float64 percent);

    // private:

      /**
       * \brief A wrapper function that calls the appropriate splitting method based on the shape of the given element.
       * 
       * \param eID The ID of the element to be split.
       * \param matOneID  The ID of the first material to use for splitting the element,
       * \param matTwoID  The ID of the second material to use for splitting the element.
       * \param tempMesh  A pointer to the intermediate mesh that is currently being processed.
       * \param out_cellData  Container to store the output of splitting the element.
       */
      void  computeClippingPoints(const int eID, const int matOneID, const int matTwoID, mir::MIRMesh& tempMesh, CellData& out_cellData);


      /**
       * \brief Splits the specified quad element into clean cells of the two given material types
       *        based on the volume fractions of the two materials.
       * 
       * \param eID  The ID of the quad element.
       * \param matOneID  The ID of the first material to use for splitting the quad,
       * \param matTwoID  The ID of the second material to use for splitting the quad.
       * \param tempMesh  A pointer to the intermediate mesh that is currently being processed.
       * \param out_cellData  Container to store the output of splitting the quad.
       */
      void  computeQuadClippingPoints(const int eID, const int matOneID, const int matTwoID, mir::MIRMesh& tempMesh, CellData& out_cellData);


      /**
       * \brief Splits the specified triangle element into clean cells of the two given material types
       *        based on the volume fractions of the two materials.
       * 
       * \param eID  The ID of the triangle element.
       * \param matOneID  The ID of the first material to use for splitting the triangle,
       * \param matTwoID  The ID of the second material to use for splitting the triangle.
       * \param tempMesh  A pointer to the intermediate mesh that is currently being processed.
       * \param out_cellData  Container to store the output of splitting the triangle.
       */
      void  computeTriangleClippingPoints(const int eID, const int matOneID, const int matTwoID, mir::MIRMesh& tempMesh, CellData& out_cellData);


      /**
       * \brief Performs linear interpolation between the two vertex positions.
       * 
       * \param vertexOnePos  The position of the first vertex.
       * \param vertexTwoPos  The position of the second vertex.
       * \param t  The percent of the distance from vertex one to vertex two to interpolate at.
       * 
       * \return  The interpolated position.
       */
      mir::Point2  interpolateVertexPosition(const mir::Point2& vertexOnePos, const mir::Point2& vertexTwoPos, const float t);


      /**
       * \brief Performs linear interpolation between the two given float values.
       * 
       * \param f0  The first float value.
       * \param f1  The second float value.
       * \param t  The percent of the distance from the first float value to the second.
       * 
       * \return  The interpolated value.
       */
      axom::float64  lerpFloat(const axom::float64 f0, const axom::float64 f1, const axom::float64 t);


      /**
       * \brief Computes the t value as a percent from vertex one to vertex two based on the materials given.
       * 
       * \param vertexOneID  The ID of the first vertex.
       * \param vertexTwoID  The ID of the second vertex.
       * \param matOneID  The ID of the first material to use for interpolating between.
       * \param matTwoID  The ID of the second material to use for interpolating between.
       * \param tempMesh  The intermediate mesh that is currently being processed.
       * 
       * \return  The t value, which is the percent distance from vertex one to vertex two where the two volume fractions are equal.
       */
      axom::float64  computeClippingPointOnEdge(const int vertexOneID, const int vertexTwoID, const int matOneID, const int matTwoID, mir::MIRMesh& tempMesh);

      
      /**
       * \brief Generates the topology of the new elements resulting from a split.
       * 
       * \param newElements  An ordered map of the generated elements' IDs to a list of the vertex IDs in the local frame of the output element.
       * \param newVertices  An order map of the vertex IDs in the local frame of the output element to the generated elements' IDs it is associated with.
       * \param out_cellData  Container to store the topology data of the generated elements.
       */
      void  generateTopologyData(const std::map<int, std::vector<int> >& newElements, const std::map<int, std::vector<int> >& newVertices, CellData& out_cellData);


      /**
       * \brief Calculates the bit map representing the clipping case for a quad.
       * 
       * \param tempMesh  The intermediate mesh that is currently being processed.
       * \param matOneID  The ID of the first material being used for splitting.
       * \param matTwoID  The ID of the second material being used for splitting.
       * \param upperLeftVertex  The ID of the upper left vertex in the global frame of the original mesh.
       * \param lowerLeftVertex  The ID of the lower left vertex in the global frame of the original mesh.
       * \param lowerRightVertex  The ID of the lower right vertex in the global frame of the original mesh.
       * \param upperRightVertex  The ID of the upper right vertex in the global frame of the original mesh.
       * 
       * \return  The bitmap representing the clipping case.
       */
      unsigned int  determineQuadClippingCase(mir::MIRMesh& tempMesh, const int matOneID, const int matTwoID, const int upperLeftVertex, const int lowerLeftVertex, const int lowerRightVertex, const int upperRightVertex);

      /**
       * \brief Generate the vertex position data for the new elements resulting from splitting a quad.
       * 
       * \param newVertices  An order map of the vertex IDs in the local frame of the output element to the generated elements' IDs it is associated with.
       * \param tempMesh  The intermediate mesh that is currently being processed.
       * \param verticesClippingTValue  The set of t values where the quad should be split on each edge.
       * \param upperLeftVertex  The ID of the upper left vertex in the global frame of the original mesh.
       * \param lowerLeftVertex  The ID of the lower left vertex in the global frame of the original mesh.
       * \param lowerRightVertex  The ID of the lower right vertex in the global frame of the original mesh.
       * \param upperRightVertex  The ID of the upper right vertex in the global frame of the original mesh.
       * \param out_cellData  Container to store the vertex position data of the generated elements.
       */
      void  generateVertexPositionsFromQuad(const std::map<int, std::vector<int> >& newVertices, mir::MIRMesh& tempMesh, axom::float64* verticesClippingTValue, const int upperLeftVertex, const int lowerLeftVertex, const int lowerRightVertex, const int upperRightVertex, CellData& out_cellData);

      /**
       * \brief Generate the vertex volume fraction data for the new vertices resulting from splitting a quad.
       * 
       * \param newVertices  An order map of the vertex IDs in the local frame of the output element to the generated elements' IDs it is associated with.
       * \param tempMesh  The intermediate mesh that is currently being processed.
       * \param verticesClippingTValue  The set of t values where the quad should be split on each edge.
       * \param upperLeftVertex  The ID of the upper left vertex in the global frame of the original mesh.
       * \param lowerLeftVertex  The ID of the lower left vertex in the global frame of the original mesh.
       * \param lowerRightVertex  The ID of the lower right vertex in the global frame of the original mesh.
       * \param upperRightVertex  The ID of the upper right vertex in the global frame of the original mesh.
       * \param out_cellData  Container to store the vertex volume fraction data of the generated elements.
       */
      void  generateVertexVolumeFractionsFromQuad(std::map<int, std::vector<int> >& newVertices, mir::MIRMesh&tempMesh, axom::float64* verticesClippingTValue, const int upperLeftVertex, const int lowerLeftVertex, const int lowerRightVertex, const int upperRightVertex, CellData& out_cellData);

      /**
       * \brief Calculates the bit map representing the clipping case for a triangle.
       * 
       * \param tempMesh  The intermediate mesh that is currently being processed.
       * \param matOneID  The ID of the first material being used for splitting.
       * \param matTwoID  The ID of the second material being used for splitting.
       * \param upperVertex  The ID of the upper vertex in the global frame of the original mesh.
       * \param lowerLeftVertex  The ID of the lower left vertex in the global frame of the original mesh.
       * \param lowerRightVertex  The ID of the lower right vertex in the global frame of the original mesh.
       * 
       * \return  The bitmap representing the clipping case.
       */
      unsigned int  determineTriangleClippingCase(mir::MIRMesh& tempMesh, const int matOneID, const int matTwoID, const int upperVertex, const int lowerLeftVertex, const int lowerRightVertex);

      /**
       * \brief Generate the vertex position data for the new elements resulting from splitting a triangle.
       * 
       * \param newVertices  An order map of the vertex IDs in the local frame of the output element to the generated elements' IDs it is associated with.
       * \param tempMesh  The intermediate mesh that is currently being processed.
       * \param verticesClippingTValue  The set of t values where the quad should be split on each edge.
       * \param upperVertex  The ID of the upper vertex in the global frame of the original mesh.
       * \param lowerLeftVertex  The ID of the lower left vertex in the global frame of the original mesh.
       * \param lowerRightVertex  The ID of the lower right vertex in the global frame of the original mesh.
       * \param out_cellData  Container to store the vertex position data of the generated elements.
       */
      void  generateVertexPositionsFromTriangle(const std::map<int, std::vector<int> >& newVertices, mir::MIRMesh& tempMesh, axom::float64* verticesClippingTValue, const int upperVertex, const int lowerLeftVertex, const int lowerRightVertex, CellData& out_cellData);

      /**
       * \brief Generate the vertex volume fraction data for the new vertices resulting from splitting a triangle.
       * 
       * \param newVertices  An order map of the vertex IDs in the local frame of the output element to the generated elements' IDs it is associated with.
       * \param tempMesh  The intermediate mesh that is currently being processed.
       * \param verticesClippingTValue  The set of t values where the quad should be split on each edge.
       * \param upperVertex  The ID of the upper vertex in the global frame of the original mesh.
       * \param lowerLeftVertex  The ID of the lower left vertex in the global frame of the original mesh.
       * \param lowerRightVertex  The ID of the lower right vertex in the global frame of the original mesh.
       * \param out_cellData  Container to store the vertex volume fraction data of the generated elements.
       */
      void  generateVertexVolumeFractionsFromTriangle(const std::map<int, std::vector<int> >& newVertices, mir::MIRMesh& tempMesh, axom::float64* verticesClippingTValue, const int upperVertex, const int lowerLeftVertex, const int lowerRightVertex, CellData& out_cellData);

      /**
       * \brief Determines the mroe dominant material of the two given for the given element.
       * 
       * \param elementShape  An enumerator denoting the element's shape.
       * \param vertexIDs  A list of vertex IDs into the vertexVF param.
       * \param matOne  The ID of the first material.
       * \param matTwo  The ID of the second material.
       * \param vertexVF  The list of volume fractions associated with the given vertices in the vertexIDs param.
       * 
       * \return  The ID of the dominant material of the element.
       * 
       * \note The dominant element for the 2D cases will be the same as the material present at one of the 
       *       original vertices that existed prior to the split. So, if you can find this vertex and its dominant
       *       material, then you know the dominant material of this new element.
       */
      int  determineDominantMaterial(const Shape elementShape, const std::vector<int>& vertexIDs, const int matOne, const int matTwo, const std::vector<std::vector<axom::float64> >& vertexVF);

      private:
        mir::MIRMesh m_originalMesh;
  };
}
}
#endif
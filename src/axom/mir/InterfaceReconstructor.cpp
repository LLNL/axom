// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "InterfaceReconstructor.hpp"

namespace axom
{
namespace mir
{

//--------------------------------------------------------------------------------

  /// Default constructor
  InterfaceReconstructor::InterfaceReconstructor()
  {
    
  }

//--------------------------------------------------------------------------------

  /// Constructor
  InterfaceReconstructor::InterfaceReconstructor(mir::MIRMesh* _mesh)
  {
    mesh = _mesh;
  }

//--------------------------------------------------------------------------------

  /// Destructor
  InterfaceReconstructor::~InterfaceReconstructor()
  {

  }

//--------------------------------------------------------------------------------

  /// Computes the points where the element should be clipped based on the volume fractions of mat one and two.
  /// TODO: This method needs to place the computed cell vertices into the intermediate mesh (pass this in as an argument)
  void InterfaceReconstructor::computeClippingPoints(const int eID, const int matOneID, const int matTwoID, mir::MIRMesh* tempMesh,
                                                std::vector<mir::PosType>& _evInds, std::vector<mir::PosType>& _evBegins, std::vector<mir::PosType>& _veInds, std::vector<mir::PosType>& _veBegins,
                                                std::vector<mir::Point2>& _vertexPositions, std::vector<axom::float64*>& _materialInCell, int* _numVerts, int* _numElements, std::vector<std::vector<axom::float64> >& _newVolumeFractionsAtVerts)
  {
    // Find the vertices associated with each element
    auto elementVertices = mesh->bdry[eID];

    // Check if one of the two currently considered materials is not present at the vertices of the current element
    // if (mesh->materialVolumeFractionsElement[matOneID][eID] != 0.0 &&  mesh->materialVolumeFractionsElement[matTwoID][eID] != 0.0) // TODO: Figure out how to best handle this case
    // {
      // Triangle Case
      if (elementVertices.size() == 3)
      {
        computeTriangleClippingPoints(eID, matOneID, matTwoID, tempMesh, _evInds, _evBegins, _veInds, _veBegins, _vertexPositions, _materialInCell, _numVerts, _numElements, _newVolumeFractionsAtVerts);
      }
      // Quad Case
      if (elementVertices.size() == 4)
      {
        computeQuadClippingPoints(eID, matOneID, matTwoID, tempMesh, _evInds, _evBegins, _veInds, _veBegins, _vertexPositions, _materialInCell, _numVerts, _numElements, _newVolumeFractionsAtVerts);    // Compute how the cell should be decomposed into new cells
      }
    // }
    // else
    // {
      // printf("No cuts in element %d.\n", eID);

      // TODO: Add the current cells vertices and data back into the mesh

    // }
  }

//--------------------------------------------------------------------------------

  /// Computes the points where the triangle element should be clipped based on the volume fractions of mat one and two.
  void InterfaceReconstructor::computeTriangleClippingPoints(const int eID, const int matOneID, const int matTwoID, mir::MIRMesh* tempMesh,
                                                std::vector<mir::PosType>& _evInds, std::vector<mir::PosType>& _evBegins, std::vector<mir::PosType>& _veInds, std::vector<mir::PosType>& _veBegins,
                                                std::vector<mir::Point2>& _vertexPositions, std::vector<axom::float64*>&  _materialInCell, int* _numVerts, int* _numElements, std::vector<std::vector<axom::float64> >& _newVolumeFractionsAtVerts)
  {
    printf("Triangle clipping case not yet implemented.\n");
  }

//--------------------------------------------------------------------------------

  /// Computes how the given cell should be decomposed into new cells baed on teh vertex volume fractions of matOne and matTwo.
  void InterfaceReconstructor::computeQuadClippingPoints(const int eID, const int matOneID, const int matTwoID, mir::MIRMesh* tempMesh,
                                                std::vector<mir::PosType>& out_evInds, std::vector<mir::PosType>& out_evBegins, std::vector<mir::PosType>& out_veInds, std::vector<mir::PosType>& out_veBegins,
                                                std::vector<mir::Point2>& out_vertexPositions, std::vector<axom::float64*>&  out_materialInCell, int* out_numVerts, int* out_numElements, std::vector<std::vector<axom::float64> >& out_newVolumeFractionsAtVerts)
  {
    printf("Processing element %d: ", eID);

    /****************************************************************
     *                  DETERMINE THE CLIPPING CASE
     ****************************************************************/
    // Get the vertices of the current element to clip
    auto elementVertices = tempMesh->bdry[eID];

    // Determine which vertices correspond to upper left, lower left, lower right, and upper right vertices of the quad.
    int upperLeftVertex = elementVertices[0];
    int lowerLeftVertex = elementVertices[1];
    int lowerRightVertex = elementVertices[2];
    int upperRightVertex = elementVertices[3];

    // Determine the dominant color at each vertex
    int upperLeftColor = tempMesh->materialVolumeFractionsVertex[matOneID][upperLeftVertex] > tempMesh->materialVolumeFractionsVertex[matTwoID][upperLeftVertex] ? matOneID : matTwoID;
    int lowerLeftColor = tempMesh->materialVolumeFractionsVertex[matOneID][lowerLeftVertex] > tempMesh->materialVolumeFractionsVertex[matTwoID][lowerLeftVertex] ? matOneID : matTwoID;
    int lowerRightColor = tempMesh->materialVolumeFractionsVertex[matOneID][lowerRightVertex] > tempMesh->materialVolumeFractionsVertex[matTwoID][lowerRightVertex] ? matOneID : matTwoID;
    int upperRightColor = tempMesh->materialVolumeFractionsVertex[matOneID][upperRightVertex] > tempMesh->materialVolumeFractionsVertex[matTwoID][upperRightVertex] ? matOneID : matTwoID;

    // Create the index into the quad clipping lookup table using the dominant colors at each vertex
    unsigned int caseIndex = 0;
    if (upperLeftColor == matOneID)  caseIndex |= 8;
    if (lowerLeftColor == matOneID)  caseIndex |= 4;
    if (lowerRightColor == matOneID) caseIndex |= 2;
    if (upperRightColor == matOneID) caseIndex |= 1;

    printf("caseIndex: %d\n", caseIndex);

    /****************************************************************
     *                    GENERATE NEW ELEMENTS
     ****************************************************************/
    // Cell information indices
    std::map<int, std::vector<int> > newElements;   // hashmap of the new elements' vertices | Note: maps are ordered sets
    std::map<int, std::vector<int> > newVertices;   // hashmap of the new vertices' elements | Note: maps are ordered sets
    std::map<int, mir::Point2> newPoints;           // hashmap of the new vertices' positions | Note: maps are ordered sets
    int verticesPresent[8];                                   // Array of flags denoting whether the vertex is present in the current case or not
    axom::float64 verticesClippingTValue[8];                  // Array of t values that denote the percent value of where the edge should be clipped
    std::map<int, std::vector<axom::float64> > vertexVolumeFractionsMap;  // index: map[vID][matID] = vertexVolumeFractionValues[numMaterials] | Note: maps are ordered sets

    int currentElementIndex = 0;    // the next available element index

    // Create the new polygons based on the clipping case
    int i = 0;
    int numVertices = quadClipTable[caseIndex][i];
    while (numVertices != -1)   // for each new element in the current clipping case
    {
      // for each vertex of the new element
      for (int j = 0; j < numVertices; ++j)
      {
        // Find the id of the next vertex of the new element
        int vID = quadClipTable[caseIndex][i + (j+1)];

        // Associate the vertex and element together
        newElements[currentElementIndex].push_back(vID);
        newVertices[vID].push_back(currentElementIndex);
        verticesPresent[vID] = 1;

        // Find t using "bilinear" interpolation method for any vertex that is not one of the original 4 vertices
        if(vID == 4)
        {
          verticesClippingTValue[vID] = computeClippingPointOnEdge(upperLeftVertex, lowerLeftVertex, matOneID, matTwoID, tempMesh);  // these might be the wrong vertex indices
        }
        else if(vID == 5)
        {
          verticesClippingTValue[vID] = computeClippingPointOnEdge(lowerLeftVertex, lowerRightVertex, matOneID, matTwoID, tempMesh);
        }
        else if(vID == 6)
        {
          verticesClippingTValue[vID] = computeClippingPointOnEdge(lowerRightVertex, upperRightVertex, matOneID, matTwoID, tempMesh);
        }
        else if(vID == 7)
        {
          verticesClippingTValue[vID] = computeClippingPointOnEdge(upperRightVertex, upperLeftVertex, matOneID, matTwoID, tempMesh);
        }
      }

      // Increment the element index counter, marking the current element as being finished processed
      currentElementIndex++;

      // Increase index into lookup table to the next element
      i += (numVertices + 1);
      numVertices = quadClipTable[caseIndex][i];
    }

    /****************************************************************
     *                CALCULATE THE NEW CELLS' DATA
     ****************************************************************/
    // Calculate the total number of elements and vertices that were generated from splitting the current element
    out_numElements[0] = (int) newElements.size();
    out_numVerts[0] = (int) newVertices.size();

    // Generate the topology of the new elements (evInds, evBegins, etc)
    // Store the evInds and evBegins data in the output vectors
    int currentEVBeginIndex = 0;
    for (auto itr = newElements.begin(); itr != newElements.end(); itr++)
    {
      // Push the start index of the next element
      out_evBegins.push_back(currentEVBeginIndex);

      // Push the next element's vertices
      for (unsigned int vIndex = 0; vIndex < itr->second.size(); ++vIndex)
      {
        out_evInds.push_back(itr->second[vIndex]);
        ++currentEVBeginIndex;
      }
    }

    // Push the index that occurs after the last vertex
    out_evBegins.push_back(currentEVBeginIndex);
    
    // Store the veInds and veBegins data in the output vectors
    int currentVEBeginIndex = 0;
    for (auto itr = newVertices.begin(); itr != newVertices.end(); itr++)
    {
      // Push the start index of the vertex's elements
      out_veBegins.push_back(currentVEBeginIndex);

      // Push the next vertex's elements
      for (unsigned int eIndex = 0; eIndex < itr->second.size(); eIndex++)
      {
        out_veInds.push_back(itr->second[eIndex]);
        ++currentVEBeginIndex;
      }
    }
    
    // Push the index that occurs after the last element
    out_veBegins.push_back(currentVEBeginIndex);

    // Find the position coordinates of the vertices and then store the positions of the vertices in the return vector
    for (auto itr = newVertices.begin(); itr != newVertices.end(); itr++)
    {
      int vID = itr->first;
      if (verticesPresent[vID] == 1)
      {
        if (vID == 0)
          newPoints[vID] = tempMesh->vertexPositions[upperLeftVertex];
        if (vID == 1)
          newPoints[vID] = tempMesh->vertexPositions[lowerLeftVertex];
        if (vID == 2)
          newPoints[vID] = tempMesh->vertexPositions[lowerRightVertex];
        if (vID == 3)
          newPoints[vID] = tempMesh->vertexPositions[upperRightVertex];
        if (vID == 4)
          newPoints[vID] = interpolateVertexPosition(tempMesh->vertexPositions[upperLeftVertex], tempMesh->vertexPositions[lowerLeftVertex], verticesClippingTValue[vID]);
        if (vID == 5)
          newPoints[vID] = interpolateVertexPosition(tempMesh->vertexPositions[lowerLeftVertex], tempMesh->vertexPositions[lowerRightVertex], verticesClippingTValue[vID]);
        if (vID == 6)
          newPoints[vID] = interpolateVertexPosition(tempMesh->vertexPositions[lowerRightVertex], tempMesh->vertexPositions[upperRightVertex], verticesClippingTValue[vID]);
        if (vID == 7)
          newPoints[vID] = interpolateVertexPosition(tempMesh->vertexPositions[upperRightVertex], tempMesh->vertexPositions[upperLeftVertex], verticesClippingTValue[vID]);
      }
    } 

    // Store the positions of the vertices in the return vector
    for (auto itr = newPoints.begin(); itr != newPoints.end(); itr++)
    {
      out_vertexPositions.push_back(itr->second);
    }

    // Calculate the vertex fractions at each vertex (use t value!)
    // Make sure the output volume fractions containers are the proper size
    out_newVolumeFractionsAtVerts.resize(tempMesh->numMaterials);
    
    for (auto itr = newVertices.begin(); itr != newVertices.end(); itr++)
    {
      int vID = itr->first;
      // Ensure that the current vertex being considered is actually present in the current clipping case
      if (verticesPresent[vID] == 1)
      {
        // vertexVolumeFractionsMap[vID].resize(tempMesh->numMaterials);

        for (int matID = 0; matID < tempMesh->numMaterials; ++matID)
        {
          if (vID == 0)
            out_newVolumeFractionsAtVerts[matID].push_back(tempMesh->materialVolumeFractionsVertex[matID][upperLeftVertex]);
          if (vID == 1)
            out_newVolumeFractionsAtVerts[matID].push_back(tempMesh->materialVolumeFractionsVertex[matID][lowerLeftVertex]);
          if (vID == 2)
            out_newVolumeFractionsAtVerts[matID].push_back(tempMesh->materialVolumeFractionsVertex[matID][lowerRightVertex]);
          if (vID == 3)
            out_newVolumeFractionsAtVerts[matID].push_back(tempMesh->materialVolumeFractionsVertex[matID][upperRightVertex]);
          if (vID == 4)
            out_newVolumeFractionsAtVerts[matID].push_back(lerpFloat(tempMesh->materialVolumeFractionsVertex[matID][upperLeftVertex], tempMesh->materialVolumeFractionsVertex[matID][lowerLeftVertex], verticesClippingTValue[vID]));
          if (vID == 5)
            out_newVolumeFractionsAtVerts[matID].push_back(lerpFloat(tempMesh->materialVolumeFractionsVertex[matID][lowerLeftVertex], tempMesh->materialVolumeFractionsVertex[matID][lowerRightVertex], verticesClippingTValue[vID]));
          if (vID == 6)
            out_newVolumeFractionsAtVerts[matID].push_back(lerpFloat(tempMesh->materialVolumeFractionsVertex[matID][lowerRightVertex], tempMesh->materialVolumeFractionsVertex[matID][upperRightVertex], verticesClippingTValue[vID]));
          if (vID == 7)
            out_newVolumeFractionsAtVerts[matID].push_back(lerpFloat(tempMesh->materialVolumeFractionsVertex[matID][upperRightVertex], tempMesh->materialVolumeFractionsVertex[matID][upperLeftVertex], verticesClippingTValue[vID]));
        }
      }
    }

    // Modify the evIndex values to account for the fact that perhaps not all 8 possible vertices are present
    int evIndexSubtract[8];
    for (int i = 0; i < 8; ++i)
    {
      if (verticesPresent[i] == 1)
        evIndexSubtract[i] = 0;
      else
        evIndexSubtract[i] = 1;
    }
    for (int i = 1; i < 8; ++i)
      evIndexSubtract[i] += evIndexSubtract[i - 1];

    for (unsigned int i = 0; i < out_evInds.size(); ++i)
    {
      out_evInds[i] -= evIndexSubtract[ out_evInds[i] ];
    }



    /****************************************************************
     *                       DEBUG PRINTING
     ****************************************************************/
    {
      // TESTING: Print out the values by which the evIndex values need to be substracted 
      //          in order to ensure there are no gaps in the vertex ids in the resulting mesh.
      // printf("  evIndexSubtract: { ");
      // for (int i = 0; i < 8; i++)
      // {
      //   printf("%d ", evIndexSubtract[i]);
      // }
      // printf("}\n");

      // // TESTING: Print out the basic cell information
      // printf("  out_numElements: %d\n", out_numElements[0]);
      // printf("  out_numVerts: %d\n", out_numVerts[0]);

      // // TESTING: Print out the mesh topology information
      // printf("  out_evInds: {");
      // for (int i = 0; i < out_evInds.size(); i++)
      // {
      //   printf("%d ", out_evInds[i]);
      // }
      // printf("}\n");

      // printf("  out_evBegins: {");
      // for (int i = 0; i < out_evBegins.size(); i++)
      // {
      //   printf("%d ", out_evBegins[i]);
      // }
      // printf("}\n");

      // printf("  out_veInds: {");
      // for (int i = 0; i < out_veInds.size(); i++)
      // {
      //   printf("%d ", out_veInds[i]);
      // }
      // printf("}\n");

      // printf("  out_veBegins: {");
      // for (int i = 0; i < out_veBegins.size(); i++)
      // {
      //   printf("%d ", out_veBegins[i]);
      // }
      // printf("}\n");

      // // TESTING: Print out the positions of the vertices
      // for (auto itr = newPoints.begin(); itr != newPoints.end(); itr++)
      // {
      //   printf("  Vertex %d: (%f, %f)\n", itr->first, itr->second.m_x, itr->second.m_y);
      // }

      // TESTING: Print out the volume fractions at the vertices from the map
      // for (auto itr = vertexVolumeFractionsMap.begin(); itr != vertexVolumeFractionsMap.end(); itr++)
      // {
      //   int vID = itr->first;
      //   printf("  Vertex %d Volume Fractions:\n", vID);
      //   for (int matID = 0; matID < tempMesh->numMaterials; ++matID)
      //   {
      //     printf("    Mat %d: %f\n", matID, vertexVolumeFractionsMap[vID][matID]);
      //   }
      // }
    
      // TESTING: Print out the volume fractions at the vertices from the vector of materials containing a vector of vertex volume fraction values
      // for (int matID = 0; matID < out_newVolumeFractionsAtVerts.size(); ++matID)
      // {
      //   printf("  Material %d:\n", matID);
      //   for (int vID = 0; vID < out_newVolumeFractionsAtVerts[matID].size(); ++vID)
      //   {
      //     printf("    Vertex %d: %f\n", vID, out_newVolumeFractionsAtVerts[matID][vID]);
      //   }
      // }
    }
  }

//--------------------------------------------------------------------------------

/// Performs linear interpolation between the two given vertex positions.
mir::Point2 InterfaceReconstructor::interpolateVertexPosition(mir::Point2 vertexOnePos, mir::Point2 vertexTwoPos, float t)
{
  mir::Point2 interpolatedPoint;
  interpolatedPoint.m_x = (1 - t) * vertexOnePos.m_x + t * vertexTwoPos.m_x;
  interpolatedPoint.m_y = (1 - t) * vertexOnePos.m_y + t * vertexTwoPos.m_y;
  return interpolatedPoint;
}

//--------------------------------------------------------------------------------

/// Performs linear interpolation between the two given float values
axom::float64 InterfaceReconstructor::lerpFloat(axom::float64 f0, axom::float64 f1, axom::float64 t)
{
  return (1 - t) * f0 + t * f1;
}

//--------------------------------------------------------------------------------

/// Computes the t value as a percent from vertexOne to vertexTwo based on the two materials given.
/// The t value is the place where this edge should be clipped based on the two materials currently being considered.
axom::float64 InterfaceReconstructor::computeClippingPointOnEdge(const int vertexOneID, const int vertexTwoID, const int matOneID, const int matTwoID, mir::MIRMesh* tempMesh)
{
  axom::float64 ret = 0.0;

  axom::float64 numerator = tempMesh->materialVolumeFractionsVertex[matTwoID][vertexOneID] - mesh->materialVolumeFractionsVertex[matOneID][vertexOneID];
  axom::float64 denominator = -mesh->materialVolumeFractionsVertex[matOneID][vertexOneID]
                                    + mesh->materialVolumeFractionsVertex[matOneID][vertexTwoID]
                                    + mesh->materialVolumeFractionsVertex[matTwoID][vertexOneID]
                                    - mesh->materialVolumeFractionsVertex[matTwoID][vertexTwoID];
        
  if (denominator != 0.0)
    ret = numerator / denominator;

  return ret;
}

//--------------------------------------------------------------------------------

mir::MIRMesh InterfaceReconstructor::computeReconstructedInterface()
{
  // Initialize the final mesh to be the same as the input mesh
  mir::MIRMesh finalMesh(mesh);                                  // TODO: This might just copy the reference to the original mesh, and not copy it. Could cause bugs when it gets overwritten.

  // For each material in the mesh, split the mesh on the current material and generate a new mesh (input into next iteration)
  for (int matID = 0; matID < mesh->numMaterials - 1; ++matID)
  {
    // Copy the mesh to be split
    mir::MIRMesh intermediateMesh(&finalMesh);

    // Update the materials upon which the split will occur
    int matOne = matID;
    int matTwo = matID + 1;

    // Process/split each element
    std::vector<mir::PosType> temp_evInds[intermediateMesh.elems.size()];           // Store the vertex indices adjacent to each element for each cell
    std::vector<mir::PosType> temp_evBegins[intermediateMesh.elems.size()];         // Store the starting indices into temp_evBewgins for each element-vertex for each cell
    std::vector<mir::PosType> temp_veInds[intermediateMesh.elems.size()];           // Store the element indices adjacent to each vertex for each cell
    std::vector<mir::PosType> temp_veBegins[intermediateMesh.elems.size()];         // Store the starting indices into temp_veInds for each vertex-element for each cell

    std::vector<mir::Point2> temp_vertexPositions[intermediateMesh.elems.size()];   // Store the positions of the vertices generated for each cell
    std::vector<axom::float64*> materialInCell[intermediateMesh.elems.size()];  // Note: after splitting, the cells should all be clean and thus will have a vf of 1.0 or 0.0 for all materials

    int temp_numVerts[intermediateMesh.elems.size()];                              // Store the number of vertices generated for each cell
    int temp_numElements[intermediateMesh.elems.size()];                           // Store the number of elements generated for each cell

    std::vector<std::vector<axom::float64> > temp_volumeFractionsVertex[intermediateMesh.elems.size()];

    for (int eID = 0; eID < intermediateMesh.elems.size(); ++eID)
    {
      computeClippingPoints(eID, matOne, matTwo, &intermediateMesh, 
                            temp_evInds[eID], temp_evBegins[eID], temp_veInds[eID], temp_veBegins[eID], temp_vertexPositions[eID], 
                            materialInCell[eID], &(temp_numVerts[eID]), &(temp_numElements[eID]), temp_volumeFractionsVertex[eID]);
    }

    // Copy the generated cell information into CellData structures
    std::vector<CellData> cellSplitInfo;
    for (int eID = 0; eID < intermediateMesh.elems.size(); ++eID)
    {
      CellData nextCellSplitInfo(temp_numVerts[eID], temp_numElements[eID], temp_evInds[eID], temp_evBegins[eID], temp_veInds[eID], temp_veBegins[eID], temp_vertexPositions[eID], temp_volumeFractionsVertex[eID]);
      cellSplitInfo.push_back(nextCellSplitInfo);
    }

    // Merge each of the cells into the first CellData struct
    for (unsigned int eID = 1; eID < cellSplitInfo.size(); ++eID)
      cellSplitInfo[0].mergeCell(cellSplitInfo[eID]);

    mir::VertSet combined_verts(cellSplitInfo[0].numVerts);
    mir::ElemSet combined_elems(cellSplitInfo[0].numElems);
    
    // TODO (later): Calculate the element volume fractions (note, they should all be either 0 or 1, since after the splits every cell should be clean)

    // Create the final, processed mesh
    mir::MIRMesh processedMesh;
    processedMesh.InitializeMesh(cellSplitInfo[0].evInds, cellSplitInfo[0].evBegins, cellSplitInfo[0].veInds, cellSplitInfo[0].veBegins, combined_verts, combined_elems, intermediateMesh.numMaterials);
    processedMesh.constructMeshRelations();
    processedMesh.constructMeshVolumeFractionsVertex(cellSplitInfo[0].vertexVolumeFractions);
    processedMesh.constructVertexPositionMap(cellSplitInfo[0].vertexPositions.data());  

    // Store the current mesh to be passed into the next iteration
    finalMesh = processedMesh;
  }

  // TODO: Run post-processing on mesh to meld the mesh vertices that are within a certain epsilon

  return finalMesh;
}

//--------------------------------------------------------------------------------

}
}

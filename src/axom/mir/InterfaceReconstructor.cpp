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

InterfaceReconstructor::InterfaceReconstructor()
{
  
}

//--------------------------------------------------------------------------------

InterfaceReconstructor::~InterfaceReconstructor()
{

}

//--------------------------------------------------------------------------------

mir::MIRMesh InterfaceReconstructor::computeReconstructedInterface(mir::MIRMesh& inputMesh)
{
  // Store a reference to the original mesh
  m_originalMesh = inputMesh;

  // Initialize the final mesh to be the same as the input mesh
  mir::MIRMesh finalMesh(m_originalMesh);

  // For each material in the mesh, split the mesh on the current material and generate a new mesh (input into next iteration)
  for (int matID = 0; matID < m_originalMesh.m_numMaterials; ++matID)
  {
    // Copy the mesh to be split
    mir::MIRMesh intermediateMesh(&finalMesh);

    // Create an array to store the output of each element being split.
    CellData temp_cellData[intermediateMesh.m_elems.size()];

    // Process/split each element
    for (int eID = 0; eID < intermediateMesh.m_elems.size(); ++eID)
    {
      // Update the materials upon which the split will occur
      int currentDominantMat = intermediateMesh.m_elementDominantMaterials[eID];
      int matOne = matID;

      computeClippingPoints(eID, currentDominantMat, matOne, intermediateMesh, temp_cellData[eID]);
    }

    // Merge each of the cells into the first CellData struct
    for (int eID = 1; eID < intermediateMesh.m_elems.size(); ++eID)
    {
      temp_cellData[0].mergeCell(temp_cellData[eID]);
    }

    mir::VertSet combined_verts(temp_cellData[0].m_numVerts);
    mir::ElemSet combined_elems(temp_cellData[0].m_numElems);

    // Create the final, processed mesh
    mir::MIRMesh processedMesh;
    processedMesh.initializeMesh(combined_verts, combined_elems, intermediateMesh.m_numMaterials, temp_cellData[0].m_topology, temp_cellData[0].m_mapData);
    processedMesh.constructMeshVolumeFractionsVertex(temp_cellData[0].m_mapData.m_vertexVolumeFractions);

    // Store the current mesh to be passed into the next iteration
    finalMesh = processedMesh;
    finalMesh.constructMeshRelations();
  }

  // TODO: Run post-processing on mesh to meld the mesh vertices that are within a certain epsilon

  return finalMesh;
}

//--------------------------------------------------------------------------------

mir::MIRMesh InterfaceReconstructor::computeReconstructedInterfaceIterative(mir::MIRMesh& inputMesh, const int numIterations, const axom::float64 percent)
{
  int numElems = inputMesh.m_elems.size();
  int numMaterials = inputMesh.m_numMaterials;

  // Make a copy of the original input mesh
  mir::MIRMesh meshToImprove(inputMesh);

  // Calculate the reconstruction on the unmodified, input mesh
  mir::MIRMesh resultingMesh = computeReconstructedInterface(meshToImprove);

  for (int it = 0; it < numIterations; ++it)
  {
    // Calculate the output element volume fractions of the resulting output mesh
    std::vector<std::vector<axom::float64> > resultingElementVF = resultingMesh.computeOriginalElementVolumeFractions();

    // Initialize the vector to store the improved element volume fractions
    std::vector<std::vector<axom::float64> > improvedElementVF;
    improvedElementVF.resize(numMaterials);
      
    // Modify the copy of the original mesh's element volume fractions by percent difference (modify meshToImprove's elementVF)
    for (int matID = 0; matID < numMaterials; ++matID)
    {
      for (int eID = 0; eID < numElems; ++eID)
      {
        axom::float64 difference = inputMesh.m_materialVolumeFractionsElement[matID][eID] - resultingElementVF[matID][eID];
        improvedElementVF[matID].push_back(   inputMesh.m_materialVolumeFractionsElement[matID][eID] + ( difference * percent ) );
      }
    }

    // Reconstruct the element AND vertex VFs
    meshToImprove.constructMeshVolumeFractionsMaps(improvedElementVF);

    // Calculate the reconstruction on the modified, input mesh
    resultingMesh = computeReconstructedInterface(meshToImprove);
  }

  return resultingMesh;
}

//--------------------------------------------------------------------------------

void InterfaceReconstructor::computeClippingPoints(const int eID, const int matOneID, const int matTwoID, mir::MIRMesh& tempMesh, CellData& out_cellData)
{
  // Find the vertices associated with each element
  auto elementVertices = tempMesh.m_bdry[eID];

  // Triangle Case
  if (elementVertices.size() == 3)
  {
    computeTriangleClippingPoints(eID, matOneID, matTwoID, tempMesh, out_cellData);
  }
  // Quad Case
  if (elementVertices.size() == 4)
  {
    computeQuadClippingPoints(eID, matOneID, matTwoID, tempMesh, out_cellData);
  }
}

//--------------------------------------------------------------------------------

void InterfaceReconstructor::computeQuadClippingPoints(const int eID, const int matOneID, const int matTwoID, mir::MIRMesh& tempMesh, CellData& out_cellData)
{
  // Determine the clipping case
  auto elementVertices = tempMesh.m_bdry[eID];

  int upperLeftVertex = elementVertices[0];
  int lowerLeftVertex = elementVertices[1];
  int lowerRightVertex = elementVertices[2];
  int upperRightVertex = elementVertices[3];

  unsigned int caseIndex = determineQuadClippingCase(tempMesh, matOneID, matTwoID, upperLeftVertex, lowerLeftVertex, lowerRightVertex, upperRightVertex);

  // Generate new elements
  std::map<int, std::vector<int> > newElements;             // hashmap of the new elements' vertices | Note: maps are ordered sets
  std::map<int, std::vector<int> > newVertices;             // hashmap of the new vertices' elements | Note: maps are ordered sets
  int verticesPresent[8] = {0,0,0,0,0,0,0,0};               // Array of flags denoting whether the vertex is present in the current case or not
  axom::float64 verticesClippingTValue[8] = {0,0,0,0,0,0,0,0}; // Array of t values that denote the percent value of where the edge should be clipped

  // Create the new polygons based on the clipping case
  int currentElementIndex = 0;    // the next available element index
  int i = 0;
  int numVertices = quadClipTable[caseIndex][i];
  
  // for each new element in the current clipping case
  while (numVertices != -1)   
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

      // Find t using linear interpolation for any vertex that is not one of the original 4 vertices
      if(vID == 4)
      {
        verticesClippingTValue[vID] = computeClippingPointOnEdge(upperLeftVertex, lowerLeftVertex, matOneID, matTwoID, tempMesh);
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

  // Calculate the total number of elements and vertices that were generated from splitting the current element
  out_cellData.m_numElems = (int) newElements.size();
  out_cellData.m_numVerts = (int) newVertices.size();

  generateTopologyData(newElements, newVertices, out_cellData);
  generateVertexPositionsFromQuad(newVertices, tempMesh, verticesClippingTValue, upperLeftVertex, lowerLeftVertex, lowerRightVertex, upperRightVertex, out_cellData);
  generateVertexVolumeFractionsFromQuad(newVertices, tempMesh, verticesClippingTValue, upperLeftVertex, lowerLeftVertex, lowerRightVertex, upperRightVertex, out_cellData);

  // Determine and store the dominant material of this element
  for (auto itr = newElements.begin(); itr != newElements.end(); itr++)
  {
    int dominantMaterial = determineDominantMaterial(Shape::Quad, itr->second, matOneID, matTwoID, out_cellData.m_mapData.m_vertexVolumeFractions);
    out_cellData.m_mapData.m_elementDominantMaterials.push_back(dominantMaterial);
  }

  // Determine and store the parent of this element
  out_cellData.m_mapData.m_elementParents.resize(newElements.size());
  for (auto itr = newElements.begin(); itr != newElements.end(); itr++)
  {
    out_cellData.m_mapData.m_elementParents[itr->first] = tempMesh.m_elementParentIDs[eID];
  }

  // Check that each element is dominated by a material that is actually present in the original parent cell
  for (auto itr = newElements.begin(); itr != newElements.end(); itr++)
  {
    int parentElementID = out_cellData.m_mapData.m_elementParents[itr->first];
    int currentDominantMaterial = out_cellData.m_mapData.m_elementDominantMaterials[itr->first];

    if (m_originalMesh.m_materialVolumeFractionsElement[currentDominantMaterial][parentElementID] == 0)
    {
      // This material is not present in the original element from which the current element comes from
      if (currentDominantMaterial == matOneID)
        out_cellData.m_mapData.m_elementDominantMaterials[itr->first] = matTwoID;
      else
        out_cellData.m_mapData.m_elementDominantMaterials[itr->first] = matOneID;
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

  for (unsigned int i = 0; i < out_cellData.m_topology.m_evInds.size(); ++i)
  {
    out_cellData.m_topology.m_evInds[i] -= evIndexSubtract[ out_cellData.m_topology.m_evInds[i] ];
  }
}

//--------------------------------------------------------------------------------

unsigned int InterfaceReconstructor::determineQuadClippingCase(mir::MIRMesh& tempMesh, const int matOneID, const int matTwoID, const int upperLeftVertex, const int lowerLeftVertex, const int lowerRightVertex, const int upperRightVertex)
{
    // Determine the dominant color at each vertex
    int upperLeftColor = matOneID, lowerLeftColor = matOneID, lowerRightColor = matOneID, upperRightColor = matOneID;

    if (matOneID == NULL_MAT)
      upperLeftColor = lowerLeftColor = lowerRightColor = upperRightColor = matTwoID;
    if (matTwoID == NULL_MAT)
      upperLeftColor = lowerLeftColor = lowerRightColor = upperRightColor = matOneID;
    if (matOneID != NULL_MAT && matTwoID != NULL_MAT)
    {
      upperLeftColor = tempMesh.m_materialVolumeFractionsVertex[matOneID][upperLeftVertex] > tempMesh.m_materialVolumeFractionsVertex[matTwoID][upperLeftVertex] ? matOneID : matTwoID;
      lowerLeftColor = tempMesh.m_materialVolumeFractionsVertex[matOneID][lowerLeftVertex] > tempMesh.m_materialVolumeFractionsVertex[matTwoID][lowerLeftVertex] ? matOneID : matTwoID;
      lowerRightColor = tempMesh.m_materialVolumeFractionsVertex[matOneID][lowerRightVertex] > tempMesh.m_materialVolumeFractionsVertex[matTwoID][lowerRightVertex] ? matOneID : matTwoID;
      upperRightColor = tempMesh.m_materialVolumeFractionsVertex[matOneID][upperRightVertex] > tempMesh.m_materialVolumeFractionsVertex[matTwoID][upperRightVertex] ? matOneID : matTwoID; 
    }
    
    // Create the index into the quad clipping lookup table using the dominant colors at each vertex
    unsigned int caseIndex = 0;
    if (upperLeftColor == matOneID)  caseIndex |= 8;
    if (lowerLeftColor == matOneID)  caseIndex |= 4;
    if (lowerRightColor == matOneID) caseIndex |= 2;
    if (upperRightColor == matOneID) caseIndex |= 1;

    return caseIndex;
}

//--------------------------------------------------------------------------------

mir::Point2 InterfaceReconstructor::interpolateVertexPosition(const mir::Point2& vertexOnePos, const mir::Point2& vertexTwoPos, const float t)
{
  mir::Point2 interpolatedPoint;
  interpolatedPoint.m_x = (1 - t) * vertexOnePos.m_x + t * vertexTwoPos.m_x;
  interpolatedPoint.m_y = (1 - t) * vertexOnePos.m_y + t * vertexTwoPos.m_y;
  return interpolatedPoint;
}

//--------------------------------------------------------------------------------

axom::float64 InterfaceReconstructor::lerpFloat(const axom::float64 f0, const axom::float64 f1, const axom::float64 t)
{
  return (1 - t) * f0 + t * f1;
}

//--------------------------------------------------------------------------------

axom::float64 InterfaceReconstructor::computeClippingPointOnEdge(const int vertexOneID, const int vertexTwoID, const int matOneID, const int matTwoID, mir::MIRMesh& tempMesh)
{
  axom::float64 ret = 0.0;

  axom::float64 vfMatOneVertexOne = tempMesh.m_materialVolumeFractionsVertex[matOneID][vertexOneID];
  axom::float64 vfMatTwoVertexOne = tempMesh.m_materialVolumeFractionsVertex[matTwoID][vertexOneID];
  axom::float64 vfMatOneVertexTwo = tempMesh.m_materialVolumeFractionsVertex[matOneID][vertexTwoID];
  axom::float64 vfMatTwoVertexTwo = tempMesh.m_materialVolumeFractionsVertex[matTwoID][vertexTwoID];  

  if (matOneID == NULL_MAT)
    vfMatOneVertexOne = vfMatOneVertexTwo = 0.0;

  if (matTwoID == NULL_MAT)
    vfMatTwoVertexOne = vfMatTwoVertexTwo = 0.0;

  axom::float64 numerator = vfMatTwoVertexOne - vfMatOneVertexOne;
  axom::float64 denominator = -vfMatOneVertexOne + vfMatOneVertexTwo + vfMatTwoVertexOne - vfMatTwoVertexTwo;
        
  if (denominator != 0.0)
    ret = numerator / denominator;

  if (ret > 1.0 || ret < 0.0)
  {
    // This shouldn't happen...
    printf("    OUT OF BOUNDS T VALUE: %f\n", ret);

    // Clamp the t value
    ret = fmin(1.0, ret);
    ret = fmax(0.0, ret);
  }

  return ret;
}

//--------------------------------------------------------------------------------

void InterfaceReconstructor::generateTopologyData(const std::map<int, std::vector<int> >& newElements, const std::map<int, std::vector<int> >& newVertices, CellData& out_cellData)
{
    // Store the evInds and evBegins data in the output vectors
    int currentEVBeginIndex = 0;
    for (auto itr = newElements.begin(); itr != newElements.end(); itr++)
    {
      // Push the start index of the next element
      out_cellData.m_topology.m_evBegins.push_back(currentEVBeginIndex);

      // Push the next element's vertices
      for (unsigned int vIndex = 0; vIndex < itr->second.size(); ++vIndex)
      {
        out_cellData.m_topology.m_evInds.push_back(itr->second[vIndex]);
        ++currentEVBeginIndex;
      }
    }

    // Push the index that occurs after the last vertex
    out_cellData.m_topology.m_evBegins.push_back(currentEVBeginIndex);
    
    // Store the veInds and veBegins data in the output vectors
    int currentVEBeginIndex = 0;
    for (auto itr = newVertices.begin(); itr != newVertices.end(); itr++)
    {
      // Push the start index of the vertex's elements
      out_cellData.m_topology.m_veBegins.push_back(currentVEBeginIndex);

      // Push the next vertex's elements
      for (unsigned int eIndex = 0; eIndex < itr->second.size(); eIndex++)
      {
        out_cellData.m_topology.m_veInds.push_back(itr->second[eIndex]);
        ++currentVEBeginIndex;
      }
    }
    
    // Push the index that occurs after the last element
    out_cellData.m_topology.m_veBegins.push_back(currentVEBeginIndex);
}

//--------------------------------------------------------------------------------

void InterfaceReconstructor::generateVertexPositionsFromQuad(const std::map<int, std::vector<int> >& newVertices, mir::MIRMesh& tempMesh, axom::float64* verticesClippingTValue, const int upperLeftVertex, const int lowerLeftVertex, const int lowerRightVertex, const int upperRightVertex, CellData& out_cellData)
{
  std::map<int, mir::Point2> newPoints;           // hashmap of the new vertices' positions | Note: maps are ordered sets

  // Find the position coordinates of the vertices and then store the positions of the vertices in the return vector
  for (auto itr = newVertices.begin(); itr != newVertices.end(); itr++)
  {
    int vID = itr->first;
    if (vID == 0)
      newPoints[vID] = tempMesh.m_vertexPositions[upperLeftVertex];
    if (vID == 1)
      newPoints[vID] = tempMesh.m_vertexPositions[lowerLeftVertex];
    if (vID == 2)
      newPoints[vID] = tempMesh.m_vertexPositions[lowerRightVertex];
    if (vID == 3)
      newPoints[vID] = tempMesh.m_vertexPositions[upperRightVertex];
    if (vID == 4)
      newPoints[vID] = interpolateVertexPosition(tempMesh.m_vertexPositions[upperLeftVertex], tempMesh.m_vertexPositions[lowerLeftVertex], verticesClippingTValue[vID]);
    if (vID == 5)
      newPoints[vID] = interpolateVertexPosition(tempMesh.m_vertexPositions[lowerLeftVertex], tempMesh.m_vertexPositions[lowerRightVertex], verticesClippingTValue[vID]);
    if (vID == 6)
      newPoints[vID] = interpolateVertexPosition(tempMesh.m_vertexPositions[lowerRightVertex], tempMesh.m_vertexPositions[upperRightVertex], verticesClippingTValue[vID]);
    if (vID == 7)
      newPoints[vID] = interpolateVertexPosition(tempMesh.m_vertexPositions[upperRightVertex], tempMesh.m_vertexPositions[upperLeftVertex], verticesClippingTValue[vID]);
  } 

  // Store the positions of the vertices in the return vector
  for (auto itr = newPoints.begin(); itr != newPoints.end(); itr++)
  {
    out_cellData.m_mapData.m_vertexPositions.push_back(itr->second);
  }
}

//--------------------------------------------------------------------------------

void InterfaceReconstructor::generateVertexVolumeFractionsFromQuad(std::map<int, std::vector<int> >& newVertices, mir::MIRMesh& tempMesh, axom::float64* verticesClippingTValue, const int upperLeftVertex, const int lowerLeftVertex, const int lowerRightVertex, const int upperRightVertex, CellData& out_cellData)
{
   // Calculate the vertex fractions at each vertex (use t value!)
    // Make sure the output volume fractions containers are the proper size
    out_cellData.m_mapData.m_vertexVolumeFractions.resize(tempMesh.m_numMaterials);
    
    for (auto itr = newVertices.begin(); itr != newVertices.end(); itr++)
    {
      int vID = itr->first;
      for (int matID = 0; matID < tempMesh.m_numMaterials; ++matID)
      {
        if (vID == 0)
          out_cellData.m_mapData.m_vertexVolumeFractions[matID].push_back(tempMesh.m_materialVolumeFractionsVertex[matID][upperLeftVertex]);
        if (vID == 1)
          out_cellData.m_mapData.m_vertexVolumeFractions[matID].push_back(tempMesh.m_materialVolumeFractionsVertex[matID][lowerLeftVertex]);
        if (vID == 2)
          out_cellData.m_mapData.m_vertexVolumeFractions[matID].push_back(tempMesh.m_materialVolumeFractionsVertex[matID][lowerRightVertex]);
        if (vID == 3)
          out_cellData.m_mapData.m_vertexVolumeFractions[matID].push_back(tempMesh.m_materialVolumeFractionsVertex[matID][upperRightVertex]);
        if (vID == 4)
          out_cellData.m_mapData.m_vertexVolumeFractions[matID].push_back(lerpFloat(tempMesh.m_materialVolumeFractionsVertex[matID][upperLeftVertex], tempMesh.m_materialVolumeFractionsVertex[matID][lowerLeftVertex], verticesClippingTValue[vID]));
        if (vID == 5)
          out_cellData.m_mapData.m_vertexVolumeFractions[matID].push_back(lerpFloat(tempMesh.m_materialVolumeFractionsVertex[matID][lowerLeftVertex], tempMesh.m_materialVolumeFractionsVertex[matID][lowerRightVertex], verticesClippingTValue[vID]));
        if (vID == 6)
          out_cellData.m_mapData.m_vertexVolumeFractions[matID].push_back(lerpFloat(tempMesh.m_materialVolumeFractionsVertex[matID][lowerRightVertex], tempMesh.m_materialVolumeFractionsVertex[matID][upperRightVertex], verticesClippingTValue[vID]));
        if (vID == 7)
          out_cellData.m_mapData.m_vertexVolumeFractions[matID].push_back(lerpFloat(tempMesh.m_materialVolumeFractionsVertex[matID][upperRightVertex], tempMesh.m_materialVolumeFractionsVertex[matID][upperLeftVertex], verticesClippingTValue[vID]));
      }
    }
}

//--------------------------------------------------------------------------------

void InterfaceReconstructor::computeTriangleClippingPoints(const int eID, const int matOneID, const int matTwoID, mir::MIRMesh& tempMesh, CellData& out_cellData)
{
  // Determine the clipping case
  auto elementVertices = tempMesh.m_bdry[eID];

  int upperVertex = elementVertices[0];
  int lowerLeftVertex = elementVertices[1];
  int lowerRightVertex = elementVertices[2];

  unsigned int caseIndex = determineTriangleClippingCase(tempMesh, matOneID, matTwoID, upperVertex, lowerLeftVertex, lowerRightVertex);

  // Generate new elements
  std::map<int, std::vector<int> > newElements;             // hashmap of the new elements' vertices | Note: maps are ordered sets
  std::map<int, std::vector<int> > newVertices;             // hashmap of the new vertices' elements | Note: maps are ordered sets
  int verticesPresent[6] = {0,0,0,0,0,0};                   // Array of flags denoting whether the vertex is present in the current case or not
  axom::float64 verticesClippingTValue[6] = {0,0,0,0,0,0};  // Array of t values that denote the percent value of where the edge should be clipped

  // Create the new polygons based on the clipping case
  int currentElementIndex = 0;    // the next available element index
  int i = 0;
  int numVertices = triangleClipTable[caseIndex][i];

  // for each new element in the current clipping case
  while (numVertices != -1)   
  {
    // for each vertex of the new element
    for (int j = 0; j < numVertices; ++j)
    {
      // Find the id of the next vertex of the new element
      int vID = triangleClipTable[caseIndex][i + (j+1)];

      // Associate the vertex and element together
      newElements[currentElementIndex].push_back(vID);
      newVertices[vID].push_back(currentElementIndex);
      verticesPresent[vID] = 1;

      // Find t using linear interpolation for any vertex that is not one of the original 3 vertices
      if(vID == 3)
      {
        verticesClippingTValue[vID] = computeClippingPointOnEdge(upperVertex, lowerLeftVertex, matOneID, matTwoID, tempMesh);
      }
      else if(vID == 4)
      {
        verticesClippingTValue[vID] = computeClippingPointOnEdge(lowerLeftVertex, lowerRightVertex, matOneID, matTwoID, tempMesh);
      }
      else if(vID == 5)
      {
        verticesClippingTValue[vID] = computeClippingPointOnEdge(lowerRightVertex, upperVertex, matOneID, matTwoID, tempMesh);
      }
    }

    // Increment the element index counter, marking the current element as being finished processed
    currentElementIndex++;

    // Increase index into lookup table to the next element
    i += (numVertices + 1);
    numVertices = triangleClipTable[caseIndex][i];
  }

  // Calculate the total number of elements and vertices that were generated from splitting the current element
  out_cellData.m_numElems = (int) newElements.size();
  out_cellData.m_numVerts = (int) newVertices.size();

  generateTopologyData(newElements, newVertices, out_cellData);
  generateVertexPositionsFromTriangle(newVertices, tempMesh, verticesClippingTValue, upperVertex, lowerLeftVertex, lowerRightVertex, out_cellData);
  generateVertexVolumeFractionsFromTriangle(newVertices, tempMesh, verticesClippingTValue, upperVertex, lowerLeftVertex, lowerRightVertex, out_cellData);

  // Determine and store the dominant material of this element
  for (auto itr = newElements.begin(); itr != newElements.end(); itr++)
  {
    int dominantMaterial = determineDominantMaterial(Shape::Triangle, itr->second, matOneID, matTwoID, out_cellData.m_mapData.m_vertexVolumeFractions);
    out_cellData.m_mapData.m_elementDominantMaterials.push_back(dominantMaterial);
  }

  // Determine and store the parent of this element
  out_cellData.m_mapData.m_elementParents.resize(newElements.size());
  for (auto itr = newElements.begin(); itr != newElements.end(); itr++)
  {
    out_cellData.m_mapData.m_elementParents[itr->first] = tempMesh.m_elementParentIDs[eID];
  }

  // Check that each element is dominated by a material that is actually present in the original parent cell
  for (auto itr = newElements.begin(); itr != newElements.end(); itr++)
  {
    int parentElementID = out_cellData.m_mapData.m_elementParents[itr->first];
    int currentDominantMaterial = out_cellData.m_mapData.m_elementDominantMaterials[itr->first];

    if (m_originalMesh.m_materialVolumeFractionsElement[currentDominantMaterial][parentElementID] == 0)
    {
      // This material is not present in the original element from which the current element comes from
      if (currentDominantMaterial == matOneID)
        out_cellData.m_mapData.m_elementDominantMaterials[itr->first] = matTwoID;
      else
        out_cellData.m_mapData.m_elementDominantMaterials[itr->first] = matOneID;
    }
  }  

  // Modify the evIndex values to account for the fact that perhaps not all 6 possible vertices are present
  int evIndexSubtract[6];
  for (int i = 0; i < 6; ++i)
  {
    if (verticesPresent[i] == 1)
      evIndexSubtract[i] = 0;
    else
      evIndexSubtract[i] = 1;
  }

  for (int i = 1; i < 6; ++i)
    evIndexSubtract[i] += evIndexSubtract[i - 1];

  for (unsigned int i = 0; i < out_cellData.m_topology.m_evInds.size(); ++i)
    out_cellData.m_topology.m_evInds[i] -= evIndexSubtract[ out_cellData.m_topology.m_evInds[i] ];

}

//--------------------------------------------------------------------------------

unsigned int InterfaceReconstructor::determineTriangleClippingCase(mir::MIRMesh& tempMesh, const int matOneID, const int matTwoID, const int upperVertex, const int lowerLeftVertex, const int lowerRightVertex)
{
 // Determine the dominant color at each vertex
    int upperColor = matOneID, lowerLeftColor = matOneID, lowerRightColor = matOneID;

    if (matOneID == NULL_MAT)
      upperColor = lowerLeftColor = lowerRightColor = matTwoID;
    if (matTwoID == NULL_MAT)
      upperColor = lowerLeftColor = lowerRightColor = matOneID;
    if (matOneID != NULL_MAT && matTwoID != NULL_MAT)
    {
      upperColor = tempMesh.m_materialVolumeFractionsVertex[matOneID][upperVertex] > tempMesh.m_materialVolumeFractionsVertex[matTwoID][upperVertex] ? matOneID : matTwoID;
      lowerLeftColor = tempMesh.m_materialVolumeFractionsVertex[matOneID][lowerLeftVertex] > tempMesh.m_materialVolumeFractionsVertex[matTwoID][lowerLeftVertex] ? matOneID : matTwoID;
      lowerRightColor = tempMesh.m_materialVolumeFractionsVertex[matOneID][lowerRightVertex] > tempMesh.m_materialVolumeFractionsVertex[matTwoID][lowerRightVertex] ? matOneID : matTwoID;
    }
    
    // Create the index into the quad clipping lookup table using the dominant colors at each vertex
    unsigned int caseIndex = 0;
    if (upperColor == matOneID)  caseIndex |= 4;
    if (lowerLeftColor == matOneID)  caseIndex |= 2;
    if (lowerRightColor == matOneID) caseIndex |= 1;

    return caseIndex;
}

//--------------------------------------------------------------------------------

void InterfaceReconstructor::generateVertexPositionsFromTriangle(const std::map<int, std::vector<int> >& newVertices, mir::MIRMesh& tempMesh, axom::float64* verticesClippingTValue, const int upperVertex, const int lowerLeftVertex, const int lowerRightVertex, CellData& out_cellData)
{
  std::map<int, mir::Point2> newPoints;           // hashmap of the new vertices' positions | Note: maps are ordered sets

  // Find the position coordinates of the vertices and then store the positions of the vertices in the return vector
  for (auto itr = newVertices.begin(); itr != newVertices.end(); itr++)
  {
    int vID = itr->first;
    if (vID == 0)
      newPoints[vID] = tempMesh.m_vertexPositions[upperVertex];
    if (vID == 1)
      newPoints[vID] = tempMesh.m_vertexPositions[lowerLeftVertex];
    if (vID == 2)
      newPoints[vID] = tempMesh.m_vertexPositions[lowerRightVertex];
    if (vID == 3)
      newPoints[vID] = interpolateVertexPosition(tempMesh.m_vertexPositions[upperVertex], tempMesh.m_vertexPositions[lowerLeftVertex], verticesClippingTValue[vID]);
    if (vID == 4)
      newPoints[vID] = interpolateVertexPosition(tempMesh.m_vertexPositions[lowerLeftVertex], tempMesh.m_vertexPositions[lowerRightVertex], verticesClippingTValue[vID]);
    if (vID == 5)
      newPoints[vID] = interpolateVertexPosition(tempMesh.m_vertexPositions[lowerRightVertex], tempMesh.m_vertexPositions[upperVertex], verticesClippingTValue[vID]);
  } 

  // Store the positions of the vertices in the return vector
  for (auto itr = newPoints.begin(); itr != newPoints.end(); itr++)
  {
    out_cellData.m_mapData.m_vertexPositions.push_back(itr->second);
  }
}

//--------------------------------------------------------------------------------

void InterfaceReconstructor::generateVertexVolumeFractionsFromTriangle(const std::map<int, std::vector<int> >& newVertices, mir::MIRMesh& tempMesh, axom::float64* verticesClippingTValue, const int upperVertex, const int lowerLeftVertex, const int lowerRightVertex, CellData& out_cellData)
{
    // Calculate the vertex fractions at each vertex (use t value!)
    // Make sure the output volume fractions containers are the proper size
    out_cellData.m_mapData.m_vertexVolumeFractions.resize(tempMesh.m_numMaterials);
    
    for (auto itr = newVertices.begin(); itr != newVertices.end(); itr++)
    {
      int vID = itr->first;
      for (int matID = 0; matID < tempMesh.m_numMaterials; ++matID)
      {
        if (vID == 0)
          out_cellData.m_mapData.m_vertexVolumeFractions[matID].push_back(tempMesh.m_materialVolumeFractionsVertex[matID][upperVertex]);
        if (vID == 1)
          out_cellData.m_mapData.m_vertexVolumeFractions[matID].push_back(tempMesh.m_materialVolumeFractionsVertex[matID][lowerLeftVertex]);
        if (vID == 2)
          out_cellData.m_mapData.m_vertexVolumeFractions[matID].push_back(tempMesh.m_materialVolumeFractionsVertex[matID][lowerRightVertex]);
        if (vID == 3)
          out_cellData.m_mapData.m_vertexVolumeFractions[matID].push_back(lerpFloat(tempMesh.m_materialVolumeFractionsVertex[matID][upperVertex], tempMesh.m_materialVolumeFractionsVertex[matID][lowerLeftVertex], verticesClippingTValue[vID]));
        if (vID == 4)
          out_cellData.m_mapData.m_vertexVolumeFractions[matID].push_back(lerpFloat(tempMesh.m_materialVolumeFractionsVertex[matID][lowerLeftVertex], tempMesh.m_materialVolumeFractionsVertex[matID][lowerRightVertex], verticesClippingTValue[vID]));
        if (vID == 5)
          out_cellData.m_mapData.m_vertexVolumeFractions[matID].push_back(lerpFloat(tempMesh.m_materialVolumeFractionsVertex[matID][lowerRightVertex], tempMesh.m_materialVolumeFractionsVertex[matID][upperVertex], verticesClippingTValue[vID]));    
      }
    }
}

//--------------------------------------------------------------------------------

int InterfaceReconstructor::determineDominantMaterial(const Shape elementShape, const std::vector<int>& vertexIDs, const int matOne, const int matTwo, const std::vector<std::vector<axom::float64> >& vertexVF)
{
  int dominantMaterial = matOne;

  axom::float64 matOneVF = -1.0;
  axom::float64 matTwoVF = -1.0;

  if (elementShape == Shape::Triangle)
  {
    for (unsigned long it = 0; it < vertexIDs.size(); ++it)
    {
      int vID = vertexIDs[it];
      if (vID == 0 || vID == 1 || vID == 2)
      {
        if (matOne != NULL_MAT)
        {
          matOneVF = vertexVF[matOne][vID];
        }
        if (matTwo != NULL_MAT)
        {
          matTwoVF = vertexVF[matTwo][vID];
        }

        dominantMaterial = (matOneVF > matTwoVF) ? matOne : matTwo;
      }
    }
  }

  if (elementShape == Shape::Quad)
  {
    for (unsigned long it = 0; it < vertexIDs.size(); ++it)
    {
      int vID = vertexIDs[it];
      if (vID == 0 || vID == 1 || vID == 2 || vID == 3)
      {
        if (matOne != NULL_MAT)
        {
          matOneVF = vertexVF[matOne][vID];
        }
        if (matTwo != NULL_MAT)
        {
          matTwoVF = vertexVF[matTwo][vID];
        }

        dominantMaterial = (matOneVF > matTwoVF) ? matOne : matTwo;
      }
    }
  }

  return dominantMaterial;
}

//--------------------------------------------------------------------------------

}
}
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

/// Destructor
InterfaceReconstructor::~InterfaceReconstructor()
{

}

//--------------------------------------------------------------------------------

/// Reconstructs the material interface and returns the resulting output mesh.
mir::MIRMesh InterfaceReconstructor::computeReconstructedInterface(mir::MIRMesh* inputMesh)
{
  // Store a pointer to the original mesh
  originalMesh = inputMesh;

  // Initialize the final mesh to be the same as the input mesh
  mir::MIRMesh finalMesh(originalMesh);

  // For each material in the mesh, split the mesh on the current material and generate a new mesh (input into next iteration)
  for (int matID = 0; matID < originalMesh->numMaterials; ++matID)
  {
    // Copy the mesh to be split
    mir::MIRMesh intermediateMesh(&finalMesh);

    // Update the materials upon which the split will occur
    int matOne = matID;

    // Create an array to store the output of each element being split.
    CellData temp_cellData[intermediateMesh.elems.size()];

    // Process/split each element
    for (int eID = 0; eID < intermediateMesh.elems.size(); ++eID)
    {
      // Update the materials upon which the split will occur (should be the currently dominant material and the next material in the list that hasn't yet been split on)
      int currentDominantMat = intermediateMesh.elementDominantMaterials[eID];
      computeClippingPoints(eID, currentDominantMat, matOne, &intermediateMesh, temp_cellData[eID]);
    }

    // Merge each of the cells into the first CellData struct
    for (int eID = 1; eID < intermediateMesh.elems.size(); ++eID)
      temp_cellData[0].mergeCell(temp_cellData[eID]);

    mir::VertSet combined_verts(temp_cellData[0].numVerts);
    mir::ElemSet combined_elems(temp_cellData[0].numElems);

    // Create the final, processed mesh
    mir::MIRMesh processedMesh;
    processedMesh.InitializeMesh(temp_cellData[0].evInds, temp_cellData[0].evBegins, temp_cellData[0].veInds, temp_cellData[0].veBegins, combined_verts, combined_elems, intermediateMesh.numMaterials);
    processedMesh.constructMeshRelations();
    processedMesh.constructMeshVolumeFractionsVertex(temp_cellData[0].vertexVolumeFractions);
    processedMesh.constructVertexPositionMap(temp_cellData[0].vertexPositions.data());  
    processedMesh.constructElementParentMap(temp_cellData[0].elementParents.data());
    processedMesh.constructElementDominantMaterialMap(temp_cellData[0].elementDominantMaterials);

    // Store the current mesh to be passed into the next iteration
    finalMesh = processedMesh;
    finalMesh.constructMeshRelations();
  }

  // TODO: Run post-processing on mesh to meld the mesh vertices that are within a certain epsilon

  return finalMesh;
}

/// Reconstructs the material interface using an iterative optimization to improve the volume fractions of the resulting mesh.
/// Based on the Meredith and Childs 2010 paper.
mir::MIRMesh InterfaceReconstructor::computeReconstructedInterfaceIterative(mir::MIRMesh* inputMesh, int numIterations, axom::float64 percent)
{
  int numElems = inputMesh->elems.size();
  int numMaterials = inputMesh->numMaterials;

  // Make a copy of the original input mesh
  mir::MIRMesh meshToImprove(inputMesh);

  // Calculate the reconstruction on the unmodified, input mesh
  mir::MIRMesh resultingMesh = computeReconstructedInterface(&meshToImprove);

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
        axom::float64 difference = inputMesh->materialVolumeFractionsElement[matID][eID] - resultingElementVF[matID][eID];
        improvedElementVF[matID].push_back(   inputMesh->materialVolumeFractionsElement[matID][eID] + ( difference * percent ) );
      }
    }

    // Reconstruct the element AND vertex VFs
    meshToImprove.constructMeshVolumeFractionsMaps(improvedElementVF);

    // Calculate the reconstruction on the modified, input mesh
    resultingMesh = computeReconstructedInterface(&meshToImprove);
  }

  return resultingMesh;
}

//--------------------------------------------------------------------------------

/// Computes the points where the element should be clipped based on the volume fractions of mat one and two.
void InterfaceReconstructor::computeClippingPoints(const int eID, const int matOneID, const int matTwoID, mir::MIRMesh* tempMesh, CellData& out_cellData)
{
  // Find the vertices associated with each element
  auto elementVertices = tempMesh->bdry[eID];

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

void InterfaceReconstructor::computeQuadClippingPoints(const int eID, const int matOneID, const int matTwoID, mir::MIRMesh* tempMesh, CellData& out_cellData)
{
  // printf("  Processing QUAD element %d: ", eID);

  // Determine the clipping case
  auto elementVertices = tempMesh->bdry[eID];

  int upperLeftVertex = elementVertices[0];
  int lowerLeftVertex = elementVertices[1];
  int lowerRightVertex = elementVertices[2];
  int upperRightVertex = elementVertices[3];

  unsigned int caseIndex = determineQuadClippingCase(tempMesh, matOneID, matTwoID, upperLeftVertex, lowerLeftVertex, lowerRightVertex, upperRightVertex);
  // printf("caseIndex: %d\n", caseIndex);

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

  /****************************************************************
   *                CALCULATE THE NEW CELLS' DATA
   ****************************************************************/
  // Calculate the total number of elements and vertices that were generated from splitting the current element
  out_cellData.numElems = (int) newElements.size();
  out_cellData.numVerts = (int) newVertices.size();

  generateTopologyData(newElements, newVertices, out_cellData);
  generateVertexPositionsFromQuad(newVertices, tempMesh, verticesClippingTValue, upperLeftVertex, lowerLeftVertex, lowerRightVertex, upperRightVertex, out_cellData);
  generateVertexVolumeFractionsFromQuad(newVertices, tempMesh, verticesClippingTValue, upperLeftVertex, lowerLeftVertex, lowerRightVertex, upperRightVertex, out_cellData);

  // Determine and store the dominant material of this element
  for (auto itr = newElements.begin(); itr != newElements.end(); itr++)
  {
    int currentDominantMat = NULL_MAT;
    for (unsigned long it = 0; it < itr->second.size(); ++it)
    {
      // int temp_vID = newElements[currentElementIndex][it];
      int temp_vID = itr->second[it];
  
      // Find the vertex that is one of the four original vertices of the quad element
      if (temp_vID == 0 || temp_vID == 1 || temp_vID == 2 || temp_vID == 3)
      {
        axom::float64 matOneVolumeFraction = -1.0, matTwoVolumeFraction = -1.0;
        if (matOneID != NULL_MAT)
          matOneVolumeFraction = out_cellData.vertexVolumeFractions[matOneID][temp_vID];
        if (matTwoID != NULL_MAT)
          matTwoVolumeFraction = out_cellData.vertexVolumeFractions[matTwoID][temp_vID];

        currentDominantMat = (matOneVolumeFraction > matTwoVolumeFraction) ? matOneID : matTwoID;
      }
    }
    out_cellData.elementDominantMaterials.push_back(currentDominantMat);  
  }

  // Determine and store the parent of this element
  out_cellData.elementParents.resize(newElements.size());
  for (auto itr = newElements.begin(); itr != newElements.end(); itr++)
  {
    out_cellData.elementParents[itr->first] = tempMesh->elementParentIDs[eID];
  }

  // Check that each element is dominated by a material that is actually present in the original parent cell
  for (auto itr = newElements.begin(); itr != newElements.end(); itr++)
  {
    int parentElementID = out_cellData.elementParents[itr->first];
    int currentDominantMaterial = out_cellData.elementDominantMaterials[itr->first];

    if (originalMesh->materialVolumeFractionsElement[currentDominantMaterial][parentElementID] == 0)
    {
      // This material is not present in the original element from which the current element comes from
      if (currentDominantMaterial == matOneID)
        out_cellData.elementDominantMaterials[itr->first] = matTwoID;
      else
        out_cellData.elementDominantMaterials[itr->first] = matOneID;
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

  for (unsigned int i = 0; i < out_cellData.evInds.size(); ++i)
  {
    out_cellData.evInds[i] -= evIndexSubtract[ out_cellData.evInds[i] ];
  }
}

//--------------------------------------------------------------------------------

/// Finds the bit map representing the clipping case for a quad.
unsigned int InterfaceReconstructor::determineQuadClippingCase(mir::MIRMesh* tempMesh, const int matOneID, const int matTwoID, const int upperLeftVertex, const int lowerLeftVertex, const int lowerRightVertex, const int upperRightVertex)
{
    // Determine the dominant color at each vertex
    int upperLeftColor = matOneID, lowerLeftColor = matOneID, lowerRightColor = matOneID, upperRightColor = matOneID;

    if (matOneID == NULL_MAT)
      upperLeftColor = lowerLeftColor = lowerRightColor = upperRightColor = matTwoID;
    if (matTwoID == NULL_MAT)
      upperLeftColor = lowerLeftColor = lowerRightColor = upperRightColor = matOneID;
    if (matOneID != NULL_MAT && matTwoID != NULL_MAT)
    {
      upperLeftColor = tempMesh->materialVolumeFractionsVertex[matOneID][upperLeftVertex] > tempMesh->materialVolumeFractionsVertex[matTwoID][upperLeftVertex] ? matOneID : matTwoID;
      lowerLeftColor = tempMesh->materialVolumeFractionsVertex[matOneID][lowerLeftVertex] > tempMesh->materialVolumeFractionsVertex[matTwoID][lowerLeftVertex] ? matOneID : matTwoID;
      lowerRightColor = tempMesh->materialVolumeFractionsVertex[matOneID][lowerRightVertex] > tempMesh->materialVolumeFractionsVertex[matTwoID][lowerRightVertex] ? matOneID : matTwoID;
      upperRightColor = tempMesh->materialVolumeFractionsVertex[matOneID][upperRightVertex] > tempMesh->materialVolumeFractionsVertex[matTwoID][upperRightVertex] ? matOneID : matTwoID; 
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

  axom::float64 vfMatOneVertexOne = tempMesh->materialVolumeFractionsVertex[matOneID][vertexOneID];   // Changed these all from originalMesh->materialVolumeFractionsVertex[][] to tempMesh->...
  axom::float64 vfMatTwoVertexOne = tempMesh->materialVolumeFractionsVertex[matTwoID][vertexOneID];
  axom::float64 vfMatOneVertexTwo = tempMesh->materialVolumeFractionsVertex[matOneID][vertexTwoID];
  axom::float64 vfMatTwoVertexTwo = tempMesh->materialVolumeFractionsVertex[matTwoID][vertexTwoID];  

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

/// Generate the topology of the new elements (evInds, evBegins, veInds, veBegins) resulting from a split.
void InterfaceReconstructor::generateTopologyData(std::map<int, std::vector<int> > newElements, std::map<int, std::vector<int> > newVertices, CellData& out_cellData)
{
    // Store the evInds and evBegins data in the output vectors
    int currentEVBeginIndex = 0;
    for (auto itr = newElements.begin(); itr != newElements.end(); itr++)
    {
      // Push the start index of the next element
      out_cellData.evBegins.push_back(currentEVBeginIndex);

      // Push the next element's vertices
      for (unsigned int vIndex = 0; vIndex < itr->second.size(); ++vIndex)
      {
        out_cellData.evInds.push_back(itr->second[vIndex]);
        ++currentEVBeginIndex;
      }
    }

    // Push the index that occurs after the last vertex
    out_cellData.evBegins.push_back(currentEVBeginIndex);
    
    // Store the veInds and veBegins data in the output vectors
    int currentVEBeginIndex = 0;
    for (auto itr = newVertices.begin(); itr != newVertices.end(); itr++)
    {
      // Push the start index of the vertex's elements
      out_cellData.veBegins.push_back(currentVEBeginIndex);

      // Push the next vertex's elements
      for (unsigned int eIndex = 0; eIndex < itr->second.size(); eIndex++)
      {
        out_cellData.veInds.push_back(itr->second[eIndex]);
        ++currentVEBeginIndex;
      }
    }
    
    // Push the index that occurs after the last element
    out_cellData.veBegins.push_back(currentVEBeginIndex);
}

//--------------------------------------------------------------------------------

/// Generate the vertex position data for the new elements resulting from spltting a quad.
void InterfaceReconstructor::generateVertexPositionsFromQuad(std::map<int, std::vector<int> > newVertices, mir::MIRMesh* tempMesh, axom::float64* verticesClippingTValue, int upperLeftVertex, int lowerLeftVertex, int lowerRightVertex, int upperRightVertex, CellData& out_cellData)
{
  std::map<int, mir::Point2> newPoints;           // hashmap of the new vertices' positions | Note: maps are ordered sets

  // Find the position coordinates of the vertices and then store the positions of the vertices in the return vector
  for (auto itr = newVertices.begin(); itr != newVertices.end(); itr++)
  {
    int vID = itr->first;
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

  // Store the positions of the vertices in the return vector
  for (auto itr = newPoints.begin(); itr != newPoints.end(); itr++)
  {
    out_cellData.vertexPositions.push_back(itr->second);
  }
}

//--------------------------------------------------------------------------------

/// Generate the vertex volume fraction data for the new vertices resulting from splitting a quad. 
void InterfaceReconstructor::generateVertexVolumeFractionsFromQuad(std::map<int, std::vector<int> > newVertices, mir::MIRMesh* tempMesh, axom::float64* verticesClippingTValue, int upperLeftVertex, int lowerLeftVertex, int lowerRightVertex, int upperRightVertex, CellData& out_cellData)
{
   // Calculate the vertex fractions at each vertex (use t value!)
    // Make sure the output volume fractions containers are the proper size
    out_cellData.vertexVolumeFractions.resize(tempMesh->numMaterials);
    
    for (auto itr = newVertices.begin(); itr != newVertices.end(); itr++)
    {
      int vID = itr->first;
      for (int matID = 0; matID < tempMesh->numMaterials; ++matID)
      {
        if (vID == 0)
          out_cellData.vertexVolumeFractions[matID].push_back(tempMesh->materialVolumeFractionsVertex[matID][upperLeftVertex]);
        if (vID == 1)
          out_cellData.vertexVolumeFractions[matID].push_back(tempMesh->materialVolumeFractionsVertex[matID][lowerLeftVertex]);
        if (vID == 2)
          out_cellData.vertexVolumeFractions[matID].push_back(tempMesh->materialVolumeFractionsVertex[matID][lowerRightVertex]);
        if (vID == 3)
          out_cellData.vertexVolumeFractions[matID].push_back(tempMesh->materialVolumeFractionsVertex[matID][upperRightVertex]);
        if (vID == 4)
          out_cellData.vertexVolumeFractions[matID].push_back(lerpFloat(tempMesh->materialVolumeFractionsVertex[matID][upperLeftVertex], tempMesh->materialVolumeFractionsVertex[matID][lowerLeftVertex], verticesClippingTValue[vID]));
        if (vID == 5)
          out_cellData.vertexVolumeFractions[matID].push_back(lerpFloat(tempMesh->materialVolumeFractionsVertex[matID][lowerLeftVertex], tempMesh->materialVolumeFractionsVertex[matID][lowerRightVertex], verticesClippingTValue[vID]));
        if (vID == 6)
          out_cellData.vertexVolumeFractions[matID].push_back(lerpFloat(tempMesh->materialVolumeFractionsVertex[matID][lowerRightVertex], tempMesh->materialVolumeFractionsVertex[matID][upperRightVertex], verticesClippingTValue[vID]));
        if (vID == 7)
          out_cellData.vertexVolumeFractions[matID].push_back(lerpFloat(tempMesh->materialVolumeFractionsVertex[matID][upperRightVertex], tempMesh->materialVolumeFractionsVertex[matID][upperLeftVertex], verticesClippingTValue[vID]));
      }
    }
}

//--------------------------------------------------------------------------------

void InterfaceReconstructor::computeTriangleClippingPoints(const int eID, const int matOneID, const int matTwoID, mir::MIRMesh* tempMesh, CellData& out_cellData)
{
  // Determine the clipping case
  auto elementVertices = tempMesh->bdry[eID];

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

  /****************************************************************
   *                CALCULATE THE NEW CELLS' DATA
   ****************************************************************/
  // Calculate the total number of elements and vertices that were generated from splitting the current element
  out_cellData.numElems = (int) newElements.size();
  out_cellData.numVerts = (int) newVertices.size();

  generateTopologyData(newElements, newVertices, out_cellData);
  generateVertexPositionsFromTriangle(newVertices, tempMesh, verticesClippingTValue, upperVertex, lowerLeftVertex, lowerRightVertex, out_cellData);
  generateVertexVolumeFractionsFromTriangle(newVertices, tempMesh, verticesClippingTValue, upperVertex, lowerLeftVertex, lowerRightVertex, out_cellData);

  // Determine and store the dominant material of this element
  for (auto itr = newElements.begin(); itr != newElements.end(); itr++)
  {
    int currentDominantMat = NULL_MAT;
    for (unsigned long it = 0; it < itr->second.size(); ++it)
    {
      int temp_vID = itr->second[it];
  
      // Find the vertex that is one of the three original vertices of the triangle element
      if (temp_vID == 0 || temp_vID == 1 || temp_vID == 2)
      {
        axom::float64 matOneVolumeFraction = -1.0, matTwoVolumeFraction = -1.0;
        if (matOneID != NULL_MAT)
          matOneVolumeFraction = out_cellData.vertexVolumeFractions[matOneID][temp_vID];
        if (matTwoID != NULL_MAT)
          matTwoVolumeFraction = out_cellData.vertexVolumeFractions[matTwoID][temp_vID];

        currentDominantMat = (matOneVolumeFraction > matTwoVolumeFraction) ? matOneID : matTwoID;
      }
    }
    out_cellData.elementDominantMaterials.push_back(currentDominantMat);  
  }

  // Determine and store the parent of this element
  out_cellData.elementParents.resize(newElements.size());
  for (auto itr = newElements.begin(); itr != newElements.end(); itr++)
  {
    out_cellData.elementParents[itr->first] = tempMesh->elementParentIDs[eID];
  }

  // Check that each element is dominated by a material that is actually present in the original parent cell
  for (auto itr = newElements.begin(); itr != newElements.end(); itr++)
  {
    int parentElementID = out_cellData.elementParents[itr->first];
    int currentDominantMaterial = out_cellData.elementDominantMaterials[itr->first];

    if (originalMesh->materialVolumeFractionsElement[currentDominantMaterial][parentElementID] == 0)
    {
      // This material is not present in the original element from which the current element comes from
      if (currentDominantMaterial == matOneID)
        out_cellData.elementDominantMaterials[itr->first] = matTwoID;
      else
        out_cellData.elementDominantMaterials[itr->first] = matOneID;
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

  for (unsigned int i = 0; i < out_cellData.evInds.size(); ++i)
    out_cellData.evInds[i] -= evIndexSubtract[ out_cellData.evInds[i] ];

}

//--------------------------------------------------------------------------------

/// Finds the bit map representing the clipping case for a triangle.
unsigned int InterfaceReconstructor::determineTriangleClippingCase(mir::MIRMesh* tempMesh, const int matOneID, const int matTwoID, const int upperVertex, const int lowerLeftVertex, const int lowerRightVertex)
{
 // Determine the dominant color at each vertex
    int upperColor = matOneID, lowerLeftColor = matOneID, lowerRightColor = matOneID;

    if (matOneID == NULL_MAT)
      upperColor = lowerLeftColor = lowerRightColor = matTwoID;
    if (matTwoID == NULL_MAT)
      upperColor = lowerLeftColor = lowerRightColor = matOneID;
    if (matOneID != NULL_MAT && matTwoID != NULL_MAT)
    {
      upperColor = tempMesh->materialVolumeFractionsVertex[matOneID][upperVertex] > tempMesh->materialVolumeFractionsVertex[matTwoID][upperVertex] ? matOneID : matTwoID;
      lowerLeftColor = tempMesh->materialVolumeFractionsVertex[matOneID][lowerLeftVertex] > tempMesh->materialVolumeFractionsVertex[matTwoID][lowerLeftVertex] ? matOneID : matTwoID;
      lowerRightColor = tempMesh->materialVolumeFractionsVertex[matOneID][lowerRightVertex] > tempMesh->materialVolumeFractionsVertex[matTwoID][lowerRightVertex] ? matOneID : matTwoID;
    }
    
    // Create the index into the quad clipping lookup table using the dominant colors at each vertex
    unsigned int caseIndex = 0;
    if (upperColor == matOneID)  caseIndex |= 4;
    if (lowerLeftColor == matOneID)  caseIndex |= 2;
    if (lowerRightColor == matOneID) caseIndex |= 1;

    return caseIndex;
}

//--------------------------------------------------------------------------------

/// Generate the vertex position data for the new elements resulting from spltting a quad.
void InterfaceReconstructor::generateVertexPositionsFromTriangle(std::map<int, std::vector<int> > newVertices, mir::MIRMesh* tempMesh, axom::float64* verticesClippingTValue, int upperVertex, int lowerLeftVertex, int lowerRightVertex, CellData& out_cellData)
{
  std::map<int, mir::Point2> newPoints;           // hashmap of the new vertices' positions | Note: maps are ordered sets

  // Find the position coordinates of the vertices and then store the positions of the vertices in the return vector
  for (auto itr = newVertices.begin(); itr != newVertices.end(); itr++)
  {
    int vID = itr->first;
    if (vID == 0)
      newPoints[vID] = tempMesh->vertexPositions[upperVertex];
    if (vID == 1)
      newPoints[vID] = tempMesh->vertexPositions[lowerLeftVertex];
    if (vID == 2)
      newPoints[vID] = tempMesh->vertexPositions[lowerRightVertex];
    if (vID == 3)
      newPoints[vID] = interpolateVertexPosition(tempMesh->vertexPositions[upperVertex], tempMesh->vertexPositions[lowerLeftVertex], verticesClippingTValue[vID]);
    if (vID == 4)
      newPoints[vID] = interpolateVertexPosition(tempMesh->vertexPositions[lowerLeftVertex], tempMesh->vertexPositions[lowerRightVertex], verticesClippingTValue[vID]);
    if (vID == 5)
      newPoints[vID] = interpolateVertexPosition(tempMesh->vertexPositions[lowerRightVertex], tempMesh->vertexPositions[upperVertex], verticesClippingTValue[vID]);
  } 

  // Store the positions of the vertices in the return vector
  for (auto itr = newPoints.begin(); itr != newPoints.end(); itr++)
  {
    out_cellData.vertexPositions.push_back(itr->second);
  }
}

//--------------------------------------------------------------------------------

/// Generate the vertex volume fraction data for the new vertices resulting from splitting a quad. 
void InterfaceReconstructor::generateVertexVolumeFractionsFromTriangle(std::map<int, std::vector<int> > newVertices, mir::MIRMesh* tempMesh, axom::float64* verticesClippingTValue, int upperVertex, int lowerLeftVertex, int lowerRightVertex, CellData& out_cellData)
{
    // Calculate the vertex fractions at each vertex (use t value!)
    // Make sure the output volume fractions containers are the proper size
    out_cellData.vertexVolumeFractions.resize(tempMesh->numMaterials);
    
    for (auto itr = newVertices.begin(); itr != newVertices.end(); itr++)
    {
      int vID = itr->first;
      for (int matID = 0; matID < tempMesh->numMaterials; ++matID)
      {
        if (vID == 0)
          out_cellData.vertexVolumeFractions[matID].push_back(tempMesh->materialVolumeFractionsVertex[matID][upperVertex]);
        if (vID == 1)
          out_cellData.vertexVolumeFractions[matID].push_back(tempMesh->materialVolumeFractionsVertex[matID][lowerLeftVertex]);
        if (vID == 2)
          out_cellData.vertexVolumeFractions[matID].push_back(tempMesh->materialVolumeFractionsVertex[matID][lowerRightVertex]);
        if (vID == 3)
          out_cellData.vertexVolumeFractions[matID].push_back(lerpFloat(tempMesh->materialVolumeFractionsVertex[matID][upperVertex], tempMesh->materialVolumeFractionsVertex[matID][lowerLeftVertex], verticesClippingTValue[vID]));
        if (vID == 4)
          out_cellData.vertexVolumeFractions[matID].push_back(lerpFloat(tempMesh->materialVolumeFractionsVertex[matID][lowerLeftVertex], tempMesh->materialVolumeFractionsVertex[matID][lowerRightVertex], verticesClippingTValue[vID]));
        if (vID == 5)
          out_cellData.vertexVolumeFractions[matID].push_back(lerpFloat(tempMesh->materialVolumeFractionsVertex[matID][lowerRightVertex], tempMesh->materialVolumeFractionsVertex[matID][upperVertex], verticesClippingTValue[vID]));    
      }
    }
}

//--------------------------------------------------------------------------------

}
}
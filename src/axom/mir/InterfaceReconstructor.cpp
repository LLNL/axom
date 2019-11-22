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

void InterfaceReconstructor::computeReconstructedInterface(mir::MIRMesh& inputMesh, 
                                                           mir::MIRMesh& outputMesh)
{
  // Store a reference to the original mesh
  m_originalMesh = inputMesh;

  // Initialize the final mesh to be the same as the input mesh
  mir::MIRMesh finalMesh(m_originalMesh);

  // For each material in the mesh, split the mesh on the current material and generate a new mesh (input into next iteration)
  for (int matID = 0; matID < m_originalMesh.m_numMaterials; ++matID)
  {
    // Copy the mesh to be split
    mir::MIRMesh intermediateMesh(finalMesh);

    // Create an array to store the output of each element being split.
    CellData temp_cellData[intermediateMesh.m_elems.size()];

    // Process/split each element
    for (int eID = 0; eID < intermediateMesh.m_elems.size(); ++eID)
    {
      // Update the materials upon which the split will occur
      int currentDominantMat = intermediateMesh.m_elementDominantMaterials[eID];
      int matOne = matID;

      generateCleanCells(eID, currentDominantMat, matOne, intermediateMesh, temp_cellData[eID]);
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
  outputMesh = finalMesh;
  // return finalMesh;
}

//--------------------------------------------------------------------------------

void InterfaceReconstructor::computeReconstructedInterfaceIterative(mir::MIRMesh& inputMesh, 
                                                                    const int numIterations, 
                                                                    const axom::float64 percent,
                                                                    mir::MIRMesh& outputMesh)
{
  int numElems = inputMesh.m_elems.size();
  int numMaterials = inputMesh.m_numMaterials;

  // Make a copy of the original input mesh
  mir::MIRMesh meshToImprove(inputMesh);

  // Calculate the reconstruction on the unmodified, input mesh
  computeReconstructedInterface(meshToImprove, outputMesh);

  for (int it = 0; it < numIterations; ++it)
  {
    printf("iteration %d\n", it);
    // Calculate the output element volume fractions of the resulting output mesh
    std::vector<std::vector<axom::float64> > resultingElementVF = outputMesh.computeOriginalElementVolumeFractions();

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
    computeReconstructedInterface(meshToImprove, outputMesh);
  }
}

//--------------------------------------------------------------------------------

void  InterfaceReconstructor::generateCleanCells(const int eID, 
                                                 const int matOne, 
                                                 const int matTwo, 
                                                 mir::MIRMesh& tempMesh, 
                                                 CellData& out_cellData)
{
  mir::Shape shapeType = (mir::Shape) tempMesh.m_shapeTypes[eID];
  auto elementVertices = tempMesh.m_bdry[eID];

  // Set up data structures needed to compute how element should be clipped
  std::map<int, std::vector<int> > newElements;             // hashmap of the new elements' vertices | Note: maps are ordered sets
  std::map<int, std::vector<int> > newVertices;             // hashmap of the new vertices' elements | Note: maps are ordered sets
  std::vector<int> verticesPresent(mir::utilities::maxPossibleNumVerts(shapeType), 0); // Vector of flags denoting whether the vertex is present in the current case or not
  axom::float64* tValues = new axom::float64[ mir::utilities::maxPossibleNumVerts(shapeType) ]{0}; // Array of t values that denote the percent value of where the edge should be clipped
  

  // Set up the volume fractions for the current element for the two materials currently being considered
  std::vector<std::vector<axom::float64> > vertexVF(2);    
  for (int vID = 0; vID < elementVertices.size(); ++vID)
  {
    int originalVID = elementVertices[vID];
    if (matOne == NULL_MAT)
    {
      vertexVF[0].push_back(-1.0);
    }
    else
    {
      vertexVF[0].push_back( tempMesh.m_materialVolumeFractionsVertex[matOne][originalVID] );
    }
    
    if (matTwo == NULL_MAT)
    {
      vertexVF[1].push_back(-1.0);
    }
    else
    {
      vertexVF[1].push_back( tempMesh.m_materialVolumeFractionsVertex[matTwo][originalVID] );
    }
  }

  // Clip the element
  CellClipper clipper;
  clipper.computeClippingPoints(shapeType, vertexVF, newElements, newVertices, tValues);

  // Determine which vertices are present
  for (auto itr = newVertices.begin(); itr != newVertices.end(); itr++)
  {
    int vID = itr->first;
    verticesPresent[vID] = 1;
  }

  // Calculate the total number of elements and vertices that were generated from splitting the current element
  out_cellData.m_numElems = (int) newElements.size();
  out_cellData.m_numVerts = (int) newVertices.size();

  // Generate the topology and connectivity of the newly split element
  CellGenerator cellGenerator;
  cellGenerator.generateTopologyData(newElements, newVertices, out_cellData);

  // Generate the vertex position values of the newly split element
  std::vector<mir::Point2> originalElementVertexPositions;
  for (int vID = 0; vID < elementVertices.size(); ++vID)
  {
    int originalVID = elementVertices[vID];
    originalElementVertexPositions.push_back( tempMesh.m_vertexPositions[originalVID] );
  }
  cellGenerator.generateVertexPositions( shapeType, newVertices, originalElementVertexPositions, tValues, out_cellData);

  // Generate the vertex volume fractions of the newly split elements
  std::vector<std::vector<axom::float64> > originalElementVertexVF;
  for (int matID = 0; matID < tempMesh.m_numMaterials; ++matID)
  {
    std::vector<axom::float64> materialVertexVF;
    for (int vID = 0; vID < elementVertices.size(); ++vID)
    {
      int originalVID = elementVertices[vID];

      materialVertexVF.push_back( tempMesh.m_materialVolumeFractionsVertex[matID][originalVID] );
    }
    originalElementVertexVF.push_back( materialVertexVF );
  }

  cellGenerator.generateVertexVolumeFractions( shapeType, newVertices, originalElementVertexVF, tValues, out_cellData);

  // Determine and store the dominant material of this element
  for (auto itr = newElements.begin(); itr != newElements.end(); itr++)
  {
    int dominantMaterial = cellGenerator.determineCleanCellMaterial(shapeType, itr->second, matOne, matTwo, out_cellData.m_mapData.m_vertexVolumeFractions);
    out_cellData.m_mapData.m_elementDominantMaterials.push_back(dominantMaterial);
  }

  // Determine and store the parent of this element
  out_cellData.m_mapData.m_elementParents.resize(newElements.size());
  for (auto itr = newElements.begin(); itr != newElements.end(); itr++)
  {
    out_cellData.m_mapData.m_elementParents[itr->first] = tempMesh.m_elementParentIDs[eID];
  }

  // Determine the generated elements' shape types
  out_cellData.m_mapData.m_shapeTypes.resize(newElements.size());
  for (auto itr = newElements.begin(); itr != newElements.end(); itr++)
  {
    out_cellData.m_mapData.m_shapeTypes[itr->first] = cellGenerator.determineElementShapeType(shapeType, itr->second.size());
  }

  // TODO: Make a function that does this
  // Check that each element is dominated by a material that is actually present in the original parent cell
  for (auto itr = newElements.begin(); itr != newElements.end(); itr++)
  {
    int parentElementID = out_cellData.m_mapData.m_elementParents[itr->first];
    int currentDominantMaterial = out_cellData.m_mapData.m_elementDominantMaterials[itr->first];

    if (m_originalMesh.m_materialVolumeFractionsElement[currentDominantMaterial][parentElementID] == 0)
    {
      // This material is not present in the original element from which the current element comes from
      if (currentDominantMaterial == matOne)
        out_cellData.m_mapData.m_elementDominantMaterials[itr->first] = matTwo;
      else
        out_cellData.m_mapData.m_elementDominantMaterials[itr->first] = matOne;
    }
  }
  
  // Modify the evIndex values to account for the fact that perhaps not all 8 possible vertices are present
  std::vector<int> evIndexSubtract;
  evIndexSubtract.resize(mir::utilities::maxPossibleNumVerts(shapeType), 0);
  for (unsigned long i = 0; i < evIndexSubtract.size(); ++i)
  {
    if (verticesPresent[i] == 1)
      evIndexSubtract[i] = 0;
    else
      evIndexSubtract[i] = 1;
  }
  for (unsigned long i = 1; i < evIndexSubtract.size(); ++i)
    evIndexSubtract[i] += evIndexSubtract[i - 1];

  for (unsigned long i = 0; i < out_cellData.m_topology.m_evInds.size(); ++i)
  {
    out_cellData.m_topology.m_evInds[i] -= evIndexSubtract[ out_cellData.m_topology.m_evInds[i] ];
  }

  // Memory management
  delete[] tValues;
}

//--------------------------------------------------------------------------------

}
}

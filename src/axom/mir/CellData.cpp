// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "CellData.hpp"

namespace axom
{
namespace mir
{

  //--------------------------------------------------------------------------------

  CellData::CellData()
  {

  }

  //--------------------------------------------------------------------------------

  CellData::~CellData()
  {

  }

  //--------------------------------------------------------------------------------

  /// Merges the cell data from the given cell into this cell
  void CellData::mergeCell(CellData cellToMerge)
  {
    // Initialize index offsets
    int evBeginsOffset = topology.evInds.size();
    int veBeginsOffset = topology.veInds.size();

    int vertexIndexOffset = numVerts;
    int elementIndexOffset = numElems;

    // Merge the cell topology information
    for (unsigned long i = 0; i < cellToMerge.topology.evInds.size(); ++i)
    {
      topology.evInds.push_back(cellToMerge.topology.evInds[i] + vertexIndexOffset);
    }

    for (unsigned long i = 1; i < cellToMerge.topology.evBegins.size(); ++i)
    {
      topology.evBegins.push_back(cellToMerge.topology.evBegins[i] + evBeginsOffset);
    }

    for (unsigned long i = 0; i < cellToMerge.topology.veInds.size(); ++i)
    {
      topology.veInds.push_back(cellToMerge.topology.veInds[i] + elementIndexOffset);
    }

    for (unsigned long i = 1; i < cellToMerge.topology.veBegins.size(); ++i)
    {
      topology.veBegins.push_back(cellToMerge.topology.veBegins[i] + veBeginsOffset);
    }

    // Merge the vertex positions
    for (unsigned long i = 0; i < cellToMerge.mapData.vertexPositions.size(); ++i)
    {
      mapData.vertexPositions.push_back(cellToMerge.mapData.vertexPositions[i]);
    }

    // Merge the vertex volume fractions
    for (unsigned long matID = 0; matID < mapData.vertexVolumeFractions.size(); ++matID)
    {
      for (unsigned long vID = 0; vID < cellToMerge.mapData.vertexVolumeFractions[matID].size(); ++vID)
      {
        mapData.vertexVolumeFractions[matID].push_back(cellToMerge.mapData.vertexVolumeFractions[matID][vID]);
      }
    }

    // Merge the elements' dominant materials
    for (unsigned long i = 0; i < cellToMerge.mapData.elementDominantMaterials.size(); ++i)
    {
      mapData.elementDominantMaterials.push_back(cellToMerge.mapData.elementDominantMaterials[i]);
    }

    // Merge the elements' parent ids
    for (unsigned long i = 0; i < cellToMerge.mapData.elementParents.size(); ++i)
    {
      mapData.elementParents.push_back(cellToMerge.mapData.elementParents[i]);
    }
    
    // Merge the total number of verts and elems in the resulting cell
    numVerts += cellToMerge.numVerts;
    numElems += cellToMerge.numElems;
  }

  //--------------------------------------------------------------------------------

}
}
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

  CellData::CellData(int _numVerts, int _numElems, std::vector<PosType> _evInds, std::vector<PosType> _evBegins, std::vector<PosType> _veInds, std::vector<PosType> _veBegins, 
                      std::vector<mir::Point2> _vertexPositions, std::vector<std::vector<axom::float64> > _vertexVolumeFractions)
  {
    numVerts = _numVerts;
    numElems = _numElems;
    evInds = _evInds;
    evBegins = _evBegins;
    veInds = _veInds;
    veBegins = _veBegins;
    vertexPositions = _vertexPositions;
    vertexVolumeFractions = _vertexVolumeFractions;
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
    int evBeginsOffset = evInds.size();
    int veBeginsOffset = veInds.size();

    int vertexIndexOffset = numVerts;
    int elementIndexOffset = numElems;

    // Merge the cell topology information
    for (auto i = 0; i < cellToMerge.evInds.size(); ++i)
    {
      evInds.push_back(cellToMerge.evInds[i] + vertexIndexOffset);
    }

    for (auto i = 1; i < cellToMerge.evBegins.size(); ++i)
    {
      evBegins.push_back(cellToMerge.evBegins[i] + evBeginsOffset);
    }

    for (auto i = 0; i < cellToMerge.veInds.size(); ++i)
    {
      veInds.push_back(cellToMerge.veInds[i] + elementIndexOffset);
    }

    for (auto i = 1; i < cellToMerge.veBegins.size(); ++i)
    {
      veBegins.push_back(cellToMerge.veBegins[i] + veBeginsOffset);
    }

    // Merge the vertex positions
    for (auto i = 0; i < cellToMerge.vertexPositions.size(); ++i)
    {
      vertexPositions.push_back(cellToMerge.vertexPositions[i]);
    }

    // Merge the vertex volume fractions
    for (auto matID = 0; matID < vertexVolumeFractions.size(); ++matID)
    {
      for (auto vID = 0; vID < cellToMerge.vertexVolumeFractions[matID].size(); ++vID)
      {
        vertexVolumeFractions[matID].push_back(cellToMerge.vertexVolumeFractions[matID][vID]);
      }
    }
    
    // Merge the total number of verts and elems in the resulting cell
    numVerts += cellToMerge.numVerts;
    numElems += cellToMerge.numElems;

  }

  //--------------------------------------------------------------------------------

}
}
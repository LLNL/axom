// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef __CELL_DATA_H__
#define __CELL_DATA_H__

#include "axom/core.hpp"  // for axom macros
#include "axom/slam.hpp"

#include "MIRMesh.hpp"

namespace numerics = axom::numerics;
namespace slam = axom::slam;

namespace axom
{
namespace mir
{

  class CellData
  {

    public:
      CellData();
      CellData(int _numVerts, int _numElems, std::vector<PosType> _evInds, std::vector<PosType> _evBegins, std::vector<PosType> _veInds, std::vector<PosType> _veBegins, 
                      std::vector<mir::Point2> _vertexPositions, std::vector<std::vector<axom::float64> > _vertexVolumeFractions);
      ~CellData();

      void mergeCell(CellData cellToMerge);

    public:
      int numVerts;
      int numElems;

      // Cell connectivity/topology
      std::vector<PosType> evInds;
      std::vector<PosType> evBegins;
      std::vector<PosType> veInds;
      std::vector<PosType> veBegins;

      std::vector<mir::Point2> vertexPositions;
      std::vector<axom::float64> materialsInCell;     // TODO: This is not currently being used.
      std::vector<std::vector<axom::float64> > vertexVolumeFractions;
  };

}
}

#endif
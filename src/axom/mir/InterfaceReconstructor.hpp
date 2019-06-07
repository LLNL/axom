// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef __INTERFACE_RECONSTRUCTOR_H__
#define __INTERFACE_RECONSTRUCTOR_H__

#include "axom/core.hpp"  // for axom macros
#include "axom/slam.hpp"

#include "MIRMesh.hpp"
#include "CellData.hpp"
#include "ZooBitMaps.hpp"
#include <map>

namespace numerics = axom::numerics;
namespace slam = axom::slam;

namespace axom
{
namespace mir
{
  class InterfaceReconstructor
  {
    public:
      InterfaceReconstructor();
      InterfaceReconstructor(mir::MIRMesh* _mesh);
      ~InterfaceReconstructor();

      mir::MIRMesh              computeReconstructedInterface();

      
      // std::vector<ScalarMap>    computeVolumeFractionAverages(mir::MIRMesh* tempMesh);

    private:
      mir::MIRMesh* mesh;
    
    public:
      std::vector<ScalarMap> materialVolumeFractionsVertex; // volume fractions for each material for each vertex

    private:
      void                computeClippingPoints(const int eID, const int matOneID, const int matTwoID, mir::MIRMesh* tempMesh,
                                                std::vector<mir::PosType>& _evInds, std::vector<mir::PosType>& _evBegins, std::vector<mir::PosType>& _veInds, std::vector<mir::PosType>& _veBegins,
                                                std::vector<mir::Point2>& _vertexPositions, std::vector<axom::float64*>& _materialInCell, int* _numVerts, int* _numElements, std::vector<std::vector<axom::float64> >& _newVolumeFractionsAtVerts);
                                                
      void                computeTriangleClippingPoints(const int eID, const int matOneID, const int matTwoID, mir::MIRMesh* tempMesh,
                                                std::vector<mir::PosType>& _evInds, std::vector<mir::PosType>& _evBegins, std::vector<mir::PosType>& _veInds, std::vector<mir::PosType>& _veBegins,
                                                std::vector<mir::Point2>& _vertexPositions, std::vector<axom::float64*>& _materialInCell, int* _numVerts, int* _numElements, std::vector<std::vector<axom::float64> >& _newVolumeFractionsAtVerts);

      void                computeQuadClippingPoints(const int eID, const int matOneID, const int matTwoID, mir::MIRMesh* tempMesh,
                                                std::vector<mir::PosType>& out_evInds, std::vector<mir::PosType>& out_evBegins, std::vector<mir::PosType>& out_veInds, std::vector<mir::PosType>& out_veBegins,
                                                std::vector<mir::Point2>& out_vertexPositions, std::vector<axom::float64*>& out_materialInCell, int* out_numVerts, int* out_numElements, std::vector<std::vector<axom::float64> >& out_newVolumeFractionsAtVerts);

      mir::Point2         interpolateVertexPosition(mir::Point2 vertexOnePos, mir::Point2 vertexTwoPos, float t);
      axom::float64       lerpFloat(axom::float64 f0, axom::float64 f1, axom::float64 t);
      axom::float64       computeClippingPointOnEdge(const int vertexOneID, const int vertexTwoID, const int matOneID, const int matTwoID, mir::MIRMesh* tempMesh);

      void                mergeCellSplitData();

  };
}
}
#endif
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
#include "ZooClippingTables.hpp"
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
      ~InterfaceReconstructor();

      mir::MIRMesh              computeReconstructedInterface(mir::MIRMesh* inputMesh);
      mir::MIRMesh              computeReconstructedInterfaceIterative(mir::MIRMesh* inputMesh, int numIterations, axom::float64 percent);

    private:
      mir::MIRMesh* originalMesh;

    private:

      // general clipping function
      void                computeQuadClippingPoints(const int eID, const int matOneID, const int matTwoID, mir::MIRMesh* tempMesh, CellData& out_cellData);

      // triangle clipping function
      void                computeTriangleClippingPoints(const int eID, const int matOneID, const int matTwoID, mir::MIRMesh* tempMesh, CellData& out_cellData);

      // quad clipping functions
      void                computeClippingPoints(const int eID, const int matOneID, const int matTwoID, mir::MIRMesh* tempMesh, CellData& out_cellData);

      // clipping interpolation functions
      mir::Point2         interpolateVertexPosition(mir::Point2 vertexOnePos, mir::Point2 vertexTwoPos, float t);
      axom::float64       lerpFloat(axom::float64 f0, axom::float64 f1, axom::float64 t);
      axom::float64       computeClippingPointOnEdge(const int vertexOneID, const int vertexTwoID, const int matOneID, const int matTwoID, mir::MIRMesh* tempMesh);

      // general helper functions
      void                generateTopologyData(std::map<int, std::vector<int> > newElements, std::map<int, std::vector<int> > newVertices, CellData& out_cellData);

      // quad clipping points helper functions
      unsigned int        determineQuadClippingCase(mir::MIRMesh* tempMesh, const int matOneID, const int matTwoID, const int upperLeftVertex, const int lowerLeftVertex, const int lowerRightVertex, const int upperRightVertex);
      void                generateVertexPositionsFromQuad(std::map<int, std::vector<int> > newVertices, mir::MIRMesh* tempMesh, axom::float64* verticesClippingTValue, int upperLeftVertex, int lowerLeftVertex, int lowerRightVertex, int upperRightVertex, CellData& out_cellData);
      void                generateVertexVolumeFractionsFromQuad(std::map<int, std::vector<int> > newVertices, mir::MIRMesh* tempMesh, axom::float64* verticesClippingTValue, int upperLeftVertex, int lowerLeftVertex, int lowerRightVertex, int upperRightVertex, CellData& out_cellData);

      // triangle clipping points helper functions
      unsigned int        determineTriangleClippingCase(mir::MIRMesh* tempMesh, const int matOneID, const int matTwoID, const int upperVertex, const int lowerLeftVertex, const int lowerRightVertex);
      void                generateVertexPositionsFromTriangle(std::map<int, std::vector<int> > newVertices, mir::MIRMesh* tempMesh, axom::float64* verticesClippingTValue, int upperVertex, int lowerLeftVertex, int lowerRightVertex, CellData& out_cellData);
      void                generateVertexVolumeFractionsFromTriangle(std::map<int, std::vector<int> > newVertices, mir::MIRMesh* tempMesh, axom::float64* verticesClippingTValue, int upperVertex, int lowerLeftVertex, int lowerRightVertex, CellData& out_cellData);
  };
}
}
#endif
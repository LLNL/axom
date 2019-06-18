// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef __MESH_TESTER_H__
#define __MESH_TESTER_H__

#include "axom/core.hpp"  // for axom macros
#include "axom/slam.hpp"  // unified header for slam classes and functions

#include "MIRMesh.hpp"

#include <algorithm>

namespace numerics = axom::numerics;
namespace slam = axom::slam;

namespace axom
{
namespace mir
{
  class MeshTester
  {
    public:
      MeshTester();
      ~MeshTester();

    public:
      MIRMesh initTestCaseOne(); // 3x3 grid, 2 materials
      mir::MIRMesh initTestCaseTwo(); // 3x3 grid, 3 materials
      mir::MIRMesh initTestCaseThree(); // triforce, 2 materials
      mir::MIRMesh initTestCaseFour(); // 3x3, 2 materials, circle of material
      mir::MIRMesh createUniformGridTestCaseMesh(int gridSize, mir::Point2 circleCenter, axom::float64 circleRadius);
      mir::MIRMesh initTestCaseFive(int gridSize, int numCircles); // multiple materials, multiple concentric circles
            
    private:
      axom::float64 distance(mir::Point2 p0, mir::Point2 p1);
      mir::CellData generateGrid(int n);
      axom::float64 calculatePercentOverlapMonteCarlo(int n, mir::Point2 circle_center, axom::float64 circle_radius, mir::Point2 q_p0, mir::Point2 q_p1, mir::Point2 q_p2, mir::Point2 q_p3);
      
      int circleQuadCornersOverlaps(mir::Point2 circle_center, axom::float64 circle_radius, mir::Point2 q_p0, mir::Point2 q_p1, mir::Point2 q_p2, mir::Point2 q_p3);
  };
}
}

#endif

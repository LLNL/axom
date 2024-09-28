// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "MeshTester.hpp"
#include <axom/mir.hpp>

namespace numerics = axom::numerics;
namespace slam = axom::slam;
namespace bputils = axom::mir::utilities::blueprint;

namespace axom
{
namespace mir
{
//--------------------------------------------------------------------------------
/**
       * \brief Calculates the percent overlap between the given circle and quad.
       * 
       * \param gridSize  The size of the uniform grid which will be sampled over to check for overlap.
       * \param circleCenter  The center point of the circle.
       * \param circleRadius  The radius of the circle.
       * \param quadP0  The upper left vertex of the quad.
       * \param quadP1  The lower left vertex of the quad.
       * \param quadP2  The lower right vertex of the quad.
       * \param quadP3  The upper right vertex of the quad.
       * 
       * /return The percent value overlap of the circle and the quad between [0, 1].
       */

template <typename PointType>
static axom::float64 calculatePercentOverlapMonteCarlo(int gridSize,
                                                       const PointType& circleCenter,
                                                       axom::float64 circleRadius,
                                                       const PointType& quadP0,
                                                       const PointType& quadP1,
                                                       const PointType& quadP2,
                                                       const PointType& quadP3)
{
  // Check if any of the quad's corners are within the circle
  auto d0Sq = primal::squared_distance(quadP0, circleCenter);
  auto d1Sq = primal::squared_distance(quadP1, circleCenter);
  auto d2Sq = primal::squared_distance(quadP2, circleCenter);
  auto d3Sq = primal::squared_distance(quadP3, circleCenter);
  auto dRSq = circleRadius * circleRadius;

  int inFlags = ((d0Sq < dRSq) ? 1 << 0 : 0) + ((d1Sq < dRSq) ? 1 << 1 : 0) +
    ((d2Sq < dRSq) ? 1 << 2 : 0) + ((d3Sq < dRSq) ? 1 << 3 : 0);
  const int allFlags = 15;
  const int noFlags = 0;

  if(inFlags == allFlags)
  {
    // The entire quad overlaps the circle
    return 1.;
  }
  else if(inFlags == noFlags)
  {
    return 0.;
  }
  else
  {
    // Some of the quad overlaps the circle, so run the Monte Carlo sampling to determine how much
    axom::float64 delta_x = axom::utilities::abs(quadP2[0] - quadP1[0]) /
      static_cast<double>(gridSize - 1);
    axom::float64 delta_y = axom::utilities::abs(quadP0[1] - quadP1[1]) /
      static_cast<double>(gridSize - 1);
    int countOverlap = 0;
    for(int y = 0; y < gridSize; ++y)
    {
      for(int x = 0; x < gridSize; ++x)
      {
        PointType samplePoint = PointType::make_point(delta_x * x + quadP1[0],
                                                      delta_y * y + quadP1[1]);
        if(primal::squared_distance(samplePoint, circleCenter) < dRSq)
          ++countOverlap;
      }
    }
    return countOverlap / static_cast<double>(gridSize * gridSize);
  }
}

//--------------------------------------------------------------------------------

MIRMesh MeshTester::initTestCaseOne()
{
  mir::CellTopologyData topoData;
  mir::CellMapData mapData;
  mir::CellData cellData;
  VolumeFractions volFracs;

  int numElements = 9;
  int numVertices = 16;
  mir::VertSet verts(numVertices);  // Construct a vertex set with 16 vertices
  mir::ElemSet elems(numElements);  // Construct an element set with 9 elements

  // Create the mesh connectivity information
  topoData.m_evInds = {
    0,  4,  5,  1,   // elem 0, card 4, start 0
    1,  5,  6,  2,   // elem 1, card 4, start 4
    2,  6,  7,  3,   // elem 2, card 4, start 8
    4,  8,  9,  5,   // elem 3, card 4, start 12
    5,  9,  10, 6,   // elem 4, card 4, start 16
    6,  10, 11, 7,   // elem 5, card 4, start 20
    8,  12, 13, 9,   // elem 6, card 4, start 24
    9,  13, 14, 10,  // elem 7, card 4, start 28
    10, 14, 15, 11   // elem 8, card 4, start 32, end 36
  };

  topoData.m_evBegins = {0, 4, 8, 12, 16, 20, 24, 28, 32, 36};
  topoData.m_veInds = {
    0,           // vert  0, card 1, start 0
    0, 1,        // vert  1, card 2, start 1
    1, 2,        // vert  2, card 2, start 3
    2,           // vert  3, card 1, start 5
    0, 3,        // vert  4, card 2, start 6
    0, 1, 3, 4,  // vert  5, card 4, start 8
    1, 2, 4, 5,  // vert  6, card 4, start 12
    2, 5,        // vert  7, card 2, start 16
    3, 6,        // vert  8, card 2, start 18
    3, 4, 6, 7,  // vert  9, card 4, start 20
    4, 5, 7, 8,  // vert  10, card 4, start 24
    5, 8,        // vert  11, card 2, start 28
    6,           // vert  12, card 1, start 30
    6, 7,        // vert  13, card 2, start 31
    7, 8,        // vert  14, card 2, start 33
    8,           // vert  15, card 1, start 35, end 36
  };
  topoData
    .m_veBegins = {0, 1, 3, 5, 6, 8, 12, 16, 18, 20, 24, 28, 30, 31, 33, 35, 36};

  int numMaterials = 2;
  enum
  {
    GREEN = 0,
    BLUE = 1
  };

  volFracs.resize(numMaterials);

  volFracs[GREEN] = {1.0, 1.0, 1.0, 1.0, 0.5, 0.2, 0.2, 0.0, 0.0};
  volFracs[BLUE] = {0.0, 0.0, 0.0, 0.0, 0.5, 0.8, 0.8, 1.0, 1.0};

  mapData.m_vertexPositions = {mir::Point2::make_point(0.0, 3.0),
                               mir::Point2::make_point(1.0, 3.0),
                               mir::Point2::make_point(2.0, 3.0),
                               mir::Point2::make_point(3.0, 3.0),

                               mir::Point2::make_point(0.0, 2.0),
                               mir::Point2::make_point(1.0, 2.0),
                               mir::Point2::make_point(2.0, 2.0),
                               mir::Point2::make_point(3.0, 2.0),

                               mir::Point2::make_point(0.0, 1.0),
                               mir::Point2::make_point(1.0, 1.0),
                               mir::Point2::make_point(2.0, 1.0),
                               mir::Point2::make_point(3.0, 1.0),

                               mir::Point2::make_point(0.0, 0.0),
                               mir::Point2::make_point(1.0, 0.0),
                               mir::Point2::make_point(2.0, 0.0),
                               mir::Point2::make_point(3.0, 0.0)};

  mapData.m_elementDominantMaterials = Vec<int>(numElements, NULL_MAT);
  mapData.m_elementParents =
    {0, 1, 2, 3, 4, 5, 6, 7, 8};  // For the base mesh, the parents are always themselves
  mapData.m_shapeTypes = Vec<mir::Shape>(numElements, mir::Shape::Quad);

  // Build the mesh
  mir::MIRMesh testMesh;
  testMesh.initializeMesh(verts, elems, numMaterials, topoData, mapData, volFracs);

  return testMesh;
}

//--------------------------------------------------------------------------------
void MeshTester::mesh3x3(conduit::Node& mesh)
{
  // clang-format off
  mesh["coordsets/coords/type"] = "explicit";
  mesh["coordsets/coords/values/x"].set(std::vector<float>{{
    0., 1., 2., 3., 0., 1., 2., 3., 0., 1., 2., 3., 0., 1., 2., 3.
  }});
  mesh["coordsets/coords/values/y"].set(std::vector<float>{{
    0., 0., 0., 0., 1., 1., 1., 1., 2., 2., 2., 2., 3., 3., 3., 3.
  }});

  mesh["topologies/mesh/type"] = "unstructured";
  mesh["topologies/mesh/coordset"] = "coords";
  mesh["topologies/mesh/elements/shape"] = "quad";
  mesh["topologies/mesh/elements/connectivity"].set(std::vector<int>{{
    0,1,5,4, 1,2,6,5, 2,3,7,6, 4,5,9,8, 5,6,10,9, 6,7,11,10, 8,9,13,12, 9,10,14,13, 10,11,15,14
  }});
  mesh["topologies/mesh/elements/sizes"].set(std::vector<int>{{
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4
  }});
  mesh["topologies/mesh/elements/offsets"].set(std::vector<int>{{
    0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60
  }});
  // clang-format on
}

//--------------------------------------------------------------------------------
void MeshTester::initTestCaseOne(conduit::Node& mesh)
{
  mesh3x3(mesh);

  // clang-format off
  mesh["matsets/mat/topology"] = "mesh";
  mesh["matsets/mat/material_map/green"] = 0;
  mesh["matsets/mat/material_map/blue"] = 1;
  mesh["matsets/mat/material_ids"].set(std::vector<int>{{
   0,
   0,
   0,
   0,
   0, 1,
   0, 1,
   0, 1,
   1,
   1
  }});
  mesh["matsets/mat/volume_fractions"].set(std::vector<float>{{
   1.,
   1.,
   1.,
   1.,
   0.5, 0.5,
   0.2, 0.8,
   0.2, 0.8,
   1.,
   1.
  }});
  mesh["matsets/mat/sizes"].set(std::vector<int>{{
    1, 1, 1, 1, 2, 2, 2, 1, 1
  }});
  mesh["matsets/mat/offsets"].set(std::vector<int>{{
    //0, 1, 2, 3, 5, 7, 9, 10, 11
    0, 1, 2, 3, 4, 6, 8, 10, 11
  }});
  mesh["matsets/mat/indices"].set(std::vector<int>{{
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11
  }});
  // clang-format on
}

//--------------------------------------------------------------------------------

mir::MIRMesh MeshTester::initTestCaseTwo()
{
  mir::CellTopologyData topoData;
  mir::CellMapData mapData;
  mir::CellData cellData;
  VolumeFractions volFracs;

  int numElements = 9;
  int numVertices = 16;
  mir::VertSet verts(numVertices);  // Construct a vertex set with 16 vertices
  mir::ElemSet elems(numElements);  // Construct an element set with 9 elements

  // Create the mesh connectivity information
  topoData.m_evInds = {
    0,  4,  5,  1,   // elem 0, card 4, start 0
    1,  5,  6,  2,   // elem 1, card 4, start 4
    2,  6,  7,  3,   // elem 2, card 4, start 8
    4,  8,  9,  5,   // elem 3, card 4, start 12
    5,  9,  10, 6,   // elem 4, card 4, start 16
    6,  10, 11, 7,   // elem 5, card 4, start 20
    8,  12, 13, 9,   // elem 6, card 4, start 24
    9,  13, 14, 10,  // elem 7, card 4, start 28
    10, 14, 15, 11   // elem 8, card 4, start 32, end 36
  };

  topoData.m_evBegins = {0, 4, 8, 12, 16, 20, 24, 28, 32, 36};
  topoData.m_veInds = {
    0,           // vert  0, card 1, start 0
    0, 1,        // vert  1, card 2, start 1
    1, 2,        // vert  2, card 2, start 3
    2,           // vert  3, card 1, start 5
    0, 3,        // vert  4, card 2, start 6
    0, 1, 3, 4,  // vert  5, card 4, start 8
    1, 2, 4, 5,  // vert  6, card 4, start 12
    2, 5,        // vert  7, card 2, start 16
    3, 6,        // vert  8, card 2, start 18
    3, 4, 6, 7,  // vert  9, card 4, start 20
    4, 5, 7, 8,  // vert  10, card 4, start 24
    5, 8,        // vert  11, card 2, start 28
    6,           // vert  12, card 1, start 30
    6, 7,        // vert  13, card 2, start 31
    7, 8,        // vert  14, card 2, start 33
    8,           // vert  15, card 1, start 35, end 36
  };
  topoData
    .m_veBegins = {0, 1, 3, 5, 6, 8, 12, 16, 18, 20, 24, 28, 30, 31, 33, 35, 36};

  int numMaterials = 3;
  enum
  {
    BLUE = 0,
    RED = 1,
    ORANGE = 2
  };

  volFracs.resize(numMaterials);
  volFracs[BLUE] = {1.0, 1.0, 1.0, 1.0, 0.5, 0.2, 0.2, 0.0, 0.0};
  volFracs[RED] = {0.0, 0.0, 0.0, 0.0, 0.3, 0.8, 0.0, 0.3, 1.0};
  volFracs[ORANGE] = {0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.8, 0.7, 0.0};

  mapData.m_vertexPositions = {mir::Point2::make_point(0.0, 3.0),
                               mir::Point2::make_point(1.0, 3.0),
                               mir::Point2::make_point(2.0, 3.0),
                               mir::Point2::make_point(3.0, 3.0),

                               mir::Point2::make_point(0.0, 2.0),
                               mir::Point2::make_point(1.0, 2.0),
                               mir::Point2::make_point(2.0, 2.0),
                               mir::Point2::make_point(3.0, 2.0),

                               mir::Point2::make_point(0.0, 1.0),
                               mir::Point2::make_point(1.0, 1.0),
                               mir::Point2::make_point(2.0, 1.0),
                               mir::Point2::make_point(3.0, 1.0),

                               mir::Point2::make_point(0.0, 0.0),
                               mir::Point2::make_point(1.0, 0.0),
                               mir::Point2::make_point(2.0, 0.0),
                               mir::Point2::make_point(3.0, 0.0)};

  mapData.m_elementDominantMaterials = Vec<int>(numElements, NULL_MAT);
  mapData.m_elementParents =
    {0, 1, 2, 3, 4, 5, 6, 7, 8};  // For the base mesh, the parents are always themselves
  mapData.m_shapeTypes = Vec<mir::Shape>(numElements, mir::Shape::Quad);

  // Build the mesh
  mir::MIRMesh testMesh;
  testMesh.initializeMesh(verts, elems, numMaterials, topoData, mapData, volFracs);

  return testMesh;
}

//--------------------------------------------------------------------------------
void MeshTester::initTestCaseTwo(conduit::Node& mesh)
{
  mesh3x3(mesh);

  // clang-format off
  constexpr int BLUE = 0;
  constexpr int RED = 1;
  constexpr int ORANGE = 2;
  mesh["matsets/mat/topology"] = "mesh";
  mesh["matsets/mat/material_map/blue"] = BLUE;
  mesh["matsets/mat/material_map/red"] = RED;
  mesh["matsets/mat/material_map/orange"] = ORANGE;
  mesh["matsets/mat/material_ids"].set(std::vector<int>{{
   BLUE,
   BLUE,
   BLUE,
   BLUE,
   BLUE, RED, ORANGE,
   BLUE, RED,
   BLUE, ORANGE,
   RED, ORANGE,
   RED
  }});
  mesh["matsets/mat/volume_fractions"].set(std::vector<float>{{
   1.,
   1.,
   1.,
   1.,
   0.5, 0.3, 0.2,
   0.2, 0.8,
   0.2, 0.8,
   0.3, 0.7,
   1.
  }});
  mesh["matsets/mat/sizes"].set(std::vector<int>{{
    1, 1, 1, 1, 3, 2, 2, 2, 1
  }});
  mesh["matsets/mat/offsets"].set(std::vector<int>{{
    0, 1, 2, 3, 4, 7, 9, 11, 13
  }});
  mesh["matsets/mat/indices"].set(std::vector<int>{{
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13
  }});
  // clang-format on
}

//--------------------------------------------------------------------------------

mir::MIRMesh MeshTester::initTestCaseThree()
{
  mir::CellTopologyData topoData;
  mir::CellMapData mapData;
  mir::CellData cellData;
  VolumeFractions volFracs;

  int numElements = 4;
  int numVertices =
    6;  // OR create a middle triangle with all of one material, and then a ring of triangles around it that are full of the other material

  mir::VertSet verts = mir::VertSet(numVertices);
  mir::ElemSet elems = mir::ElemSet(numElements);

  // Create the mesh connectivity information
  topoData.m_evInds = {
    0,
    1,
    2,  // elem 0, card 3, start 0
    1,
    3,
    4,  // elem 1, card 3, start 3
    1,
    4,
    2,  // elem 2, card 3, start 6
    2,
    4,
    5  // elem 3, card 3, start 9, end 12
  };

  topoData.m_evBegins = {0, 3, 6, 9, 12};
  topoData.m_veInds = {
    0,  // vert  0, card 1, start 0
    0,
    1,
    2,  // vert  1, card 3, start 1
    0,
    2,
    3,  // vert  2, card 3, start 4
    1,  // vert  3, card 1, start 7
    1,
    2,
    3,  // vert  4, card 3, start 8
    3   // vert  5, card 1, start 11, end 12
  };
  topoData.m_veBegins = {0, 1, 4, 7, 8, 11, 12};

  int numMaterials = 2;
  enum
  {
    BLUE = 0,
    RED = 1,
  };

  volFracs.resize(numMaterials);
  volFracs[BLUE] = {0.0, 0.5, 0.8, 0.5};
  volFracs[RED] = {1.0, 0.5, 0.2, 0.5};

  mapData.m_vertexPositions = {mir::Point2::make_point(1.0, 2.0),
                               mir::Point2::make_point(0.5, 1.0),
                               mir::Point2::make_point(1.5, 1.0),
                               mir::Point2::make_point(0.0, 0.0),
                               mir::Point2::make_point(1.0, 0.0),
                               mir::Point2::make_point(2.0, 0.0)};

  mapData.m_elementDominantMaterials = Vec<int>(numElements, NULL_MAT);
  mapData.m_elementParents = {0, 1, 2, 3};  // For the base mesh, the parents are always themselves
  mapData.m_shapeTypes = Vec<mir::Shape>(numElements, mir::Shape::Triangle);

  // Build the mesh
  mir::MIRMesh testMesh;
  testMesh.initializeMesh(verts, elems, numMaterials, topoData, mapData, volFracs);

  return testMesh;
}

void MeshTester::initTestCaseThree(conduit::Node& mesh)
{
  // clang-format off
  mesh["coordsets/coords/type"] = "explicit";
  mesh["coordsets/coords/values/x"].set(std::vector<float>{{
    1., 0.5, 1.5, 0., 1., 2.,
  }});
  mesh["coordsets/coords/values/y"].set(std::vector<float>{{
    2., 1., 1., 0., 0., 0.
  }});

  mesh["topologies/mesh/type"] = "unstructured";
  mesh["topologies/mesh/coordset"] = "coords";
  mesh["topologies/mesh/elements/shape"] = "tri";
  mesh["topologies/mesh/elements/connectivity"].set(std::vector<int>{{
    0,1,2, 1,3,4, 1,4,2, 2,4,5
  }});
  mesh["topologies/mesh/elements/sizes"].set(std::vector<int>{{3,3,3,3}});
  mesh["topologies/mesh/elements/offsets"].set(std::vector<int>{{0,3,6,9}});

  constexpr int BLUE = 0;
  constexpr int RED = 1;
  mesh["matsets/mat/topology"] = "mesh";
  mesh["matsets/mat/material_map/blue"] = BLUE;
  mesh["matsets/mat/material_map/red"] = RED;
  mesh["matsets/mat/material_ids"].set(std::vector<int>{{
   RED,
   BLUE, RED,
   BLUE, RED,
   BLUE, RED
  }});
  mesh["matsets/mat/volume_fractions"].set(std::vector<float>{{
   1.,
   0.5, 0.5,
   0.8, 0.2,
   0.5, 0.5
  }});
  mesh["matsets/mat/sizes"].set(std::vector<int>{{
    1, 2, 2, 2
  }});
  mesh["matsets/mat/offsets"].set(std::vector<int>{{
    0, 1, 3, 5
  }});
  mesh["matsets/mat/indices"].set(std::vector<int>{{
    0, 1, 2, 3, 4, 5, 6
  }});
  // clang-format on
}

//--------------------------------------------------------------------------------

mir::MIRMesh MeshTester::initTestCaseFour()
{
  mir::CellTopologyData topoData;
  mir::CellMapData mapData;
  mir::CellData cellData;
  VolumeFractions volFracs;

  int numElements = 9;
  int numVertices = 16;
  mir::VertSet verts =
    mir::VertSet(numVertices);  // Construct a vertex set with 16 vertices
  mir::ElemSet elems =
    mir::ElemSet(numElements);  // Construct an element set with 9 elements

  // Create the mesh connectivity information
  topoData.m_evInds = {
    0,  4,  5,  1,   // elem 0, card 4, start 0
    1,  5,  6,  2,   // elem 1, card 4, start 4
    2,  6,  7,  3,   // elem 2, card 4, start 8
    4,  8,  9,  5,   // elem 3, card 4, start 12
    5,  9,  10, 6,   // elem 4, card 4, start 16
    6,  10, 11, 7,   // elem 5, card 4, start 20
    8,  12, 13, 9,   // elem 6, card 4, start 24
    9,  13, 14, 10,  // elem 7, card 4, start 28
    10, 14, 15, 11   // elem 8, card 4, start 32, end 36
  };

  topoData.m_evBegins = {0, 4, 8, 12, 16, 20, 24, 28, 32, 36};
  topoData.m_veInds = {
    0,           // vert  0, card 1, start 0
    0, 1,        // vert  1, card 2, start 1
    1, 2,        // vert  2, card 2, start 3
    2,           // vert  3, card 1, start 5
    0, 3,        // vert  4, card 2, start 6
    0, 1, 3, 4,  // vert  5, card 4, start 8
    1, 2, 4, 5,  // vert  6, card 4, start 12
    2, 5,        // vert  7, card 2, start 16
    3, 6,        // vert  8, card 2, start 18
    3, 4, 6, 7,  // vert  9, card 4, start 20
    4, 5, 7, 8,  // vert  10, card 4, start 24
    5, 8,        // vert  11, card 2, start 28
    6,           // vert  12, card 1, start 30
    6, 7,        // vert  13, card 2, start 31
    7, 8,        // vert  14, card 2, start 33
    8,           // vert  15, card 1, start 35, end 36
  };
  topoData
    .m_veBegins = {0, 1, 3, 5, 6, 8, 12, 16, 18, 20, 24, 28, 30, 31, 33, 35, 36};

  mapData.m_vertexPositions = {mir::Point2::make_point(0.0, 3.0),
                               mir::Point2::make_point(1.0, 3.0),
                               mir::Point2::make_point(2.0, 3.0),
                               mir::Point2::make_point(3.0, 3.0),

                               mir::Point2::make_point(0.0, 2.0),
                               mir::Point2::make_point(1.0, 2.0),
                               mir::Point2::make_point(2.0, 2.0),
                               mir::Point2::make_point(3.0, 2.0),

                               mir::Point2::make_point(0.0, 1.0),
                               mir::Point2::make_point(1.0, 1.0),
                               mir::Point2::make_point(2.0, 1.0),
                               mir::Point2::make_point(3.0, 1.0),

                               mir::Point2::make_point(0.0, 0.0),
                               mir::Point2::make_point(1.0, 0.0),
                               mir::Point2::make_point(2.0, 0.0),
                               mir::Point2::make_point(3.0, 0.0)};

  int numMaterials = 2;
  enum
  {
    GREEN = 0,
    BLUE = 1
  };

  volFracs.resize(numMaterials);

  auto& greenVolumeFractions = volFracs[GREEN];
  auto& blueVolumeFractions = volFracs[BLUE];
  const auto& points = mapData.m_vertexPositions;
  const auto& evInds = topoData.m_evInds;

  greenVolumeFractions.resize(numElements);
  blueVolumeFractions.resize(numElements);

  // Generate the element volume fractions for the circle
  auto circleCenter = mir::Point2::make_point(1.5, 1.5);
  axom::float64 circleRadius = 1.25;
  int gridSize = 1000;
  for(int i = 0; i < numElements; ++i)
  {
    auto vf = calculatePercentOverlapMonteCarlo(gridSize,
                                                circleCenter,
                                                circleRadius,
                                                points[evInds[i * 4 + 0]],
                                                points[evInds[i * 4 + 1]],
                                                points[evInds[i * 4 + 2]],
                                                points[evInds[i * 4 + 3]]);
    greenVolumeFractions[i] = vf;
    blueVolumeFractions[i] = 1.0 - vf;
  }

  mapData.m_elementDominantMaterials = Vec<int>(numElements, NULL_MAT);
  mapData.m_elementParents =
    {0, 1, 2, 3, 4, 5, 6, 7, 8};  // For the base mesh, the parents are always themselves
  mapData.m_shapeTypes = Vec<mir::Shape>(numElements, mir::Shape::Quad);

  // Build the mesh
  mir::MIRMesh testMesh;
  testMesh.initializeMesh(verts, elems, numMaterials, topoData, mapData, volFracs);

  return testMesh;
}

//--------------------------------------------------------------------------------

template <typename TopoView, typename CoordsetView>
static void addCircleMaterial(const TopoView& topoView,
                              const CoordsetView& coordsetView,
                              conduit::Node& mesh,
                              const mir::Point2& circleCenter,
                              axom::float64 circleRadius,
                              int numSamples)
{
  constexpr int GREEN = 0;
  constexpr int BLUE = 1;

  int numElements = topoView.numberOfZones();
  std::vector<float> volFracs[2];
  volFracs[GREEN].resize(numElements);
  volFracs[BLUE].resize(numElements);

  // Generate the element volume fractions for the circle
  axom::ArrayView<float> greenView(volFracs[GREEN].data(), numElements);
  axom::ArrayView<float> blueView(volFracs[BLUE].data(), numElements);
  typename CoordsetView::PointType center;
  center[0] = circleCenter[0];
  center[1] = circleCenter[1];
#if 1
  const TopoView deviceTopologyView(topoView);
  axom::for_all<axom::SEQ_EXEC>(topoView.numberOfZones(), AXOM_LAMBDA(axom::IndexType zoneIndex)
  {
    const auto zone = deviceTopologyView.zone(zoneIndex);
    auto vf = calculatePercentOverlapMonteCarlo(numSamples,
                                                center,
                                                circleRadius,
                                                coordsetView[zone.getId(0)],
                                                coordsetView[zone.getId(1)],
                                                coordsetView[zone.getId(2)],
                                                coordsetView[zone.getId(3)]);
    greenView[zoneIndex] = vf;
    blueView[zoneIndex] = 1.0 - vf;
  });
#else
  topoView.template for_all_zones<axom::SEQ_EXEC>(
    AXOM_LAMBDA(auto zoneIndex, const auto& zone) {
      auto vf = calculatePercentOverlapMonteCarlo(numSamples,
                                                  center,
                                                  circleRadius,
                                                  coordsetView[zone.getId(0)],
                                                  coordsetView[zone.getId(1)],
                                                  coordsetView[zone.getId(2)],
                                                  coordsetView[zone.getId(3)]);
      greenView[zoneIndex] = vf;
      blueView[zoneIndex] = 1.0 - vf;
    });
#endif
  // Figure out the material buffers from the volume fractions.
  std::vector<int> material_ids, sizes, offsets, indices;
  std::vector<int> volume_fractions;
  for(int i = 0; i < numElements; ++i)
  {
    int nmats = 0;
    offsets.push_back(indices.size());
    if(volFracs[GREEN][i] > 0.)
    {
      material_ids.push_back(GREEN);
      volume_fractions.push_back(volFracs[GREEN][i]);
      indices.push_back(indices.size());
      nmats++;
    }
    if(volFracs[BLUE][i] > 0.)
    {
      material_ids.push_back(BLUE);
      volume_fractions.push_back(volFracs[BLUE][i]);
      indices.push_back(indices.size());
      nmats++;
    }
    sizes.push_back(nmats);
  }

  mesh["matsets/mat/topology"] = "mesh";
  mesh["matsets/mat/material_map/green"] = GREEN;
  mesh["matsets/mat/material_map/blue"] = BLUE;
  mesh["matsets/mat/material_ids"].set(material_ids);
  mesh["matsets/mat/volume_fractions"].set(volume_fractions);
  mesh["matsets/mat/sizes"].set(sizes);
  mesh["matsets/mat/offsets"].set(offsets);
  mesh["matsets/mat/indices"].set(indices);
}

//--------------------------------------------------------------------------------
void MeshTester::initTestCaseFour(conduit::Node& mesh)
{
  mesh3x3(mesh);

  // Make views
  using CoordsetView = axom::mir::views::ExplicitCoordsetView<float, 2>;
  CoordsetView coordsetView(
    bputils::make_array_view<float>(mesh["coordsets/coords/values/x"]),
    bputils::make_array_view<float>(mesh["coordsets/coords/values/y"]));
  using TopoView = axom::mir::views::UnstructuredTopologySingleShapeView<
    axom::mir::views::QuadShape<int>>;
  TopoView topoView(bputils::make_array_view<int>(
    mesh["topologies/mesh/elements/connectivity"]));

  // Add material
  const auto circleCenter = mir::Point2::make_point(1.5, 1.5);
  const axom::float64 circleRadius = 1.25;
  const int numSamples = 100;
  addCircleMaterial<TopoView, CoordsetView>(topoView,
                                            coordsetView,
                                            mesh,
                                            circleCenter,
                                            circleRadius,
                                            numSamples);
}

//--------------------------------------------------------------------------------

mir::MIRMesh MeshTester::createUniformGridTestCaseMesh(
  int gridSize,
  const mir::Point2& circleCenter,
  axom::float64 circleRadius)
{
  // Generate the mesh topology
  mir::CellData cellData = generateGrid(gridSize);

  mir::VertSet verts =
    mir::VertSet(cellData.m_numVerts);  // Construct the vertex set
  mir::ElemSet elems =
    mir::ElemSet(cellData.m_numElems);  // Construct the element set

  int numMaterials = 2;
  enum
  {
    GREEN = 0,
    BLUE = 1
  };

  VolumeFractions volFracs;
  volFracs.resize(numMaterials);
  volFracs[GREEN].resize(cellData.m_numElems);
  volFracs[BLUE].resize(cellData.m_numElems);

  // Generate the element volume fractions for the circle
  const int numMonteCarloSamples = 100;
  auto& pos = cellData.m_mapData.m_vertexPositions;
  const auto& evInds = cellData.m_topology.m_evInds;
  for(int i = 0; i < cellData.m_numElems; ++i)
  {
    auto vf = calculatePercentOverlapMonteCarlo(numMonteCarloSamples,
                                                circleCenter,
                                                circleRadius,
                                                pos[evInds[i * 4 + 0]],
                                                pos[evInds[i * 4 + 1]],
                                                pos[evInds[i * 4 + 2]],
                                                pos[evInds[i * 4 + 3]]);
    volFracs[GREEN][i] = vf;
    volFracs[BLUE][i] = 1.0 - vf;
  }

  cellData.m_mapData.m_elementDominantMaterials =
    Vec<int>(cellData.m_numVerts, NULL_MAT);
  cellData.m_mapData.m_shapeTypes =
    Vec<mir::Shape>(cellData.m_numVerts, mir::Shape::Quad);
  cellData.m_mapData.m_elementParents.resize(cellData.m_numVerts);
  for(auto i : elems.positions())
  {
    cellData.m_mapData.m_elementParents[i] = i;
  }

  // Build the mesh
  mir::MIRMesh testMesh;
  testMesh.initializeMesh(verts,
                          elems,
                          numMaterials,
                          cellData.m_topology,
                          cellData.m_mapData,
                          volFracs);

  return testMesh;
}

//--------------------------------------------------------------------------------
void MeshTester::createUniformGridTestCaseMesh(int gridSize,
                                               const mir::Point2& circleCenter,
                                               axom::float64 circleRadius,
                                               conduit::Node& mesh)
{
  // Generate the mesh
  generateGrid(gridSize, mesh);

  // Make views
  using CoordsetView = axom::mir::views::ExplicitCoordsetView<float, 2>;
  CoordsetView coordsetView(
    bputils::make_array_view<float>(mesh["coordsets/coords/values/x"]),
    bputils::make_array_view<float>(mesh["coordsets/coords/values/y"]));
  using TopoView = axom::mir::views::UnstructuredTopologySingleShapeView<
    axom::mir::views::QuadShape<int>>;
  TopoView topoView(bputils::make_array_view<int>(
    mesh["topologies/mesh/elements/connectivity"]));

  // Add material
  int numSamples = 100;
  addCircleMaterial<TopoView, CoordsetView>(topoView,
                                            coordsetView,
                                            mesh,
                                            circleCenter,
                                            circleRadius,
                                            numSamples);
}

//--------------------------------------------------------------------------------

mir::CellData MeshTester::generateGrid(int gridSize)
{
  // Generate the topology for a uniform quad mesh with n x n elements automatically
  int numElements = gridSize * gridSize;
  int numVertices = (gridSize + 1) * (gridSize + 1);

  mir::CellData data;

  data.m_numVerts = numVertices;
  data.m_numElems = numElements;

  // Generate the evInds
  auto& evInds = data.m_topology.m_evInds;
  for(int eID = 0; eID < numElements; ++eID)
  {
    int row = eID / gridSize;  // note the integer division
    int vertsPerRow = gridSize + 1;
    int elemsPerRow = gridSize;

    evInds.push_back((eID % elemsPerRow) + row * vertsPerRow + 0);
    evInds.push_back((eID % elemsPerRow) + (row + 1) * vertsPerRow + 0);
    evInds.push_back((eID % elemsPerRow) + (row + 1) * vertsPerRow + 1);
    evInds.push_back((eID % elemsPerRow) + row * vertsPerRow + 1);
  }

  // Generate the evBegins
  auto& evBegins = data.m_topology.m_evBegins;
  evBegins.push_back(0);
  for(int i = 0; i < numElements; ++i)
  {
    evBegins.push_back((i + 1) * 4);
  }

  // Generate the veInds
  auto& veInds = data.m_topology.m_veInds;
  auto& veBegins = data.m_topology.m_veBegins;
  std::map<int, std::vector<int>> veInds_data;
  for(int evInd_itr = 0; evInd_itr < numElements * 4; ++evInd_itr)
  {
    int currentElementID = evInd_itr / 4;  // note the integer division
    veInds_data[evInds[evInd_itr]].push_back(currentElementID);
  }

  for(auto itr = veInds_data.begin(); itr != veInds_data.end(); itr++)
  {
    // Sort the vector
    std::sort(itr->second.begin(), itr->second.end());

    // Add the elements associated with the current vertex to veInds
    for(unsigned long i = 0; i < itr->second.size(); ++i)
      veInds.push_back(itr->second[i]);
  }

  // Generate the veBegins
  veBegins.push_back(0);
  int currentIndexCount = 0;
  for(auto itr = veInds_data.begin(); itr != veInds_data.end(); itr++)
  {
    currentIndexCount += itr->second.size();
    veBegins.push_back(currentIndexCount);
  }

  // Generate the vertex positions
  auto& points = data.m_mapData.m_vertexPositions;
  for(int y = gridSize; y > -1; --y)
  {
    for(int x = 0; x < gridSize + 1; ++x)
    {
      points.push_back(mir::Point2::make_point(x, y));
    }
  }

  // // Print out the results
  // printf("evInds: { ");
  // for (int i = 0; i < evInds.size(); i++)
  // {
  //   printf("%d ", evInds[i]);
  //   if ((i+1) % 4 == 0 && i != 0)
  //     printf("\n");
  // }
  // printf("}\n");

  // printf("evBegins: { ");
  // for (int i = 0; i < evBegins.size(); i++)
  // {
  //   printf("%d ", evBegins[i]);
  // }
  // printf("}\n");

  // printf("veInds: { ");
  // for (int i = 0; i < veInds.size(); i++)
  // {
  //   printf("%d ", veInds[i]);
  // }
  // printf("}\n");

  // printf("veBegins: { ");
  // for (int i = 0; i < veBegins.size(); i++)
  // {
  //   printf("%d ", veBegins[i]);
  // }
  // printf("}\n");

  // printf("points: { ");
  // for (int i = 0; i < numVertices; ++i)
  // {
  //   printf("{%.2f, %.2f} ", points[i][0], points[i][1]);
  // }
  // printf("}\n");

  return data;
}

void MeshTester::generateGrid(int gridSize, conduit::Node& mesh)
{
  int nx = gridSize + 1;
  int ny = gridSize + 1;
  int nzones = gridSize * gridSize;
  int nnodes = nx * ny;

  std::vector<float> xc, yc;
  xc.reserve(nnodes);
  yc.reserve(nnodes);
  for(int j = 0; j < ny; j++)
  {
    for(int i = 0; i < nx; i++)
    {
      xc.push_back(i);
      yc.push_back(j);
    }
  }

  std::vector<int> conn, sizes, offsets;
  conn.reserve(nzones * 4);
  sizes.reserve(nzones);
  offsets.reserve(nzones);
  for(int j = 0; j < gridSize; j++)
  {
    for(int i = 0; i < gridSize; i++)
    {
      offsets.push_back(offsets.size() * 4);
      sizes.push_back(4);
      conn.push_back(j * nx + i);
      conn.push_back(j * nx + i + 1);
      conn.push_back((j + 1) * nx + i + 1);
      conn.push_back((j + 1) * nx + i);
    }
  }

  mesh["coordsets/coords/type"] = "explicit";
  mesh["coordsets/coords/values/x"].set(xc);
  mesh["coordsets/coords/values/y"].set(yc);
  mesh["topologies/mesh/type"] = "unstructured";
  mesh["topologies/mesh/coordset"] = "coords";
  mesh["topologies/mesh/elements/shape"] = "quad";
  mesh["topologies/mesh/elements/connectivity"].set(conn);
  mesh["topologies/mesh/elements/sizes"].set(sizes);
  mesh["topologies/mesh/elements/offsets"].set(offsets);
}

void MeshTester::generateGrid3D(int gridSize, conduit::Node& mesh)
{
  int nx = gridSize + 1;
  int ny = gridSize + 1;
  int nz = gridSize + 1;
  int nzones = gridSize * gridSize * gridSize;
  int nnodes = nx * ny * nz;

  std::vector<float> xc, yc, zc;
  xc.reserve(nnodes);
  yc.reserve(nnodes);
  zc.reserve(nnodes);
  for(int k = 0; k < nz; k++)
  {
    for(int j = 0; j < ny; j++)
    {
      for(int i = 0; i < nx; i++)
      {
        xc.push_back(i);
        yc.push_back(j);
        zc.push_back(k);
      }
    }
  }

  std::vector<int> conn, sizes, offsets;
  conn.reserve(nzones * 8);
  sizes.reserve(nzones);
  offsets.reserve(nzones);
  for(int k = 0; k < gridSize; k++)
  {
    for(int j = 0; j < gridSize; j++)
    {
      for(int i = 0; i < gridSize; i++)
      {
        offsets.push_back(offsets.size() * 8);
        sizes.push_back(8);
        conn.push_back((k * nx * ny) + (j * nx) + i);
        conn.push_back((k * nx * ny) + (j * nx) + i + 1);
        conn.push_back((k * nx * ny) + ((j + 1) * nx) + i + 1);
        conn.push_back((k * nx * ny) + ((j + 1) * nx) + i);
        conn.push_back(((k + 1) * nx * ny) + (j * nx) + i);
        conn.push_back(((k + 1) * nx * ny) + (j * nx) + i + 1);
        conn.push_back(((k + 1) * nx * ny) + ((j + 1) * nx) + i + 1);
        conn.push_back(((k + 1) * nx * ny) + ((j + 1) * nx) + i);
      }
    }
  }

  mesh["coordsets/coords/type"] = "explicit";
  mesh["coordsets/coords/values/x"].set(xc);
  mesh["coordsets/coords/values/y"].set(yc);
  mesh["coordsets/coords/values/z"].set(zc);
  mesh["topologies/mesh/type"] = "unstructured";
  mesh["topologies/mesh/coordset"] = "coords";
  mesh["topologies/mesh/elements/shape"] = "hex";
  mesh["topologies/mesh/elements/connectivity"].set(conn);
  mesh["topologies/mesh/elements/sizes"].set(sizes);
  mesh["topologies/mesh/elements/offsets"].set(offsets);
}
//--------------------------------------------------------------------------------

mir::MIRMesh MeshTester::initTestCaseFive(int gridSize, int numCircles)
{
  // Generate the mesh topology
  mir::CellData cellData = generateGrid(gridSize);

  mir::VertSet verts =
    mir::VertSet(cellData.m_numVerts);  // Construct the vertex set
  mir::ElemSet elems =
    mir::ElemSet(cellData.m_numElems);  // Construct the element set

  // Generate the element volume fractions with concentric circles
  int numMaterials = numCircles + 1;
  int defaultMaterialID =
    numMaterials - 1;  // default material is always the last index

  mir::Point2 circleCenter = mir::Point2::make_point(
    gridSize / 2.0,
    gridSize / 2.0);  // all circles are centered around the same point

  // Initialize the radii of the circles
  std::vector<axom::float64> circleRadii;
  axom::float64 maxRadius =
    gridSize / 2.4;  // Note: The choice of divisor is arbitrary
  axom::float64 minRadius =
    gridSize / 8;  // Note: The choice of divisor is arbitrary

  axom::float64 radiusDelta;
  if(numCircles <= 1)
    radiusDelta = (maxRadius - minRadius);
  else
    radiusDelta = (maxRadius - minRadius) / (double)(numCircles - 1);

  for(int i = 0; i < numCircles; ++i)
  {
    circleRadii.push_back(minRadius + (i * radiusDelta));
  }

  // Initialize all material volume fractions to 0
  std::vector<std::vector<axom::float64>> materialVolumeFractionsData;
  for(int i = 0; i < numMaterials; ++i)
  {
    std::vector<axom::float64> tempVec;
    tempVec.resize(cellData.m_numElems);
    materialVolumeFractionsData.push_back(tempVec);
  }

  // Use the uniform sampling method to generate volume fractions for each material
  // Note: Assumes that the cell is a parallelogram. This could be modified via biliear interpolation
  for(int eID = 0; eID < cellData.m_numElems; ++eID)
  {
    mir::Point2& v0 =
      cellData.m_mapData
        .m_vertexPositions[cellData.m_topology.m_evInds[eID * 4 + 0]];
    mir::Point2& v1 =
      cellData.m_mapData
        .m_vertexPositions[cellData.m_topology.m_evInds[eID * 4 + 1]];
    mir::Point2& v2 =
      cellData.m_mapData
        .m_vertexPositions[cellData.m_topology.m_evInds[eID * 4 + 2]];
    //mir::Point2& v3 = cellData.m_mapData.m_vertexPositions[cellData.m_topology.m_evInds[eID * 4 + 3]];

    // Run the uniform sampling to determine how much of the current cell is composed of each material
    int materialCount[numMaterials];
    for(int i = 0; i < numMaterials; ++i) materialCount[i] = 0;

    for(int matID = 0; matID < numMaterials; ++matID)
    {
      materialVolumeFractionsData[matID][eID] =
        materialCount[matID] / (double)(gridSize * gridSize);
    }

    axom::float64 delta_x =
      axom::utilities::abs(v2[0] - v1[1]) / (double)(gridSize - 1);
    axom::float64 delta_y =
      axom::utilities::abs(v0[1] - v1[1]) / (double)(gridSize - 1);

    for(int y = 0; y < gridSize; ++y)
    {
      for(int x = 0; x < gridSize; ++x)
      {
        mir::Point2 samplePoint =
          mir::Point2::make_point(delta_x * x + v1[0], delta_y * y + v1[1]);
        bool isPointSampled = false;
        for(int cID = 0; cID < numCircles && !isPointSampled; ++cID)
        {
          const auto r = circleRadii[cID];
          if(primal::squared_distance(samplePoint, circleCenter) < r * r)
          {
            materialCount[cID]++;
            isPointSampled = true;
          }
        }
        if(!isPointSampled)
        {
          // The point was not within any of the circles, so increment the count for the default material
          materialCount[defaultMaterialID]++;
        }
      }
    }

    // Assign the element volume fractions based on the count of the samples in each circle
    for(int matID = 0; matID < numMaterials; ++matID)
    {
      materialVolumeFractionsData[matID][eID] =
        materialCount[matID] / (double)(gridSize * gridSize);
    }
  }

  std::vector<int> elementParents;  // For the base mesh, the parents are always themselves
  std::vector<int> elementDominantMaterials;
  std::vector<mir::Shape> elementShapeTypes;
  for(int i = 0; i < cellData.m_numElems; ++i)
  {
    elementParents.push_back(i);
    elementDominantMaterials.push_back(NULL_MAT);
    elementShapeTypes.push_back(mir::Shape::Quad);
  }

  CellTopologyData topology;
  topology.m_evInds = cellData.m_topology.m_evInds;
  topology.m_evBegins = cellData.m_topology.m_evBegins;
  topology.m_veInds = cellData.m_topology.m_veInds;
  topology.m_veBegins = cellData.m_topology.m_veBegins;

  CellMapData mapData;
  mapData.m_elementDominantMaterials = elementDominantMaterials;
  mapData.m_elementParents = elementParents;
  mapData.m_vertexPositions = cellData.m_mapData.m_vertexPositions;
  mapData.m_shapeTypes = elementShapeTypes;

  // Build the mesh
  mir::MIRMesh testMesh;
  testMesh.initializeMesh(verts,
                          elems,
                          numMaterials,
                          topology,
                          mapData,
                          materialVolumeFractionsData);

  return testMesh;
}

template <typename TopoView, typename CoordsetView>
void addConcentricCircleMaterial(const TopoView& topoView,
                                 const CoordsetView& coordsetView,
                                 const mir::Point2& circleCenter,
                                 std::vector<axom::float64>& circleRadii,
                                 int numSamples,
                                 conduit::Node& mesh)
{
  // Generate the element volume fractions with concentric circles
  int numMaterials = circleRadii.size() + 1;
  int defaultMaterialID =
    numMaterials - 1;  // default material is always the last index

  // Initialize all material volume fractions to 0
  std::vector<std::vector<axom::float64>> materialVolumeFractionsData(
    numMaterials);
  constexpr int MAXMATERIALS = 100;
  axom::StackArray<axom::ArrayView<axom::float64>, MAXMATERIALS> matvfViews;
  for(int i = 0; i < numMaterials; ++i)
  {
    const auto len = topoView.numberOfZones();
    materialVolumeFractionsData[i].resize(len, 0.);
    matvfViews[i] =
      axom::ArrayView<axom::float64>(materialVolumeFractionsData[i].data(), len);
  }
  auto circleRadiiView =
    axom::ArrayView<axom::float64>(circleRadii.data(), circleRadii.size());
  const int numCircles = circleRadii.size();

  // Use the uniform sampling method to generate volume fractions for each material
  // Note: Assumes that the cell is a parallelogram. This could be modified via biliear interpolation
#if 1
  const TopoView deviceTopologyView(topoView);
  axom::for_all<axom::SEQ_EXEC>(topoView.numberOfZones(),
    AXOM_LAMBDA(axom::IndexType eID) {
      const auto zone = deviceTopologyView.zone(eID);
#else
  topoView.template for_all_zones<axom::SEQ_EXEC>(
    AXOM_LAMBDA(auto eID, const auto& zone) {
#endif
      auto v0 = coordsetView[zone.getId(0)];
      auto v1 = coordsetView[zone.getId(1)];
      auto v2 = coordsetView[zone.getId(2)];

      // Run the uniform sampling to determine how much of the current cell is composed of each material
      int materialCount[numMaterials];
      for(int i = 0; i < numMaterials; ++i) materialCount[i] = 0;

      axom::float64 delta_x =
        axom::utilities::abs(v1[0] - v0[0]) / (axom::float64)(numSamples - 1);
      axom::float64 delta_y =
        axom::utilities::abs(v2[1] - v1[1]) / (axom::float64)(numSamples - 1);

      for(int y = 0; y < numSamples; ++y)
      {
        for(int x = 0; x < numSamples; ++x)
        {
          mir::Point2 samplePoint =
            mir::Point2::make_point(delta_x * x + v0[0], delta_y * y + v0[1]);
          bool isPointSampled = false;
          for(int cID = 0; cID < numCircles && !isPointSampled; ++cID)
          {
            const auto r = circleRadiiView[cID];
            if(primal::squared_distance(samplePoint, circleCenter) < r * r)
            {
              materialCount[cID]++;
              isPointSampled = true;
            }
          }
          if(!isPointSampled)
          {
            // The point was not within any of the circles, so increment the count for the default material
            materialCount[defaultMaterialID]++;
          }
        }
      }

      // Assign the element volume fractions based on the count of the samples in each circle
      for(int matID = 0; matID < numMaterials; ++matID)
      {
        matvfViews[matID][eID] =
          materialCount[matID] / (axom::float64)(numSamples * numSamples);
      }
    });

  // Figure out the material buffers from the volume fractions.
  std::vector<int> material_ids, sizes, offsets, indices;
  std::vector<float> volume_fractions;
  const axom::IndexType numElements = topoView.numberOfZones();
  std::vector<float> sums(numMaterials, 0.);
  for(axom::IndexType i = 0; i < numElements; ++i)
  {
    int nmats = 0;
    offsets.push_back(indices.size());
    for(int mat = 0; mat < numMaterials; mat++)
    {
      if(materialVolumeFractionsData[mat][i] > 0.)
      {
        material_ids.push_back(mat);
        volume_fractions.push_back(materialVolumeFractionsData[mat][i]);
        indices.push_back(indices.size());
        nmats++;

        // Keep a total of the VFs for each material.
        sums[mat] += materialVolumeFractionsData[mat][i];
      }
    }
    sizes.push_back(nmats);
  }

  // Add the material
  mesh["matsets/mat/topology"] = "mesh";
  for(int mat = 0; mat < numMaterials; mat++)
  {
    if(sums[mat] > 0.)
    {
      std::stringstream ss;
      ss << "matsets/mat/material_map/mat" << mat;
      mesh[ss.str()] = mat;
    }
  }
  mesh["matsets/mat/material_ids"].set(material_ids);
  mesh["matsets/mat/volume_fractions"].set(volume_fractions);
  mesh["matsets/mat/sizes"].set(sizes);
  mesh["matsets/mat/offsets"].set(offsets);
  mesh["matsets/mat/indices"].set(indices);
}

void MeshTester::initTestCaseFive(int gridSize, int numCircles, conduit::Node& mesh)
{
  // Generate the mesh topology
  generateGrid(gridSize, mesh);

  mir::Point2 circleCenter = mir::Point2::make_point(
    gridSize / 2.0,
    gridSize / 2.0);  // all circles are centered around the same point

  // Initialize the radii of the circles
  std::vector<axom::float64> circleRadii;
  axom::float64 maxRadius =
    gridSize / 2.4;  // Note: The choice of divisor is arbitrary
  axom::float64 minRadius =
    gridSize / 8;  // Note: The choice of divisor is arbitrary

  axom::float64 radiusDelta;
  if(numCircles <= 1)
    radiusDelta = (maxRadius - minRadius);
  else
    radiusDelta = (maxRadius - minRadius) / (double)(numCircles - 1);

  for(int i = 0; i < numCircles; ++i)
  {
    circleRadii.push_back(minRadius + (i * radiusDelta));
  }

  // Make views
  using CoordsetView = axom::mir::views::ExplicitCoordsetView<float, 2>;
  CoordsetView coordsetView(
    bputils::make_array_view<float>(mesh["coordsets/coords/values/x"]),
    bputils::make_array_view<float>(mesh["coordsets/coords/values/y"]));
  using TopoView = axom::mir::views::UnstructuredTopologySingleShapeView<
    axom::mir::views::QuadShape<int>>;
  TopoView topoView(bputils::make_array_view<int>(
    mesh["topologies/mesh/elements/connectivity"]));

  // Add the material
  addConcentricCircleMaterial<TopoView, CoordsetView>(topoView,
                                                      coordsetView,
                                                      circleCenter,
                                                      circleRadii,
                                                      100,
                                                      mesh);
}

//--------------------------------------------------------------------------------

int MeshTester::circleQuadCornersOverlaps(const mir::Point2& circleCenter,
                                          axom::float64 circleRadius,
                                          const mir::Point2& quadP0,
                                          const mir::Point2& quadP1,
                                          const mir::Point2& quadP2,
                                          const mir::Point2& quadP3)
{
  // Check if any of the quad's corners are within the circle
  auto d0Sq = primal::squared_distance(quadP0, circleCenter);
  auto d1Sq = primal::squared_distance(quadP1, circleCenter);
  auto d2Sq = primal::squared_distance(quadP2, circleCenter);
  auto d3Sq = primal::squared_distance(quadP3, circleCenter);
  auto dRSq = circleRadius * circleRadius;

  int numCorners = 0;

  if(d0Sq < dRSq) numCorners++;
  if(d1Sq < dRSq) numCorners++;
  if(d2Sq < dRSq) numCorners++;
  if(d3Sq < dRSq) numCorners++;

  return numCorners;
}

//--------------------------------------------------------------------------------

mir::MIRMesh MeshTester::initQuadClippingTestMesh()
{
  // Generate the mesh topology
  int gridSize = 3;
  mir::CellData cellData = generateGrid(gridSize);

  mir::VertSet verts =
    mir::VertSet(cellData.m_numVerts);  // Construct the vertex set
  mir::ElemSet elems =
    mir::ElemSet(cellData.m_numElems);  // Construct the element set

  int numMaterials = 2;

  std::vector<std::vector<axom::float64>> elementVF;
  elementVF.resize(numMaterials);
  elementVF[0] = {1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.0, 0.0, 0.0};
  elementVF[1] = {0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0};

  std::vector<int> elementParents;
  std::vector<int> elementDominantMaterials;
  std::vector<mir::Shape> elementShapeTypes;
  for(int i = 0; i < cellData.m_numElems; ++i)
  {
    elementParents.push_back(i);
    elementDominantMaterials.push_back(NULL_MAT);
    elementShapeTypes.push_back(mir::Shape::Quad);
  }

  cellData.m_mapData.m_elementDominantMaterials = elementDominantMaterials;
  cellData.m_mapData.m_elementParents = elementParents;
  cellData.m_mapData.m_shapeTypes = elementShapeTypes;

  // Build the mesh
  mir::MIRMesh testMesh;
  testMesh.initializeMesh(verts,
                          elems,
                          numMaterials,
                          cellData.m_topology,
                          cellData.m_mapData,
                          elementVF);

  return testMesh;
}

void MeshTester::initQuadClippingTestMesh(conduit::Node& mesh)
{
  // Generate the mesh topology
  constexpr int gridSize = 3;
  generateGrid(gridSize, mesh);

  // clang-format off
  mesh["matsets/mat/topology"] = "mesh";
  mesh["matsets/mat/material_map/green"] = 0;
  mesh["matsets/mat/material_map/blue"] = 1;
  mesh["matsets/mat/material_ids"].set(std::vector<int>{{
   0,
   0,
   0,
   0, 1,
   0, 1,
   0, 1,
   1,
   1,
   1
  }});
  mesh["matsets/mat/volume_fractions"].set(std::vector<float>{{
   1.,
   1.,
   1.,
   0.5, 0.5,
   0.5, 0.5,
   0.5, 0.5,
   1.,
   1.,
   1.
  }});
  mesh["matsets/mat/sizes"].set(std::vector<int>{{
    1, 1, 1, 2, 2, 2, 1, 1, 1
  }});
  mesh["matsets/mat/offsets"].set(std::vector<int>{{
    0, 1, 2, 3, 5, 7, 9, 10, 11
  }});
  mesh["matsets/mat/indices"].set(std::vector<int>{{
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11
  }});
  // clang-format on
}
//--------------------------------------------------------------------------------

mir::MIRMesh MeshTester::initTestCaseSix(int gridSize, int numSpheres)
{
  // Generate the mesh topology
  mir::CellData cellData = generateGrid3D(gridSize);

  mir::VertSet verts =
    mir::VertSet(cellData.m_numVerts);  // Construct the vertex set
  mir::ElemSet elems =
    mir::ElemSet(cellData.m_numElems);  // Construct the element set

  // // Generate the element volume fractions with concentric spheres
  int numMaterials = numSpheres + 1;
  int defaultMaterialID =
    numMaterials - 1;  // default material is always the last index

  mir::Point2 sphereCenter = mir::Point2::make_point(
    gridSize / 2.0,
    gridSize / 2.0,
    gridSize / 2.0);  // all spheres are centered around the same point

  // Initialize the radii of the circles
  std::vector<axom::float64> sphereRadii;
  axom::float64 maxRadius =
    gridSize / 2.0;  // Note: The choice of divisor is arbitrary
  axom::float64 minRadius =
    gridSize / 4.0;  // Note: The choice of divisor is arbitrary

  axom::float64 radiusDelta;
  if(numSpheres <= 1)
    radiusDelta = (maxRadius - minRadius);
  else
    radiusDelta = (maxRadius - minRadius) / (double)(numSpheres - 1);

  for(int i = 0; i < numSpheres; ++i)
  {
    auto rad = minRadius + (i * radiusDelta);
    sphereRadii.push_back(rad * rad);
  }

  // Initialize all material volume fractions to 0
  std::vector<std::vector<axom::float64>> materialVolumeFractionsData;
  for(int i = 0; i < numMaterials; ++i)
  {
    std::vector<axom::float64> tempVec;
    tempVec.resize(cellData.m_numElems, 0);
    materialVolumeFractionsData.push_back(tempVec);
  }

  // Use the uniform sampling method to generate volume fractions for each material
  for(int eID = 0; eID < cellData.m_numElems; ++eID)
  {
    mir::Point2 v0 =
      cellData.m_mapData
        .m_vertexPositions[cellData.m_topology.m_evInds[eID * 8 + 0]];
    mir::Point2 v1 =
      cellData.m_mapData
        .m_vertexPositions[cellData.m_topology.m_evInds[eID * 8 + 1]];
    mir::Point2 v2 =
      cellData.m_mapData
        .m_vertexPositions[cellData.m_topology.m_evInds[eID * 8 + 2]];
    mir::Point2 v3 =
      cellData.m_mapData
        .m_vertexPositions[cellData.m_topology.m_evInds[eID * 8 + 3]];

    mir::Point2 v4 =
      cellData.m_mapData
        .m_vertexPositions[cellData.m_topology.m_evInds[eID * 8 + 4]];
    mir::Point2 v5 =
      cellData.m_mapData
        .m_vertexPositions[cellData.m_topology.m_evInds[eID * 8 + 5]];
    mir::Point2 v6 =
      cellData.m_mapData
        .m_vertexPositions[cellData.m_topology.m_evInds[eID * 8 + 6]];
    mir::Point2 v7 =
      cellData.m_mapData
        .m_vertexPositions[cellData.m_topology.m_evInds[eID * 8 + 7]];

    // Run the uniform sampling to determine how much of the current cell is composed of each material
    int materialCount[numMaterials];
    for(int i = 0; i < numMaterials; ++i) materialCount[i] = 0;

    axom::float64 delta_x =
      axom::utilities::abs(v2[0] - v1[0]) / (double)(gridSize - 1);
    axom::float64 delta_y =
      axom::utilities::abs(v0[1] - v1[1]) / (double)(gridSize - 1);
    axom::float64 delta_z =
      axom::utilities::abs(v5[2] - v1[2]) / (double)(gridSize - 1);

    for(int z = 0; z < gridSize; ++z)
    {
      for(int y = 0; y < gridSize; ++y)
      {
        for(int x = 0; x < gridSize; ++x)
        {
          mir::Point2 samplePoint = mir::Point2::make_point(delta_x * x + v1[0],
                                                            delta_y * y + v1[1],
                                                            delta_z * z + v1[2]);
          bool isPointSampled = false;
          for(int cID = 0; cID < numSpheres && !isPointSampled; ++cID)
          {
            if(primal::squared_distance(samplePoint, sphereCenter) <
               sphereRadii[cID])
            {
              materialCount[cID]++;
              isPointSampled = true;
            }
          }
          if(!isPointSampled)
          {
            // The point was not within any of the circles, so increment the count for the default material
            materialCount[defaultMaterialID]++;
          }
        }
      }
    }

    // Assign the element volume fractions based on the count of the samples in each circle
    for(int matID = 0; matID < numMaterials; ++matID)
    {
      materialVolumeFractionsData[matID][eID] =
        materialCount[matID] / (double)(gridSize * gridSize * gridSize);
    }
  }

  std::vector<int> elementParents;  // For the base mesh, the parents are always themselves
  std::vector<int> elementDominantMaterials;
  std::vector<mir::Shape> elementShapeTypes;
  for(int i = 0; i < cellData.m_numElems; ++i)
  {
    elementParents.push_back(i);
    elementDominantMaterials.push_back(NULL_MAT);
    elementShapeTypes.push_back(mir::Shape::Hexahedron);
  }

  CellTopologyData topology;
  topology.m_evInds = cellData.m_topology.m_evInds;
  topology.m_evBegins = cellData.m_topology.m_evBegins;
  topology.m_veInds = cellData.m_topology.m_veInds;
  topology.m_veBegins = cellData.m_topology.m_veBegins;

  CellMapData mapData;
  mapData.m_elementDominantMaterials = elementDominantMaterials;
  mapData.m_elementParents = elementParents;
  mapData.m_vertexPositions = cellData.m_mapData.m_vertexPositions;
  mapData.m_shapeTypes = elementShapeTypes;

  // Build the mesh
  mir::MIRMesh testMesh;
  testMesh.initializeMesh(verts,
                          elems,
                          numMaterials,
                          topology,
                          mapData,
                          materialVolumeFractionsData);

  return testMesh;
}

//--------------------------------------------------------------------------------

mir::CellData MeshTester::generateGrid3D(int gridSize)
{
  // Generate the topology for a uniform quad mesh with n x n elements automatically
  int numElements = gridSize * gridSize * gridSize;
  int numVertices = (gridSize + 1) * (gridSize + 1) * (gridSize + 1);

  // Generate the evInds
  std::vector<mir::PosType> evInds;
  for(int eID = 0; eID < numElements; ++eID)
  {
    int vertsPerLayer = (gridSize + 1) * (gridSize + 1);
    int elemsPerLayer = gridSize * gridSize;
    int layer = eID / (gridSize * gridSize);  // note the integer division

    int vertsPerRow = gridSize + 1;
    int elemsPerRow = gridSize;
    int row = (eID % elemsPerLayer) / gridSize;  // note the integer division

    // front face of the grid cube
    evInds.push_back((eID % elemsPerLayer % elemsPerRow) + (vertsPerRow * row) +
                     0 + (vertsPerLayer * layer));
    evInds.push_back((eID % elemsPerLayer % elemsPerRow) +
                     (vertsPerRow * (row + 1)) + 0 + (vertsPerLayer * layer));
    evInds.push_back((eID % elemsPerLayer % elemsPerRow) +
                     (vertsPerRow * (row + 1)) + 1 + (vertsPerLayer * layer));
    evInds.push_back((eID % elemsPerLayer % elemsPerRow) + (vertsPerRow * row) +
                     1 + (vertsPerLayer * layer));

    // back face of the grid cube
    evInds.push_back((eID % elemsPerLayer % elemsPerRow) + (vertsPerRow * row) +
                     0 + (vertsPerLayer * (layer + 1)));
    evInds.push_back((eID % elemsPerLayer % elemsPerRow) +
                     (vertsPerRow * (row + 1)) + 0 +
                     (vertsPerLayer * (layer + 1)));
    evInds.push_back((eID % elemsPerLayer % elemsPerRow) +
                     (vertsPerRow * (row + 1)) + 1 +
                     (vertsPerLayer * (layer + 1)));
    evInds.push_back((eID % elemsPerLayer % elemsPerRow) + (vertsPerRow * row) +
                     1 + (vertsPerLayer * (layer + 1)));
  }

  // Generate the evBegins
  std::vector<mir::PosType> evBegins;
  evBegins.push_back(0);
  for(int i = 0; i < numElements; ++i)
  {
    evBegins.push_back((i + 1) * 8);
  }

  // Generate the veInds
  std::map<int, std::vector<int>> veInds_data;
  std::vector<mir::PosType> veInds;
  for(int evInd_itr = 0; evInd_itr < numElements * 8; ++evInd_itr)
  {
    int currentElementID = evInd_itr / 8;  // note the integer division
    veInds_data[evInds[evInd_itr]].push_back(currentElementID);
  }

  for(auto itr = veInds_data.begin(); itr != veInds_data.end(); itr++)
  {
    // Sort the vector
    std::sort(itr->second.begin(), itr->second.end());

    // Add the elements associated with the current vertex to veInds
    for(unsigned long i = 0; i < itr->second.size(); ++i)
      veInds.push_back(itr->second[i]);
  }

  // Generate the veBegins
  std::vector<mir::PosType> veBegins;
  veBegins.push_back(0);
  int currentIndexCount = 0;
  for(auto itr = veInds_data.begin(); itr != veInds_data.end(); itr++)
  {
    currentIndexCount += itr->second.size();
    veBegins.push_back(currentIndexCount);
  }

  // Generate the vertex positions
  std::vector<mir::Point2> points;
  for(int z = 0; z < gridSize + 1; ++z)
  {
    for(int y = gridSize; y > -1; --y)
    {
      for(int x = 0; x < gridSize + 1; ++x)
      {
        points.push_back(mir::Point2::make_point(x, y, z));
      }
    }
  }

  mir::CellData data;
  data.m_numVerts = numVertices;
  data.m_numElems = numElements;
  data.m_topology.m_evInds = evInds;
  data.m_topology.m_evBegins = evBegins;
  data.m_topology.m_veInds = veInds;
  data.m_topology.m_veBegins = veBegins;
  data.m_mapData.m_vertexPositions = points;

  return data;
}

//--------------------------------------------------------------------------------

}  // namespace mir
}  // namespace axom

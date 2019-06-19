// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "MeshTester.hpp"

namespace axom
{
namespace mir
{
      
//--------------------------------------------------------------------------------

MeshTester::MeshTester()
{

}

//--------------------------------------------------------------------------------

MeshTester::~MeshTester()
{

}

//--------------------------------------------------------------------------------

/// Initialize mesh for Test Case 1 (from Meredith 2004)
MIRMesh MeshTester::initTestCaseOne()
{
  int numElements = 9;
  int numVertices = 16;

  // Create the mesh connectivity information
  std::vector<mir::PosType>  evInds = {
      0,4,5,1,     // elem 0, card 4, start 0
      1,5,6,2,     // elem 1, card 4, start 4
      2,6,7,3,     // elem 2, card 4, start 8
      4,8,9,5,     // elem 3, card 4, start 12
      5,9,10,6,    // elem 4, card 4, start 16
      6,10,11,7,   // elem 5, card 4, start 20
      8,12,13,9,   // elem 6, card 4, start 24
      9,13,14,10,  // elem 7, card 4, start 28
      10,14,15,11  // elem 8, card 4, start 32, end 36
    };

  std::vector<mir::PosType>  evBegins = {
      0,4,8,12,16,20,24,28,32,36
    };
  std::vector<mir::PosType>  veInds = {
      0,          // vert  0, card 1, start 0
      0,1,        // vert  1, card 2, start 1
      1,2,        // vert  2, card 2, start 3
      2,          // vert  3, card 1, start 5
      0,3,        // vert  4, card 2, start 6
      0,1,3,4,    // vert  5, card 4, start 8
      1,2,4,5,    // vert  6, card 4, start 12
      2,5,        // vert  7, card 2, start 16
      3,6,        // vert  8, card 2, start 18
      3,4,6,7,    // vert  9, card 4, start 20
      4,5,7,8,    // vert  10, card 4, start 24
      5,8,        // vert  11, card 2, start 28
      6,          // vert  12, card 1, start 30
      6,7,        // vert  13, card 2, start 31
      7,8,        // vert  14, card 2, start 33
      8,          // vert  15, card 1, start 35, end 36
    };
  std::vector<mir::PosType>  veBegins = {
      0,1,3,5,6,8,12,16,18,20,24,28,30,31,33,35,36
    };

  mir::VertSet  verts = mir::VertSet(numVertices);  // Construct a vertex set with 16 vertices
  mir::ElemSet  elems = mir::ElemSet(numElements);   // Construct an element set with 9 elements


  int numMaterials = 2;
  enum { GREEN = 0, BLUE = 1 };

  std::vector<std::vector<axom::float64> > elementVF;   elementVF.resize(numMaterials);

  std::vector<axom::float64> greenVolumeFractions = {1.0, 1.0, 1.0, 1.0, 0.5, 0.2, 0.2, 0.0, 0.0};
  elementVF[GREEN] = greenVolumeFractions;
  std::vector<axom::float64> blueVolumeFractions = {0.0, 0.0, 0.0, 0.0, 0.5, 0.8, 0.8, 1.0, 1.0};
  elementVF[BLUE] = blueVolumeFractions;

  std::vector<mir::Point2> points =
  {
    mir::Point2( 0.0, 3.0 ),
    mir::Point2( 1.0, 3.0 ),
    mir::Point2( 2.0, 3.0 ),
    mir::Point2( 3.0, 3.0 ),

    mir::Point2( 0.0, 2.0 ),
    mir::Point2( 1.0, 2.0 ),
    mir::Point2( 2.0, 2.0 ),
    mir::Point2( 3.0, 2.0 ),

    mir::Point2( 0.0, 1.0 ),
    mir::Point2( 1.0, 1.0 ),
    mir::Point2( 2.0, 1.0 ),
    mir::Point2( 3.0, 1.0 ),

    mir::Point2( 0.0, 0.0 ),
    mir::Point2( 1.0, 0.0 ),
    mir::Point2( 2.0, 0.0 ),
    mir::Point2( 3.0, 0.0 )
  };

  CellTopologyData topology;
  topology.evInds = evInds;
  topology.evBegins = evBegins;
  topology.veInds = veInds;
  topology.veBegins = veBegins;

  CellMapData mapData;
  mapData.elementDominantMaterials = {NULL_MAT, NULL_MAT, NULL_MAT, NULL_MAT, NULL_MAT, NULL_MAT, NULL_MAT, NULL_MAT, NULL_MAT};
  mapData.elementParents = { 0,1,2,3,4,5,6,7,8 }; // For the base mesh, the parents are always themselves
  mapData.vertexPositions = points;

  // Build the mesh
  mir::MIRMesh testMesh;
  testMesh.initializeMesh(verts, elems, numMaterials, topology, mapData, elementVF);

  return testMesh;
}

//--------------------------------------------------------------------------------

/// Initialize mesh for Test Case 2 (from Meredith and Childs 2010)
mir::MIRMesh MeshTester::initTestCaseTwo()
{
  int numElements = 9;
  int numVertices = 16;

  // Create the mesh connectivity information
  std::vector<mir::PosType>  evInds = {
      0,4,5,1,     // elem 0, card 4, start 0
      1,5,6,2,     // elem 1, card 4, start 4
      2,6,7,3,     // elem 2, card 4, start 8
      4,8,9,5,     // elem 3, card 4, start 12
      5,9,10,6,    // elem 4, card 4, start 16
      6,10,11,7,   // elem 5, card 4, start 20
      8,12,13,9,   // elem 6, card 4, start 24
      9,13,14,10,  // elem 7, card 4, start 28
      10,14,15,11  // elem 8, card 4, start 32, end 36
    };

  std::vector<mir::PosType>  evBegins = {
      0,4,8,12,16,20,24,28,32,36
    };
  std::vector<mir::PosType>  veInds = {
      0,          // vert  0, card 1, start 0
      0,1,        // vert  1, card 2, start 1
      1,2,        // vert  2, card 2, start 3
      2,          // vert  3, card 1, start 5
      0,3,        // vert  4, card 2, start 6
      0,1,3,4,    // vert  5, card 4, start 8
      1,2,4,5,    // vert  6, card 4, start 12
      2,5,        // vert  7, card 2, start 16
      3,6,        // vert  8, card 2, start 18
      3,4,6,7,    // vert  9, card 4, start 20
      4,5,7,8,    // vert  10, card 4, start 24
      5,8,        // vert  11, card 2, start 28
      6,          // vert  12, card 1, start 30
      6,7,        // vert  13, card 2, start 31
      7,8,        // vert  14, card 2, start 33
      8,          // vert  15, card 1, start 35, end 36
    };
  std::vector<mir::PosType>  veBegins = {
      0,1,3,5,6,8,12,16,18,20,24,28,30,31,33,35,36
    };

  mir::VertSet  verts = mir::VertSet(numVertices);  // Construct a vertex set with 16 vertices
  mir::ElemSet  elems = mir::ElemSet(numElements);   // Construct an element set with 9 elements

  int numMaterials = 3;
  enum { BLUE = 0, RED = 1, ORANGE = 2 };

  std::vector<std::vector<axom::float64> > elementVF;   elementVF.resize(numMaterials);

  std::vector<axom::float64> blueVolumeFractions = {1.0, 1.0, 1.0, 1.0, 0.5, 0.2, 0.2, 0.0, 0.0};
  elementVF[BLUE] = blueVolumeFractions;
  std::vector<axom::float64>  redVolumeFractions = {0.0, 0.0, 0.0, 0.0, 0.3, 0.8, 0.0, 0.3, 1.0};
  elementVF[RED] = redVolumeFractions;
  std::vector<axom::float64> orangeVolumeFractions = {0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.8, 0.7, 0.0};
  elementVF[ORANGE] = orangeVolumeFractions;

  std::vector<mir::Point2> points =
  {
    mir::Point2( 0.0, 3.0 ),
    mir::Point2( 1.0, 3.0 ),
    mir::Point2( 2.0, 3.0 ),
    mir::Point2( 3.0, 3.0 ),

    mir::Point2( 0.0, 2.0 ),
    mir::Point2( 1.0, 2.0 ),
    mir::Point2( 2.0, 2.0 ),
    mir::Point2( 3.0, 2.0 ),

    mir::Point2( 0.0, 1.0 ),
    mir::Point2( 1.0, 1.0 ),
    mir::Point2( 2.0, 1.0 ),
    mir::Point2( 3.0, 1.0 ),

    mir::Point2( 0.0, 0.0 ),
    mir::Point2( 1.0, 0.0 ),
    mir::Point2( 2.0, 0.0 ),
    mir::Point2( 3.0, 0.0 )
  };

  CellTopologyData topology;
  topology.evInds = evInds;
  topology.evBegins = evBegins;
  topology.veInds = veInds;
  topology.veBegins = veBegins;

  CellMapData mapData;
  mapData.elementDominantMaterials = {NULL_MAT, NULL_MAT, NULL_MAT, NULL_MAT, NULL_MAT, NULL_MAT, NULL_MAT, NULL_MAT, NULL_MAT};
  mapData.elementParents = { 0,1,2,3,4,5,6,7,8 }; // For the base mesh, the parents are always themselves
  mapData.vertexPositions = points;

  // Build the mesh
  mir::MIRMesh testMesh;
  testMesh.initializeMesh(verts, elems, numMaterials, topology, mapData, elementVF);

  return testMesh;
}

//--------------------------------------------------------------------------------

/// Initialize mesh for Test Case 3, which tests all the triangle clipping cases.
mir::MIRMesh MeshTester::initTestCaseThree()
{
  int numElements = 4;
  int numVertices = 6;      // OR create a middle triangle with all of one material, and then a ring of triangles around it that are full of the other material

  // Create the mesh connectivity information
  std::vector<mir::PosType>  evInds = {
      0,1,2,     // elem 0, card 3, start 0
      1,3,4,     // elem 1, card 3, start 3
      1,4,2,     // elem 2, card 3, start 6
      2,4,5      // elem 3, card 3, start 9, end 12
    };

  std::vector<mir::PosType>  evBegins = {
      0,3,6,9,12
    };
  std::vector<mir::PosType>  veInds = {
      0,          // vert  0, card 1, start 0
      0,1,2,      // vert  1, card 3, start 1
      0,2,3,      // vert  2, card 3, start 4
      1,          // vert  3, card 1, start 7
      1,2,3,      // vert  4, card 3, start 8
      3           // vert  5, card 1, start 11, end 12
    };
  std::vector<mir::PosType>  veBegins = {
      0,1,4,7,8,11,12
    };

  mir::VertSet  verts = mir::VertSet(numVertices);  // Construct a vertex set with 24 vertices
  mir::ElemSet  elems = mir::ElemSet(numElements);   // Construct an element set with 8 elements

  int numMaterials = 2;
  enum { BLUE = 0, RED = 1, };

  std::vector<std::vector<axom::float64> > elementVF;   elementVF.resize(numMaterials);

  std::vector<axom::float64> blueVolumeFractions = {0.0, 0.5, 0.8, 0.5};
  elementVF[BLUE] = blueVolumeFractions;
  std::vector<axom::float64> redVolumeFractions = {1.0, 0.5, 0.2, 0.5};
  elementVF[RED] = redVolumeFractions;

  std::vector<mir::Point2> points =
  {
    mir::Point2( 1.0, 2.0 ),
    mir::Point2( 0.5, 1.0 ),
    mir::Point2( 1.5, 1.0 ),
    mir::Point2( 0.0, 0.0 ),
    mir::Point2( 1.0, 0.0 ),
    mir::Point2( 2.0, 0.0 )
  };


  CellTopologyData topology;
  topology.evInds = evInds;
  topology.evBegins = evBegins;
  topology.veInds = veInds;
  topology.veBegins = veBegins;

  CellMapData mapData;
  mapData.elementDominantMaterials = {NULL_MAT, NULL_MAT, NULL_MAT, NULL_MAT};
  mapData.elementParents = { 0,1,2,3 }; // For the base mesh, the parents are always themselves
  mapData.vertexPositions = points;

  // Build the mesh
  mir::MIRMesh testMesh;
  testMesh.initializeMesh(verts, elems, numMaterials, topology, mapData, elementVF);

  return testMesh;
}

//--------------------------------------------------------------------------------

/// Initialize mesh for Test Case 4, a 3x3 grid with a circle of one material in the middle
mir::MIRMesh MeshTester::initTestCaseFour()
{
  int numElements = 9;
  int numVertices = 16;

  // Create the mesh connectivity information
  std::vector<mir::PosType>  evInds = {
      0,4,5,1,     // elem 0, card 4, start 0
      1,5,6,2,     // elem 1, card 4, start 4
      2,6,7,3,     // elem 2, card 4, start 8
      4,8,9,5,     // elem 3, card 4, start 12
      5,9,10,6,    // elem 4, card 4, start 16
      6,10,11,7,   // elem 5, card 4, start 20
      8,12,13,9,   // elem 6, card 4, start 24
      9,13,14,10,  // elem 7, card 4, start 28
      10,14,15,11  // elem 8, card 4, start 32, end 36
    };

  std::vector<mir::PosType>  evBegins = {
      0,4,8,12,16,20,24,28,32,36
    };
  std::vector<mir::PosType>  veInds = {
      0,          // vert  0, card 1, start 0
      0,1,        // vert  1, card 2, start 1
      1,2,        // vert  2, card 2, start 3
      2,          // vert  3, card 1, start 5
      0,3,        // vert  4, card 2, start 6
      0,1,3,4,    // vert  5, card 4, start 8
      1,2,4,5,    // vert  6, card 4, start 12
      2,5,        // vert  7, card 2, start 16
      3,6,        // vert  8, card 2, start 18
      3,4,6,7,    // vert  9, card 4, start 20
      4,5,7,8,    // vert  10, card 4, start 24
      5,8,        // vert  11, card 2, start 28
      6,          // vert  12, card 1, start 30
      6,7,        // vert  13, card 2, start 31
      7,8,        // vert  14, card 2, start 33
      8,          // vert  15, card 1, start 35, end 36
    };
  std::vector<mir::PosType>  veBegins = {
      0,1,3,5,6,8,12,16,18,20,24,28,30,31,33,35,36
    };

  mir::VertSet  verts = mir::VertSet(numVertices);  // Construct a vertex set with 16 vertices
  mir::ElemSet  elems = mir::ElemSet(numElements);   // Construct an element set with 9 elements

  std::vector<mir::Point2> points =
  {
    mir::Point2( 0.0, 3.0 ),
    mir::Point2( 1.0, 3.0 ),
    mir::Point2( 2.0, 3.0 ),
    mir::Point2( 3.0, 3.0 ),

    mir::Point2( 0.0, 2.0 ),
    mir::Point2( 1.0, 2.0 ),
    mir::Point2( 2.0, 2.0 ),
    mir::Point2( 3.0, 2.0 ),

    mir::Point2( 0.0, 1.0 ),
    mir::Point2( 1.0, 1.0 ),
    mir::Point2( 2.0, 1.0 ),
    mir::Point2( 3.0, 1.0 ),

    mir::Point2( 0.0, 0.0 ),
    mir::Point2( 1.0, 0.0 ),
    mir::Point2( 2.0, 0.0 ),
    mir::Point2( 3.0, 0.0 )
  };

  int numMaterials = 2;
  enum { GREEN = 0, BLUE = 1 };

  std::vector<std::vector<axom::float64> > elementVF;   elementVF.resize(numMaterials);

  std::vector<axom::float64> greenVolumeFractions;      greenVolumeFractions.resize(numElements);
  std::vector<axom::float64> blueVolumeFractions;       blueVolumeFractions.resize(numElements);

  // Generate the element volume fractions for the circle
  mir::Point2 circleCenter(1.5, 1.5);
  axom::float64 circleRadius = 1.25;
  int gridSize = 1000;
  for (int i = 0; i < numElements; ++i)
  {
    greenVolumeFractions[i] = calculatePercentOverlapMonteCarlo(gridSize, circleCenter, circleRadius, points[evInds[i * 4 + 0]], points[evInds[i * 4 + 1]], points[evInds[i * 4 + 2]], points[evInds[i * 4 + 3]]);
    blueVolumeFractions[i] = 1.0 - greenVolumeFractions[i];
  }

  elementVF[GREEN] = greenVolumeFractions;
  elementVF[BLUE] = blueVolumeFractions;

  CellTopologyData topology;
  topology.evInds = evInds;
  topology.evBegins = evBegins;
  topology.veInds = veInds;
  topology.veBegins = veBegins;

  CellMapData mapData;
  mapData.elementDominantMaterials = {NULL_MAT, NULL_MAT, NULL_MAT, NULL_MAT, NULL_MAT, NULL_MAT, NULL_MAT, NULL_MAT, NULL_MAT};
  mapData.elementParents = { 0,1,2,3,4,5,6,7,8 }; // For the base mesh, the parents are always themselves
  mapData.vertexPositions = points;

  // Build the mesh
  mir::MIRMesh testMesh;
  testMesh.initializeMesh(verts, elems, numMaterials, topology, mapData, elementVF);

  return testMesh;
}

//--------------------------------------------------------------------------------

/// Intializes a uniform grid with a circle of one material surrounded by another material.
mir::MIRMesh MeshTester::createUniformGridTestCaseMesh(int gridSize, mir::Point2 circleCenter, axom::float64 circleRadius)
{
  // Generate the mesh topology
  mir::CellData cellData = generateGrid(gridSize);

  mir::VertSet  verts = mir::VertSet(cellData.numVerts);  // Construct the vertex set
  mir::ElemSet  elems = mir::ElemSet(cellData.numElems);   // Construct the element set

  int numMaterials = 2;
  enum { GREEN = 0, BLUE = 1 };

  std::vector<std::vector<axom::float64> > elementVF;   elementVF.resize(numMaterials);

  std::vector<axom::float64> greenVolumeFractions;      greenVolumeFractions.resize(cellData.numElems);
  std::vector<axom::float64> blueVolumeFractions;       blueVolumeFractions.resize(cellData.numElems);

  // Generate the element volume fractions for the circle
  int numMonteCarloSamples = 100;
  for (int i = 0; i < cellData.numElems; ++i)
  {
    greenVolumeFractions[i] = calculatePercentOverlapMonteCarlo(numMonteCarloSamples, circleCenter, circleRadius, 
                                                                cellData.mapData.vertexPositions[cellData.topology.evInds[i * 4 + 0]], 
                                                                cellData.mapData.vertexPositions[cellData.topology.evInds[i * 4 + 1]], 
                                                                cellData.mapData.vertexPositions[cellData.topology.evInds[i * 4 + 2]], 
                                                                cellData.mapData.vertexPositions[cellData.topology.evInds[i * 4 + 3]]);
    blueVolumeFractions[i] = 1.0 - greenVolumeFractions[i];
  }

  elementVF[GREEN] = greenVolumeFractions;
  elementVF[BLUE] = blueVolumeFractions;

  std::vector<int> elementParents;// For the base mesh, the parents are always themselves
  std::vector<int> elementDominantMaterials;
  for (int i = 0; i < cellData.numElems; ++i)
  {
    elementParents.push_back(i);
    elementDominantMaterials.push_back(NULL_MAT);
  }

  CellTopologyData topology;
  topology.evInds = cellData.topology.evInds;
  topology.evBegins = cellData.topology.evBegins;
  topology.veInds = cellData.topology.veInds;
  topology.veBegins = cellData.topology.veBegins;

  CellMapData mapData;
  mapData.elementDominantMaterials = elementDominantMaterials;
  mapData.elementParents = elementParents;
  mapData.vertexPositions = cellData.mapData.vertexPositions;

  // Build the mesh
  mir::MIRMesh testMesh;
  testMesh.initializeMesh(verts, elems, numMaterials, topology, mapData, elementVF);

  return testMesh;
}

//--------------------------------------------------------------------------------

// Assumes the quad vertices are ordered the same as the clipping case.
axom::float64 MeshTester::calculatePercentOverlapMonteCarlo(int gridSize, mir::Point2 circle_center, axom::float64 circle_radius, mir::Point2 q_p0, mir::Point2 q_p1, mir::Point2 q_p2, mir::Point2 q_p3)
{
  // Check if any of the quad's corners are within the circle
  axom::float64 distP0 = distance(q_p0, circle_center);
  axom::float64 distP1 = distance(q_p1, circle_center);
  axom::float64 distP2 = distance(q_p2, circle_center);
  axom::float64 distP3 = distance(q_p3, circle_center);

  if (distP0 < circle_radius && distP1 < circle_radius && distP2 < circle_radius && distP3 < circle_radius)
  {
    // The entire quad overlaps the circle
    return 1.0;
  }
  else if (distP0 < circle_radius || distP1 < circle_radius || distP2 < circle_radius || distP3 < circle_radius)
  {
    // Some of the quad overlaps the circle, so run the Monte Carlo sampling to determine how much
    axom::float64 delta_x = abs(q_p2.m_x - q_p1.m_x) / (double) (gridSize - 1);
    axom::float64 delta_y = abs(q_p0.m_y - q_p1.m_y) / (double) (gridSize - 1);
    int countOverlap = 0;
    for (int y = 0; y < gridSize; ++y)
    {
      for (int x = 0; x < gridSize; ++x)
      {
        mir::Point2 samplePoint(delta_x * x + q_p1.m_x, delta_y * y + q_p1.m_y);
        if (distance(samplePoint, circle_center) < circle_radius)
          ++countOverlap;
      }
    }
    return countOverlap / (double) (gridSize * gridSize);
  }
  else
  {
    // None of the quad overlaps the circle
    return 0;
  }
}

//--------------------------------------------------------------------------------

/// Generates a 2D uniform grid with n x n elements.
mir::CellData MeshTester::generateGrid(int n)
{
  // Generate the topology for a uniform quad mesh with n x n elements automatically
  int numElements = n * n;
  int numVertices = (n + 1) * (n + 1);

  // Generate the evInds
  std::vector<mir::PosType> evInds;
  for (int eID = 0; eID < numElements; ++eID)
  {
    int row = eID / n;  // note the integer division
    int vertsPerRow = n + 1;
    int elemsPerRow = n;

    evInds.push_back( (eID % elemsPerRow) + row * vertsPerRow + 0);
    evInds.push_back( (eID % elemsPerRow) + (row + 1) * vertsPerRow + 0);
    evInds.push_back( (eID % elemsPerRow) + (row + 1) * vertsPerRow + 1);
    evInds.push_back( (eID % elemsPerRow) + row * vertsPerRow + 1);
  }

  // Generate the evBegins
  std::vector<mir::PosType> evBegins;
  evBegins.push_back(0);
  for (int i = 0; i < n * n; ++i)
  {
    evBegins.push_back((i + 1) * 4);
  }

  // Generate the veInds
  std::map<int, std::vector<int> > veInds_data;
  std::vector<mir::PosType> veInds;
  for (int evInd_itr = 0; evInd_itr < numElements * 4; ++evInd_itr)
  {
    int currentElementID = evInd_itr / 4; // note the integer division
    veInds_data[evInds[evInd_itr]].push_back(currentElementID);
  }
  
  for (auto itr = veInds_data.begin(); itr != veInds_data.end(); itr++)
  {
    // Sort the vector
    std::sort(itr->second.begin(), itr->second.end());

    // Add the elements associated with the current vertex to veInds
    for (unsigned long i = 0; i < itr->second.size(); ++i)
      veInds.push_back(itr->second[i]);
  }

  // Generate the veBegins
  std::vector<mir::PosType> veBegins;
  veBegins.push_back(0);
  int currentIndexCount = 0;
  for (auto itr = veInds_data.begin(); itr != veInds_data.end(); itr++)
  {
    currentIndexCount += itr->second.size();
    veBegins.push_back(currentIndexCount);
  }

  // Generate the vertex positions
  std::vector<mir::Point2> points;
  for (int y = n; y > -1; --y)
  {
    for (int x = 0; x < n + 1; ++x)
    {
      points.push_back(mir::Point2(x, y));
    }
  }

  mir::CellData data;
  data.numVerts = numVertices;
  data.numElems = numElements;
  data.topology.evInds = evInds;
  data.topology.evBegins = evBegins;
  data.topology.veInds = veInds;
  data.topology.veBegins = veBegins;
  data.mapData.vertexPositions = points;

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
  //   printf("{%.2f, %.2f} ", points[i].m_x, points[i].m_y);
  // }
  // printf("}\n");

  return data;
}

//--------------------------------------------------------------------------------

/// Calculate the distance between the two given points.
axom::float64 MeshTester::distance(mir::Point2 p0, mir::Point2 p1)
{
  return sqrt( ((p1.m_x - p0.m_x) * (p1.m_x - p0.m_x)) + ((p1.m_y - p0.m_y) * (p1.m_y - p0.m_y)) );
}

//--------------------------------------------------------------------------------

/// Multiple materials, multiple concentric circles.
/// Note: Assumes each circle has a unique material.
mir::MIRMesh MeshTester::initTestCaseFive(int gridSize, int numCircles)
{

  // Generate the mesh topology
  mir::CellData cellData = generateGrid(gridSize);

  mir::VertSet  verts = mir::VertSet(cellData.numVerts);  // Construct the vertex set
  mir::ElemSet  elems = mir::ElemSet(cellData.numElems);   // Construct the element set

  // Generate the element volume fractions with concentric circles
  int numMaterials = numCircles + 1;
  int defaultMaterialID = numMaterials - 1;   // default material is always the last index

  mir::Point2 circleCenter(gridSize / 2.0, gridSize / 2.0);   // all circles are centered around the same point

  // Initialize the radii of the circles
  std::vector<axom::float64> circleRadii;   
  axom::float64 maxRadius = gridSize / 2.4;   // Note: The choice of divisor is arbitrary
  axom::float64 minRadius = gridSize / 8;     // Note: The choice of divisor is arbitrary

  axom::float64 radiusDelta;
  if (numCircles <= 1)
    radiusDelta = (maxRadius - minRadius);
  else
    radiusDelta = (maxRadius - minRadius) / (double) (numCircles - 1);

  for (int i = 0; i < numCircles; ++i)
  {
    circleRadii.push_back( minRadius + (i * radiusDelta) );
  }

  // Initialize all material volume fractions to 0
  std::vector<std::vector<axom::float64> > materialVolumeFractionsData;
  for (int i = 0; i < numMaterials; ++i)
  {
    std::vector<axom::float64> tempVec;
    tempVec.resize(cellData.numElems);
    materialVolumeFractionsData.push_back(tempVec);
  }

  // Use the uniform sampling method to generate volume fractions for each material
  for (int eID = 0; eID < cellData.numElems; ++eID)
  {
    mir::Point2 v0 = cellData.mapData.vertexPositions[cellData.topology.evInds[eID * 4 + 0]]; 
    mir::Point2 v1 = cellData.mapData.vertexPositions[cellData.topology.evInds[eID * 4 + 1]]; 
    mir::Point2 v2 = cellData.mapData.vertexPositions[cellData.topology.evInds[eID * 4 + 2]]; 
    mir::Point2 v3 = cellData.mapData.vertexPositions[cellData.topology.evInds[eID * 4 + 3]];

    // Run the uniform sampling to determine how much of the current cell is composed of each material
    int materialCount[numMaterials];  for (int i = 0; i < numMaterials; ++i) materialCount[i] = 0;

    for (int matID = 0; matID < numMaterials; ++matID)
    {
      materialVolumeFractionsData[matID][eID] = materialCount[matID] / (double) (gridSize * gridSize);
    }

    axom::float64 delta_x = abs(v2.m_x - v1.m_x) / (double) (gridSize - 1);
    axom::float64 delta_y = abs(v0.m_y - v1.m_y) / (double) (gridSize - 1);

    for (int y = 0; y < gridSize; ++y)
    {
      for (int x = 0; x < gridSize; ++x)
      {
        mir::Point2 samplePoint(delta_x * x + v1.m_x, delta_y * y + v1.m_y);
        bool isPointSampled = false;
        for (int cID = 0; cID < numCircles && !isPointSampled; ++cID)
        {
          if (distance(samplePoint, circleCenter) < circleRadii[cID])
          {
            materialCount[cID]++;
            isPointSampled = true;
          }
        }
        if (!isPointSampled)
        {
          // The point was not within any of the circles, so increment the count for the default material
          materialCount[defaultMaterialID]++;
        }
      }
    }

    // Assign the element volume fractions based on the count of the samples in each circle
    for (int matID = 0; matID < numMaterials; ++matID)
    {
      materialVolumeFractionsData[matID][eID] = materialCount[matID] / (double) (gridSize * gridSize);
    }
  }

  std::vector<int> elementParents; // For the base mesh, the parents are always themselves
  std::vector<int> elementDominantMaterials;
  for (int i = 0; i < cellData.numElems; ++i)
  {
    elementParents.push_back(i);
    elementDominantMaterials.push_back(NULL_MAT);
  }

  CellTopologyData topology;
  topology.evInds = cellData.topology.evInds;
  topology.evBegins = cellData.topology.evBegins;
  topology.veInds = cellData.topology.veInds;
  topology.veBegins = cellData.topology.veBegins;

  CellMapData mapData;
  mapData.elementDominantMaterials = elementDominantMaterials;
  mapData.elementParents = elementParents;
  mapData.vertexPositions = cellData.mapData.vertexPositions;

  // Build the mesh
  mir::MIRMesh testMesh;
  testMesh.initializeMesh(verts, elems, numMaterials, topology, mapData, materialVolumeFractionsData);

  return testMesh;
}

//--------------------------------------------------------------------------------

/// Calculate the number of corners of the quad that are within the circle
int MeshTester::circleQuadCornersOverlaps(mir::Point2 circle_center, axom::float64 circle_radius, mir::Point2 q_p0, mir::Point2 q_p1, mir::Point2 q_p2, mir::Point2 q_p3)
{
  // Check if any of the quad's corners are within the circle
  axom::float64 distP0 = distance(q_p0, circle_center);
  axom::float64 distP1 = distance(q_p1, circle_center);
  axom::float64 distP2 = distance(q_p2, circle_center);
  axom::float64 distP3 = distance(q_p3, circle_center);

  int numCorners = 0;

  if (distP0 < circle_radius)
    numCorners++;
  if (distP1 < circle_radius)
    numCorners++;
  if (distP2 < circle_radius)
    numCorners++;
  if (distP3 < circle_radius)
    numCorners++;

  return numCorners;
}

//--------------------------------------------------------------------------------

}
}

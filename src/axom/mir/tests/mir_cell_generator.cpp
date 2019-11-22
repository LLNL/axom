// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MIR_CELL_GENERATOR_TEST_H_
#define MIR_CELL_GENERATOR_TEST_H_

#include "gtest/gtest.h"

#include "axom/slic.hpp"
#include "axom/mir.hpp"

using namespace axom;

//----------------------------------------------------------------------

TEST(mir_cell_generator, generate_quad_topology)
{
  EXPECT_EQ(true, 1);
  // this function generates the evInds, evBegins, ... etc given the map of elements -> verts, and verts -> elements
  std::map<int, std::vector<int> > elementMap;
  elementMap[0] = { 4, 1, 5 };
  elementMap[1] = { 0, 4, 5, 2 };
  elementMap[2] = { 0, 2, 6, 7 };
  elementMap[3] = { 7, 6, 3 };
  
  std::map<int, std::vector<int> > vertexMap;
  vertexMap[0] = { 1, 2 };
  vertexMap[1] = { 0 };
  vertexMap[2] = { 1, 2 };
  vertexMap[3] = { 3 };
  vertexMap[4] = { 0, 1 };
  vertexMap[5] = { 0, 1 };
  vertexMap[6] = { 2, 3 };
  vertexMap[7] = { 2, 3 };

  mir::CellData cellData;
  mir::CellGenerator generator;
  generator.generateTopologyData(elementMap, vertexMap, cellData);

  EXPECT_EQ( 4, cellData.m_topology.m_evInds[0] );
  EXPECT_EQ( 1, cellData.m_topology.m_evInds[1] );
  EXPECT_EQ( 5, cellData.m_topology.m_evInds[2] );
  EXPECT_EQ( 0, cellData.m_topology.m_evInds[3] );
  EXPECT_EQ( 4, cellData.m_topology.m_evInds[4] );
  EXPECT_EQ( 5, cellData.m_topology.m_evInds[5] );
  EXPECT_EQ( 2, cellData.m_topology.m_evInds[6] );
  EXPECT_EQ( 0, cellData.m_topology.m_evInds[7] );
  EXPECT_EQ( 2, cellData.m_topology.m_evInds[8] );
  EXPECT_EQ( 6, cellData.m_topology.m_evInds[9] );
  EXPECT_EQ( 7, cellData.m_topology.m_evInds[10] );
  EXPECT_EQ( 7, cellData.m_topology.m_evInds[11] );
  EXPECT_EQ( 6, cellData.m_topology.m_evInds[12] );
  EXPECT_EQ( 3, cellData.m_topology.m_evInds[13] );

  EXPECT_EQ( 1, cellData.m_topology.m_veInds[0] );
  EXPECT_EQ( 2, cellData.m_topology.m_veInds[1] );
  EXPECT_EQ( 0, cellData.m_topology.m_veInds[2] );
  EXPECT_EQ( 1, cellData.m_topology.m_veInds[3] );
  EXPECT_EQ( 2, cellData.m_topology.m_veInds[4] );
  EXPECT_EQ( 3, cellData.m_topology.m_veInds[5] );
  EXPECT_EQ( 0, cellData.m_topology.m_veInds[6] );
  EXPECT_EQ( 1, cellData.m_topology.m_veInds[7] );
  EXPECT_EQ( 0, cellData.m_topology.m_veInds[8] );
  EXPECT_EQ( 1, cellData.m_topology.m_veInds[9] );
  EXPECT_EQ( 2, cellData.m_topology.m_veInds[10] );
  EXPECT_EQ( 3, cellData.m_topology.m_veInds[11] );
  EXPECT_EQ( 2, cellData.m_topology.m_veInds[12] );
  EXPECT_EQ( 3, cellData.m_topology.m_veInds[13] );

  EXPECT_EQ(  0, cellData.m_topology.m_evBegins[0] );
  EXPECT_EQ(  3, cellData.m_topology.m_evBegins[1] );
  EXPECT_EQ(  7, cellData.m_topology.m_evBegins[2] );
  EXPECT_EQ( 11, cellData.m_topology.m_evBegins[3] );

  EXPECT_EQ(  0, cellData.m_topology.m_veBegins[0] );
  EXPECT_EQ(  2, cellData.m_topology.m_veBegins[1] );
  EXPECT_EQ(  3, cellData.m_topology.m_veBegins[2] );
  EXPECT_EQ(  5, cellData.m_topology.m_veBegins[3] );
  EXPECT_EQ(  6, cellData.m_topology.m_veBegins[4] );
  EXPECT_EQ(  8, cellData.m_topology.m_veBegins[5] );
  EXPECT_EQ( 10, cellData.m_topology.m_veBegins[6] );
  EXPECT_EQ( 12, cellData.m_topology.m_veBegins[7] );
}

//----------------------------------------------------------------------

TEST(mir_cell_generator, generate_vertex_positions)
{
  mir::Shape shapeType = mir::Shape::Quad;

  std::map<int, std::vector<int> > vertexMap;
  vertexMap[0] = { 1, 2 };
  vertexMap[1] = { 0 };
  vertexMap[2] = { 1, 2 };
  vertexMap[3] = { 3 };
  vertexMap[4] = { 0, 1 };
  vertexMap[5] = { 0, 1 };
  vertexMap[6] = { 2, 3 };
  vertexMap[7] = { 2, 3 };

  std::vector<mir::Point2> originalVertexPositions;
  originalVertexPositions.push_back( mir::Point2::make_point(0.0, 1.0) );
  originalVertexPositions.push_back( mir::Point2::make_point(0.0, 0.0) );
  originalVertexPositions.push_back( mir::Point2::make_point(1.0, 0.0) );
  originalVertexPositions.push_back( mir::Point2::make_point(1.0, 1.0) );

  axom::float64 tValues[8] = { 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5 };

  mir::CellData cellData;
  mir::CellGenerator cellGenerator;
  cellGenerator.generateVertexPositions(shapeType, vertexMap, originalVertexPositions, tValues, cellData);

  const auto& positions = cellData.m_mapData.m_vertexPositions;
  EXPECT_NEAR(  positions[0][0], 0.0, 0.00001 );
  EXPECT_NEAR(  positions[0][1], 1.0, 0.00001 );

  EXPECT_NEAR(  positions[1][0], 0.0, 0.00001 );
  EXPECT_NEAR(  positions[1][1], 0.0, 0.00001 );

  EXPECT_NEAR(  positions[2][0], 1.0, 0.00001 );
  EXPECT_NEAR(  positions[2][1], 0.0, 0.00001 );

  EXPECT_NEAR(  positions[3][0], 1.0, 0.00001 );
  EXPECT_NEAR(  positions[3][1], 1.0, 0.00001 );

  EXPECT_NEAR(  positions[4][0], 0.0, 0.00001 );
  EXPECT_NEAR(  positions[4][1], 0.5, 0.00001 );

  EXPECT_NEAR(  positions[5][0], 0.5, 0.00001 );
  EXPECT_NEAR(  positions[5][1], 0.0, 0.00001 );

  EXPECT_NEAR(  positions[6][0], 1.0, 0.00001 );
  EXPECT_NEAR(  positions[6][1], 0.5, 0.00001 );

  EXPECT_NEAR(  positions[7][0], 0.5, 0.00001 );
  EXPECT_NEAR(  positions[7][1], 1.0, 0.00001 );
}

//----------------------------------------------------------------------

TEST(mir_cell_generator, generate_vertex_volume_fractions)
{
  mir::Shape shapeType = mir::Shape::Quad;

  std::map<int, std::vector<int> > vertexMap;
  vertexMap[0] = { 1, 2 };
  vertexMap[1] = { 0 };
  vertexMap[2] = { 1, 2 };
  vertexMap[3] = { 3 };
  vertexMap[4] = { 0, 1 };
  vertexMap[5] = { 0, 1 };
  vertexMap[6] = { 2, 3 };
  vertexMap[7] = { 2, 3 };

  std::vector<std::vector<axom::float64> > originalVertexVF(2);
  originalVertexVF[0] = { 0.0, 0.33, 0.67, 1.0 };
  originalVertexVF[1] = { 1.0, 0.67, 0.33, 0.0 };

  axom::float64 tValues[8] = { 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5 };

  mir::CellData cellData;
  mir::CellGenerator cellGenerator;
  cellGenerator.generateVertexVolumeFractions(shapeType, vertexMap, originalVertexVF, tValues, cellData);


  EXPECT_NEAR(  cellData.m_mapData.m_vertexVolumeFractions[0][0], 0.0, 0.00001 );
  EXPECT_NEAR(  cellData.m_mapData.m_vertexVolumeFractions[1][0], 1.0, 0.00001 );

  EXPECT_NEAR(  cellData.m_mapData.m_vertexVolumeFractions[0][1], 0.33, 0.00001 );
  EXPECT_NEAR(  cellData.m_mapData.m_vertexVolumeFractions[1][1], 0.67, 0.00001 );

  EXPECT_NEAR(  cellData.m_mapData.m_vertexVolumeFractions[0][2], 0.67, 0.00001 );
  EXPECT_NEAR(  cellData.m_mapData.m_vertexVolumeFractions[1][2], 0.33, 0.00001 );

  EXPECT_NEAR(  cellData.m_mapData.m_vertexVolumeFractions[0][3], 1.0, 0.00001 );
  EXPECT_NEAR(  cellData.m_mapData.m_vertexVolumeFractions[1][3], 0.0, 0.00001 );


  EXPECT_NEAR(  cellData.m_mapData.m_vertexVolumeFractions[0][4], 0.165, 0.00001 );
  EXPECT_NEAR(  cellData.m_mapData.m_vertexVolumeFractions[1][4], 0.835, 0.00001 );

  EXPECT_NEAR(  cellData.m_mapData.m_vertexVolumeFractions[0][5], 0.5, 0.00001 );
  EXPECT_NEAR(  cellData.m_mapData.m_vertexVolumeFractions[1][5], 0.5, 0.00001 );

  EXPECT_NEAR(  cellData.m_mapData.m_vertexVolumeFractions[0][6], 0.835, 0.00001 );
  EXPECT_NEAR(  cellData.m_mapData.m_vertexVolumeFractions[1][6], 0.165, 0.00001 );

  EXPECT_NEAR(  cellData.m_mapData.m_vertexVolumeFractions[0][7], 0.5, 0.00001 );
  EXPECT_NEAR(  cellData.m_mapData.m_vertexVolumeFractions[1][7], 0.5, 0.00001 );
}

//----------------------------------------------------------------------

TEST(mir_cell_generator, determine_clean_cell_material)
{
  mir::Shape shapeType = mir::Shape::Quad;

  std::vector<int> vertexIDs = { 0, 1, 5, 7 };

  int matOne = 0;
  int matTwo = 1;

  std::vector<std::vector<axom::float64> > originalVertexVF(2);
  originalVertexVF[0] = { 0.0, 0.33, 0.67, 1.0 };
  originalVertexVF[1] = { 1.0, 0.67, 0.33, 0.0 };

  mir::CellData cellData;
  mir::CellGenerator cellGenerator;

  int dominantMaterial = cellGenerator.determineCleanCellMaterial(shapeType, vertexIDs, matOne, matTwo, originalVertexVF);
  
  EXPECT_EQ( dominantMaterial, 1 );
}

//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::UnitTestLogger logger;  // create & initialize test logger,

  result = RUN_ALL_TESTS();
  return result;
}


#endif //  MIR_CELL_GENERATOR_TEST_H_

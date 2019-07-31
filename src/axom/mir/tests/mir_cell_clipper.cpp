// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MIR_CELL_CLIPPER_TEST_H_
#define MIR_CELL_CLIPPER_TEST_H_

#include "gtest/gtest.h"

#include "axom/slic.hpp"
#include "axom/mir.hpp"

#include <cmath>

using namespace axom;

//----------------------------------------------------------------------

TEST(mir_clipping_case, all_triangle_cases)
{
  mir::Shape shape = mir::Shape::Triangle;
  int numVerts = mir::utilities::numVerts(shape);
  int numCases = pow(2, numVerts);

  for (auto actualClippingCase = 0; actualClippingCase < numCases; ++actualClippingCase)
  {
    std::vector<axom::float64> matOneVF;
    std::vector<axom::float64> matTwoVF;

    for (auto bitIndex = 0; bitIndex < numVerts; ++bitIndex)
    {
      unsigned int shiftedBit = 1;
      shiftedBit = shiftedBit << (numVerts - 1 - bitIndex);

      axom::float64 matOneValue = 0.0;
      if (actualClippingCase & shiftedBit)
      {
        matOneValue = 1.0;
      }

      matOneVF.push_back( matOneValue );
      matTwoVF.push_back( 1.0 - matOneValue );
    }

    mir::CellClipper clipper;
    unsigned int computedClippingCase = clipper.determineClippingCase(shape, matOneVF, matTwoVF);

    EXPECT_EQ ( computedClippingCase, actualClippingCase );
  }
}

//----------------------------------------------------------------------

TEST(mir_clipping_case, all_quad_cases)
{
  mir::Shape shape = mir::Shape::Quad;
  int numVerts = mir::utilities::numVerts(shape);
  int numCases = pow(2, numVerts);

  for (auto actualClippingCase = 0; actualClippingCase < numCases; ++actualClippingCase)
  {
    std::vector<axom::float64> matOneVF;
    std::vector<axom::float64> matTwoVF;

    for (auto bitIndex = 0; bitIndex < numVerts; ++bitIndex)
    {
      unsigned int shiftedBit = 1;
      shiftedBit = shiftedBit << (numVerts - 1 - bitIndex);

      axom::float64 matOneValue = 0.0;
      if (actualClippingCase & shiftedBit)
      {
        matOneValue = 1.0;
      }

      matOneVF.push_back( matOneValue );
      matTwoVF.push_back( 1.0 - matOneValue );
    }

    mir::CellClipper clipper;
    unsigned int computedClippingCase = clipper.determineClippingCase(shape, matOneVF, matTwoVF);

    EXPECT_EQ ( computedClippingCase, actualClippingCase );
  }
}

//----------------------------------------------------------------------

TEST(mir_clipping_case, all_tetrahedron_cases)
{
  mir::Shape shape = mir::Shape::Tetrahedron;
  int numVerts = mir::utilities::numVerts(shape);
  int numCases = pow(2, numVerts);

  for (unsigned int actualClippingCase = 0; actualClippingCase < (unsigned int) numCases; ++actualClippingCase)
  {
    std::vector<axom::float64> matOneVF;
    std::vector<axom::float64> matTwoVF;

    for (unsigned int bitIndex = 0; bitIndex < (unsigned int) numVerts; ++bitIndex)
    {
      unsigned int shiftedBit = 1;
      shiftedBit = shiftedBit << (numVerts - 1 - bitIndex);

      axom::float64 matOneValue = 0.0;
      if (actualClippingCase & shiftedBit)
      {
        matOneValue = 1.0;
      }

      matOneVF.push_back( matOneValue );
      matTwoVF.push_back( 1.0 - matOneValue );
    }

    mir::CellClipper clipper;
    unsigned int computedClippingCase = clipper.determineClippingCase(shape, matOneVF, matTwoVF);

    EXPECT_EQ ( computedClippingCase, actualClippingCase );
  }
}

//----------------------------------------------------------------------

TEST(mir_clipping_case, all_pyramid_cases)
{
  mir::Shape shape = mir::Shape::Pyramid;
  int numVerts = mir::utilities::numVerts(shape);
  int numCases = pow(2, numVerts);

  for (unsigned int actualClippingCase = 0; actualClippingCase < (unsigned int) numCases; ++actualClippingCase)
  {
    std::vector<axom::float64> matOneVF;
    std::vector<axom::float64> matTwoVF;

    for (unsigned int bitIndex = 0; bitIndex < (unsigned int) numVerts; ++bitIndex)
    {
      unsigned int shiftedBit = 1;
      shiftedBit = shiftedBit << (numVerts - 1 - bitIndex);

      axom::float64 matOneValue = 0.0;
      if (actualClippingCase & shiftedBit)
      {
        matOneValue = 1.0;
      }

      matOneVF.push_back( matOneValue );
      matTwoVF.push_back( 1.0 - matOneValue );
    }

    mir::CellClipper clipper;
    unsigned int computedClippingCase = clipper.determineClippingCase(shape, matOneVF, matTwoVF);

    EXPECT_EQ ( computedClippingCase, actualClippingCase );
  }
}

//----------------------------------------------------------------------

TEST(mir_clipping_case, all_triangular_prism_cases)
{
  mir::Shape shape = mir::Shape::Triangular_Prism;
  int numVerts = mir::utilities::numVerts(shape);
  int numCases = pow(2, numVerts);

  for (unsigned int actualClippingCase = 0; actualClippingCase < (unsigned int) numCases; ++actualClippingCase)
  {
    std::vector<axom::float64> matOneVF;
    std::vector<axom::float64> matTwoVF;

    for (unsigned int bitIndex = 0; bitIndex < (unsigned int) numVerts; ++bitIndex)
    {
      unsigned int shiftedBit = 1;
      shiftedBit = shiftedBit << (numVerts - 1 - bitIndex);

      axom::float64 matOneValue = 0.0;
      if (actualClippingCase & shiftedBit)
      {
        matOneValue = 1.0;
      }

      matOneVF.push_back( matOneValue );
      matTwoVF.push_back( 1.0 - matOneValue );
    }

    mir::CellClipper clipper;
    unsigned int computedClippingCase = clipper.determineClippingCase(shape, matOneVF, matTwoVF);

    EXPECT_EQ ( computedClippingCase, actualClippingCase );
  }
}

//----------------------------------------------------------------------

TEST(mir_clipping_case, all_hexahedron_cases)
{
  mir::Shape shape = mir::Shape::Hexahedron;
  int numVerts = mir::utilities::numVerts(shape);
  int numCases = pow(2, numVerts);

  for (unsigned int actualClippingCase = 0; actualClippingCase < (unsigned int) numCases; ++actualClippingCase)
  {
    std::vector<axom::float64> matOneVF;
    std::vector<axom::float64> matTwoVF;

    for (unsigned int bitIndex = 0; bitIndex < (unsigned int) numVerts; ++bitIndex)
    {
      unsigned int shiftedBit = 1;
      shiftedBit = shiftedBit << (numVerts - 1 - bitIndex);

      axom::float64 matOneValue = 0.0;
      if (actualClippingCase & shiftedBit)
      {
        matOneValue = 1.0;
      }

      matOneVF.push_back( matOneValue );
      matTwoVF.push_back( 1.0 - matOneValue );
    }

    mir::CellClipper clipper;
    unsigned int computedClippingCase = clipper.determineClippingCase(shape, matOneVF, matTwoVF);

    EXPECT_EQ ( computedClippingCase, actualClippingCase );
  }
}

//----------------------------------------------------------------------

TEST(mir_clipping_place, clip_edge_when_mat_one_dominates)
{
  axom::float64 vfMatOneVertexOne = 0.0;
  axom::float64 vfMatTwoVertexOne = 1.0;

  axom::float64 vfMatOneVertexTwo = 0.0;
  axom::float64 vfMatTwoVertexTwo = 1.0;

  mir::CellClipper clipper;
  axom::float64 tValue = clipper.computeTValueOnEdge(vfMatOneVertexOne, vfMatTwoVertexOne, vfMatOneVertexTwo, vfMatTwoVertexTwo);

  EXPECT_DOUBLE_EQ(  tValue, 0.0  );
}

//----------------------------------------------------------------------

TEST(mir_clipping_place, clip_edge_when_mat_two_dominates)
{
  axom::float64 vfMatOneVertexOne = 1.0;
  axom::float64 vfMatTwoVertexOne = 0.0;

  axom::float64 vfMatOneVertexTwo = 1.0;
  axom::float64 vfMatTwoVertexTwo = 0.0;

  mir::CellClipper clipper;
  axom::float64 tValue = clipper.computeTValueOnEdge(vfMatOneVertexOne, vfMatTwoVertexOne, vfMatOneVertexTwo, vfMatTwoVertexTwo);

  EXPECT_DOUBLE_EQ(  tValue, 0.0  );
}

// //----------------------------------------------------------------------


TEST(mir_clipping_place, clip_edge_in_middle)
{
  axom::float64 vfMatOneVertexOne = 1.0;
  axom::float64 vfMatTwoVertexOne = 0.0;

  axom::float64 vfMatOneVertexTwo = 0.0;
  axom::float64 vfMatTwoVertexTwo = 1.0;

  mir::CellClipper clipper;
  axom::float64 tValue = clipper.computeTValueOnEdge(vfMatOneVertexOne, vfMatTwoVertexOne, vfMatOneVertexTwo, vfMatTwoVertexTwo);

  EXPECT_DOUBLE_EQ(  tValue, 0.5  );
}

// //----------------------------------------------------------------------

TEST(mir_clipping_place, clip_edge_with_one_null_material)
{
  axom::float64 vfMatOneVertexOne = -1.0;
  axom::float64 vfMatTwoVertexOne = 0.0;

  axom::float64 vfMatOneVertexTwo = -1.0;
  axom::float64 vfMatTwoVertexTwo = 0.0;

  mir::CellClipper clipper;
  axom::float64 tValue = clipper.computeTValueOnEdge(vfMatOneVertexOne, vfMatTwoVertexOne, vfMatOneVertexTwo, vfMatTwoVertexTwo);

  EXPECT_DOUBLE_EQ(  tValue, 0.0  ); 
}

// //----------------------------------------------------------------------

TEST(mir_clipping_place, clip_edge_with_two_null_material)
{
  axom::float64 vfMatOneVertexOne = -1.0;
  axom::float64 vfMatTwoVertexOne = -1.0;

  axom::float64 vfMatOneVertexTwo = -1.0;
  axom::float64 vfMatTwoVertexTwo = -1.0;

  mir::CellClipper clipper;
  axom::float64 tValue = clipper.computeTValueOnEdge(vfMatOneVertexOne, vfMatTwoVertexOne, vfMatOneVertexTwo, vfMatTwoVertexTwo);

  EXPECT_DOUBLE_EQ(  tValue, 0.0  ); 
}

//----------------------------------------------------------------------

TEST(mir_clipping_cell_and_vertex_output, clip_quad_case_zero)
{
  mir::Shape shape = mir::Shape::Quad;
  std::vector<axom::float64> matOneVF = {0.0, 0.0, 0.0, 0.0};
  std::vector<axom::float64> matTwoVF = {1.0, 1.0, 1.0, 1.0};

  std::vector<std::vector<axom::float64> > vertexVF;
  vertexVF.push_back(matOneVF);
  vertexVF.push_back(matTwoVF);

  std::map<int, std::vector<int> > newElements;
  std::map<int, std::vector<int> > newVertices;
  axom::float64 tValues[8] = { 0 };

  mir::CellClipper clipper;
  clipper.computeClippingPoints(shape, vertexVF, newElements, newVertices, tValues);
  
  EXPECT_EQ( 1, newElements.size() );
  EXPECT_EQ( 4, newVertices.size() );

  EXPECT_EQ( 0, newElements[0][0] );
  EXPECT_EQ( 1, newElements[0][1] );
  EXPECT_EQ( 2, newElements[0][2] );
  EXPECT_EQ( 3, newElements[0][3] );

  EXPECT_EQ( 0, newVertices[0][0] );
  EXPECT_EQ( 0, newVertices[1][0] );
  EXPECT_EQ( 0, newVertices[2][0] );
  EXPECT_EQ( 0, newVertices[3][0] );

  for (int i = 0; i < 8; ++i)
  {
    EXPECT_DOUBLE_EQ( 0, tValues[i] );
  }
}

//----------------------------------------------------------------------

TEST(mir_clipping_cell_and_vertex_output, clip_quad_case_one)
{
  mir::Shape shape = mir::Shape::Quad;
  std::vector<axom::float64> matOneVF = {0.0, 0.0, 0.0, 1.0};
  std::vector<axom::float64> matTwoVF = {1.0, 1.0, 1.0, 0.0};

  std::vector<std::vector<axom::float64> > vertexVF;
  vertexVF.push_back(matOneVF);
  vertexVF.push_back(matTwoVF);

  std::map<int, std::vector<int> > newElements;
  std::map<int, std::vector<int> > newVertices;
  axom::float64 tValues[8] = { 0 };

  mir::CellClipper clipper;
  clipper.computeClippingPoints(shape, vertexVF, newElements, newVertices, tValues);
  
  EXPECT_EQ( 3, newElements.size() );
  EXPECT_EQ( 6, newVertices.size() );

  // Check the first element
  EXPECT_EQ( 3, newElements[0][0] );
  EXPECT_EQ( 7, newElements[0][1] );
  EXPECT_EQ( 6, newElements[0][2] );

  // Check the second element
  EXPECT_EQ( 7, newElements[1][0] );
  EXPECT_EQ( 0, newElements[1][1] );
  EXPECT_EQ( 2, newElements[1][2] );
  EXPECT_EQ( 6, newElements[1][3] );

  // Check the third element
  EXPECT_EQ( 0, newElements[2][0] );
  EXPECT_EQ( 1, newElements[2][1] );
  EXPECT_EQ( 2, newElements[2][2] );

  // Check each vertex's associated elements
  EXPECT_EQ( 1, newVertices[0][0] );
  EXPECT_EQ( 2, newVertices[0][1] );

  EXPECT_EQ( 2, newVertices[1][0] );
  
  EXPECT_EQ( 1, newVertices[2][0] );
  EXPECT_EQ( 2, newVertices[2][1] );

  EXPECT_EQ( 0, newVertices[3][0] );
  
  EXPECT_EQ( 0, newVertices[6][0] );
  EXPECT_EQ( 1, newVertices[6][1] );

  EXPECT_EQ( 0, newVertices[7][0] );
  EXPECT_EQ( 1, newVertices[7][1] );
}

//----------------------------------------------------------------------

TEST(mir_clipping_cell_and_vertex_output, clip_quad_case_three)
{
  mir::Shape shape = mir::Shape::Quad;
  std::vector<axom::float64> matOneVF = {0.0, 0.0, 1.0, 1.0};
  std::vector<axom::float64> matTwoVF = {1.0, 1.0, 0.0, 0.0};

  std::vector<std::vector<axom::float64> > vertexVF;
  vertexVF.push_back(matOneVF);
  vertexVF.push_back(matTwoVF);

  std::map<int, std::vector<int> > newElements;
  std::map<int, std::vector<int> > newVertices;
  axom::float64 tValues[8] = { 0 };

  mir::CellClipper clipper;
  clipper.computeClippingPoints(shape, vertexVF, newElements, newVertices, tValues);
  
  EXPECT_EQ( 2, newElements.size() );
  EXPECT_EQ( 6, newVertices.size() );

  EXPECT_EQ( 2, newElements[0][0] );
  EXPECT_EQ( 3, newElements[0][1] );
  EXPECT_EQ( 7, newElements[0][2] );
  EXPECT_EQ( 5, newElements[0][3] );

  EXPECT_EQ( 5, newElements[1][0] );
  EXPECT_EQ( 7, newElements[1][1] );
  EXPECT_EQ( 0, newElements[1][2] );
  EXPECT_EQ( 1, newElements[1][3] );

  EXPECT_EQ( 1, newVertices[0][0] );
  EXPECT_EQ( 1, newVertices[1][0] );
  EXPECT_EQ( 0, newVertices[2][0] );
  EXPECT_EQ( 0, newVertices[3][0] );
  EXPECT_EQ( 0, newVertices[5][0] );
  EXPECT_EQ( 1, newVertices[5][1] );
  EXPECT_EQ( 0, newVertices[7][0] );
  EXPECT_EQ( 1, newVertices[7][1] );

}

//----------------------------------------------------------------------

TEST(mir_clipping_cell_and_vertex_output, clip_quad_case_five)
{
  mir::Shape shape = mir::Shape::Quad;
  std::vector<axom::float64> matOneVF = {0.0, 1.0, 0.0, 1.0};
  std::vector<axom::float64> matTwoVF = {1.0, 0.0, 1.0, 0.0};

  std::vector<std::vector<axom::float64> > vertexVF;
  vertexVF.push_back(matOneVF);
  vertexVF.push_back(matTwoVF);

  std::map<int, std::vector<int> > newElements;
  std::map<int, std::vector<int> > newVertices;
  axom::float64 tValues[8] = { 0 };

  mir::CellClipper clipper;
  clipper.computeClippingPoints(shape, vertexVF, newElements, newVertices, tValues);
  
  EXPECT_EQ( 4, newElements.size() );
  EXPECT_EQ( 8, newVertices.size() );

  EXPECT_EQ( 1, newElements[0][0] );
  EXPECT_EQ( 5, newElements[0][1] );
  EXPECT_EQ( 4, newElements[0][2] );
  
  EXPECT_EQ( 5, newElements[1][0] );
  EXPECT_EQ( 2, newElements[1][1] );
  EXPECT_EQ( 0, newElements[1][2] );
  EXPECT_EQ( 4, newElements[1][3] );

  EXPECT_EQ( 2, newElements[2][0] );
  EXPECT_EQ( 6, newElements[2][1] );
  EXPECT_EQ( 7, newElements[2][2] );
  EXPECT_EQ( 0, newElements[2][3] );

  EXPECT_EQ( 6, newElements[3][0] );
  EXPECT_EQ( 3, newElements[3][1] );
  EXPECT_EQ( 7, newElements[3][2] );

  EXPECT_EQ( 1, newVertices[0][0] );
  EXPECT_EQ( 2, newVertices[0][1] );
  EXPECT_EQ( 0, newVertices[1][0] );
  EXPECT_EQ( 1, newVertices[2][0] );
  EXPECT_EQ( 2, newVertices[2][1] );
  EXPECT_EQ( 3, newVertices[3][0] );
  EXPECT_EQ( 0, newVertices[4][0] );
  EXPECT_EQ( 1, newVertices[4][1] );
  EXPECT_EQ( 0, newVertices[5][0] );
  EXPECT_EQ( 1, newVertices[5][1] );
  EXPECT_EQ( 2, newVertices[6][0] );
  EXPECT_EQ( 3, newVertices[6][1] );
  EXPECT_EQ( 2, newVertices[7][0] );
  EXPECT_EQ( 3, newVertices[7][1] );

}

//----------------------------------------------------------------------

// TEST(clipping_table_mesh_generation, triangle_meshes)
// {
//   mir::Shape shape = mir::Shape::Triangle;
//   int numVerts = mir::utilities::numVerts(shape);

//   for (auto actualClippingCase = 0u; actualClippingCase < mir::triangleClipTableVec.size(); ++actualClippingCase)
//   {
//     // Initialize the mesh
//     int numElements = 1;
//     int numVertices = 3;

//     // Create the mesh connectivity information
//     mir::CellTopologyData topology;
//     topology.m_evInds = { 0,1,2 };
//     topology.m_evBegins = { 0,3 };
//     topology.m_veInds = { 0,0,0 };
//     topology.m_veBegins = { 0,1,2,3 };

//     mir::VertSet  verts = mir::VertSet(numVertices);
//     mir::ElemSet  elems = mir::ElemSet(numElements);

//     // Calculate the vertex volume fractions needed to clip with the current case
//     std::vector<axom::float64> matOneVF;
//     std::vector<axom::float64> matTwoVF;
//     std::string bitString("");
//     for (auto bitIndex = 0; bitIndex < numVerts; ++bitIndex)
//     {
//       unsigned int shiftedBit = 1;
//       shiftedBit = shiftedBit << (numVerts - 1 - bitIndex);

//       axom::float64 matOneValue;
//       if (actualClippingCase & shiftedBit)
//       {
//         matOneValue = 1.0;
//         bitString += "1";
//       }
//       else
//       {
//         matOneValue = 0.0;
//         bitString += "0";
//       }

//       matOneVF.push_back( matOneValue );
//       matTwoVF.push_back( 1.0 - matOneValue );
//     }

//     std::vector<std::vector<axom::float64> > vertexVF = {
//       matOneVF,
//       matTwoVF
//     };

//     std::vector<std::vector<axom::float64> > elementVF = {
//       { 0.5 },
//       { 0.5 }
//     };


//     std::vector<mir::Point2> points =
//     {
//       mir::Point2::make_point( 0.5, 0.717 ),
//       mir::Point2::make_point( 0.0, 0.0 ),
//       mir::Point2::make_point( 1.0, 0.0 )
//     };

//     mir::CellMapData mapData;
//     mapData.m_elementDominantMaterials = { mir::NULL_MAT };
//     mapData.m_elementParents = { 0 };
//     mapData.m_vertexPositions = points;
//     mapData.m_shapeTypes = { mir::Shape::Triangle };

//     // Build the mesh
//     mir::MIRMesh testMesh;
//     testMesh.initializeMesh(verts, elems, 2, topology, mapData, elementVF);
//     testMesh.constructMeshVolumeFractionsVertex(vertexVF);

//     // Clip the mesh using the actualClippingCase index
//     mir::MIRMesh outputMesh;
//     mir::InterfaceReconstructor interfaceReconstructor;
//     interfaceReconstructor.computeReconstructedInterface(testMesh, outputMesh);

//     // Write out the processed mesh
//     std::string dirName = std::string(AXOM_BIN_DIR) + "/meshes";
//     std::string fileName = "mir_clippingcase_triangle_" + std::to_string(actualClippingCase) + "_" + bitString + ".vtk";
//     outputMesh.writeMeshToFile(dirName, fileName, "/");

//     // Ensure that all elements have a positive volume
//     for (int eID = 0; eID < outputMesh.m_elems.size(); ++eID)
//     {
//       // Get the element's shape
//       mir::Shape shapeType = (mir::Shape) outputMesh.m_shapeTypes[eID];

//       // Get the element's points
//       std::vector<int> elementVertices;
//       for (int vID = 0; vID < outputMesh.m_bdry[eID].size(); ++vID)
//       {
//         elementVertices.push_back(outputMesh.m_bdry[eID][vID]);
//       }
//       std::vector<mir::Point2> points;
//       for (int vID = 0; vID < mir::utilities::numVerts(shapeType); ++vID)
//       {
//         int originalVID = elementVertices[vID];
//         points.push_back(outputMesh.m_vertexPositions[originalVID] );
//       }

//       axom::float64 volume = mir::utilities::computeShapeVolume( shapeType, points.data());

//       EXPECT_TRUE(  volume > 0.0  );
//     }
//   }    
// }

// //----------------------------------------------------------------------

// TEST(clipping_table_mesh_generation, quad_meshes)
// {
//   mir::Shape shape = mir::Shape::Quad;
//   int numVerts = mir::utilities::numVerts(shape);

//   for (unsigned int actualClippingCase = 0; actualClippingCase < mir::quadClipTableVec.size(); ++actualClippingCase)
//   {
//     // Initialize the mesh
//     int numElements = 1;
//     int numVertices = 4;

//     // Create the mesh connectivity information
//     mir::CellTopologyData topology;
//     topology.m_evInds = { 0,1,2,3 };
//     topology.m_evBegins = { 0,4 };
//     topology.m_veInds = { 0,0,0,0 };
//     topology.m_veBegins = { 0,1,2,3,4 };

//     mir::VertSet  verts = mir::VertSet(numVertices);
//     mir::ElemSet  elems = mir::ElemSet(numElements);

//     // Calculate the vertex volume fractions needed to clip with the current case
//     std::vector<axom::float64> matOneVF;
//     std::vector<axom::float64> matTwoVF;
//     std::string bitString("");
//     for (unsigned int bitIndex = 0; bitIndex < (unsigned int) numVerts; ++bitIndex)
//     {
//       unsigned int shiftedBit = 1;
//       shiftedBit = shiftedBit << (numVerts - 1 - bitIndex);

//       axom::float64 matOneValue;
//       if (actualClippingCase & shiftedBit)
//       {
//         matOneValue = 1.0;
//         bitString += "1";
//       }
//       else
//       {
//         matOneValue = 0.0;
//         bitString += "0";
//       }

//       matOneVF.push_back( matOneValue );
//       matTwoVF.push_back( 1.0 - matOneValue );
//     }

//     std::vector<std::vector<axom::float64> > vertexVF = {
//       matOneVF,
//       matTwoVF
//     };

//     std::vector<std::vector<axom::float64> > elementVF = {
//       { 0.5 },
//       { 0.5 }
//     };


//     std::vector<mir::Point2> points =
//     {
//       mir::Point2::make_point( 0.0, 1.0 ),
//       mir::Point2::make_point( 0.0, 0.0 ),
//       mir::Point2::make_point( 1.0, 0.0 ),
//       mir::Point2::make_point( 1.0, 1.0 )
//     };

//     mir::CellMapData mapData;
//     mapData.m_elementDominantMaterials = { mir::NULL_MAT };
//     mapData.m_elementParents = { 0 };
//     mapData.m_vertexPositions = points;
//     mapData.m_shapeTypes = { mir::Shape::Quad };

//     // Build the mesh
//     mir::MIRMesh testMesh;
//     testMesh.initializeMesh(verts, elems, 2, topology, mapData, elementVF);
//     testMesh.constructMeshVolumeFractionsVertex(vertexVF);

//     // Clip the mesh using the actualClippingCase index
//     mir::MIRMesh outputMesh;
//     mir::InterfaceReconstructor interfaceReconstructor;
//     interfaceReconstructor.computeReconstructedInterface(testMesh, outputMesh);

//     // Write out the processed mesh
//     std::string dirName = std::string(AXOM_BIN_DIR) + "/meshes";
//     std::string fileName = "mir_clippingcase_quad_" + std::to_string(actualClippingCase) + "_" + bitString + ".vtk";
//     outputMesh.writeMeshToFile(dirName, fileName, "/");

//     // Ensure that all elements have a positive volume
//     for (int eID = 0; eID < outputMesh.m_elems.size(); ++eID)
//     {
//       // Get the element's shape
//       mir::Shape shapeType = (mir::Shape) outputMesh.m_shapeTypes[eID];

//       // Get the element's points
//       std::vector<int> elementVertices;
//       for (int vID = 0; vID < outputMesh.m_bdry[eID].size(); ++vID)
//       {
//         elementVertices.push_back(outputMesh.m_bdry[eID][vID]);
//       }
//       std::vector<mir::Point2> points;
//       for (int vID = 0; vID < mir::utilities::numVerts(shapeType); ++vID)
//       {
//         int originalVID = elementVertices[vID];
//         points.push_back(outputMesh.m_vertexPositions[originalVID] );
//       }

//       axom::float64 volume = mir::utilities::computeShapeVolume( shapeType, points.data());

//       EXPECT_TRUE(  volume > 0.0  );
//     }
//   }    
// }

// //----------------------------------------------------------------------

// TEST(clipping_table_mesh_generation, tetrahedron_meshes)
// {
//   mir::Shape shape = mir::Shape::Tetrahedron;
//   int numVerts = mir::utilities::numVerts(shape);

//   for (unsigned int actualClippingCase = 0; actualClippingCase < mir::tetrahedronClipTableVec.size(); ++actualClippingCase)
//   {
//     // Initialize the mesh
//     int numElements = 1;
//     int numVertices = 4;

//     // Create the mesh connectivity information
//     mir::CellTopologyData topology;
//     topology.m_evInds = { 0,1,2,3 };
//     topology.m_evBegins = { 0,4 };
//     topology.m_veInds = { 0,0,0,0 };
//     topology.m_veBegins = { 0,1,2,3,4 };

//     mir::VertSet  verts = mir::VertSet(numVertices);
//     mir::ElemSet  elems = mir::ElemSet(numElements);

//     // Calculate the vertex volume fractions needed to clip with the current case
//     std::vector<axom::float64> matOneVF;
//     std::vector<axom::float64> matTwoVF;
//     std::string bitString("");
//     for (auto bitIndex = 0; bitIndex < numVerts; ++bitIndex)
//     {
//       unsigned int shiftedBit = 1;
//       shiftedBit = shiftedBit << (numVerts - 1 - bitIndex);

//       axom::float64 matOneValue;
//       if (actualClippingCase & shiftedBit)
//       {
//         matOneValue = 1.0;
//         bitString += "1";
//       }
//       else
//       {
//         matOneValue = 0.0;
//         bitString += "0";
//       }

//       matOneVF.push_back( matOneValue );
//       matTwoVF.push_back( 1.0 - matOneValue );
//     }

//     std::vector<std::vector<axom::float64> > vertexVF = {
//       matOneVF,
//       matTwoVF
//     };

//     std::vector<std::vector<axom::float64> > elementVF = {
//       { 0.5 },
//       { 0.5 }
//     };


//     std::vector<mir::Point2> points =
//     {
//       mir::Point2::make_point( 0.5, 0.717, 0.0 ),
//       mir::Point2::make_point( 0.0, 0.0, 0.0 ),
//       mir::Point2::make_point( 1.0, 0.0, 0.0 ),
//       mir::Point2::make_point( 0.5, 0.3585, 0.717 )
//     };

//     mir::CellMapData mapData;
//     mapData.m_elementDominantMaterials = { mir::NULL_MAT };
//     mapData.m_elementParents = { 0 };
//     mapData.m_vertexPositions = points;
//     mapData.m_shapeTypes = { mir::Shape::Tetrahedron };

//     // Build the mesh
//     mir::MIRMesh testMesh;
//     testMesh.initializeMesh(verts, elems, 2, topology, mapData, elementVF);
//     testMesh.constructMeshVolumeFractionsVertex(vertexVF);

//     // Clip the mesh using the actualClippingCase index
//     mir::MIRMesh outputMesh;
//     mir::InterfaceReconstructor interfaceReconstructor;
//     interfaceReconstructor.computeReconstructedInterface(testMesh, outputMesh);

//     // Write out the processed mesh
//     std::string dirName = std::string(AXOM_BIN_DIR) + "/meshes";
//     std::string fileName = "mir_clippingcase_tetrahedron_" + std::to_string(actualClippingCase) + "_" + bitString + ".vtk";
//     outputMesh.writeMeshToFile(dirName, fileName, "/");

//     // Ensure that all elements have a positive volume
//     axom::float64 totalVolume = 0.0;
//     for (int eID = 0; eID < outputMesh.m_elems.size(); ++eID)
//     {
//       // Get the element's shape
//       mir::Shape shapeType = (mir::Shape) outputMesh.m_shapeTypes[eID];

//       // Get the element's points
//       std::vector<int> elementVertices;
//       for (int vID = 0; vID < outputMesh.m_bdry[eID].size(); ++vID)
//       {
//         elementVertices.push_back(outputMesh.m_bdry[eID][vID]);
//       }
//       std::vector<mir::Point2> points;
//       for (int vID = 0; vID < mir::utilities::numVerts(shapeType); ++vID)
//       {
//         int originalVID = elementVertices[vID];
//         points.push_back(outputMesh.m_vertexPositions[originalVID] );
//       }

//       axom::float64 volume = mir::utilities::computeShapeVolume( shapeType, points.data() );

//       EXPECT_TRUE(  volume > 0.0  );

//       totalVolume += volume;
//     }

//     // Ensure that the total volume of the generated elements equals the original
//     EXPECT_DOUBLE_EQ( totalVolume, 0.0856815 );
//   }    
// }

//----------------------------------------------------------------------

TEST(clipping_table_mesh_generation, pyramid_meshes)
{
  mir::Shape shape = mir::Shape::Pyramid;
  int numVerts = mir::utilities::numVerts(shape);

  for (unsigned int actualClippingCase = 0; actualClippingCase < mir::pyramidClipTableVec.size(); ++actualClippingCase)
  {
    // if (actualClippingCase == 18 || actualClippingCase == 24 || actualClippingCase == 12 || actualClippingCase == 6 
    //    || actualClippingCase == 7 || actualClippingCase == 19 || actualClippingCase == 25 || actualClippingCase == 13 )
    if (actualClippingCase == 12)
    {
      // Initialize the mesh
      int numElements = 1;
      int numVertices = 5;

      // Create the mesh connectivity information
      mir::CellTopologyData topology;
      topology.m_evInds = { 0,1,2,3,4 };
      topology.m_evBegins = { 0,5 };
      topology.m_veInds = { 0,0,0,0,0 };
      topology.m_veBegins = { 0,1,2,3,4,5 };

      mir::VertSet  verts = mir::VertSet(numVertices);
      mir::ElemSet  elems = mir::ElemSet(numElements);

      // Calculate the vertex volume fractions needed to clip with the current case
      std::vector<axom::float64> matOneVF;
      std::vector<axom::float64> matTwoVF;
      std::string bitString("");
      for (unsigned int bitIndex = 0; bitIndex < (unsigned int) numVerts; ++bitIndex)
      {
        unsigned int shiftedBit = 1;
        shiftedBit = shiftedBit << (numVerts - 1 - bitIndex);

        axom::float64 matOneValue;
        if (actualClippingCase & shiftedBit)
        {
          matOneValue = 1.0;
          bitString += "1";
        }
        else
        {
          matOneValue = 0.0;
          bitString += "0";
        }

        matOneVF.push_back( matOneValue );
        matTwoVF.push_back( 1.0 - matOneValue );
      }

      std::vector<std::vector<axom::float64> > vertexVF = {
        matOneVF,
        matTwoVF
      };

      std::vector<std::vector<axom::float64> > elementVF = {
        { 0.5 },
        { 0.5 }
      };


      std::vector<mir::Point2> points =
      {
        mir::Point2::make_point( 0.0, 0.0, 0.0 ),
        mir::Point2::make_point( 1.0, 0.0, 0.0 ),
        mir::Point2::make_point( 1.0, 1.0, 0.0 ),
        mir::Point2::make_point( 0.0, 1.0, 0.0 ),
        mir::Point2::make_point( 0.5, 0.5, 0.717)
      };

      mir::CellMapData mapData;
      mapData.m_elementDominantMaterials = { mir::NULL_MAT };
      mapData.m_elementParents = { 0 };
      mapData.m_vertexPositions = points;
      mapData.m_shapeTypes = { mir::Shape::Pyramid };

      // Build the mesh
      mir::MIRMesh testMesh;
      testMesh.initializeMesh(verts, elems, 2, topology, mapData, elementVF);
      testMesh.constructMeshVolumeFractionsVertex(vertexVF);

      // Clip the mesh using the actualClippingCase index
      mir::MIRMesh outputMesh;
      mir::InterfaceReconstructor interfaceReconstructor;
      interfaceReconstructor.computeReconstructedInterface(testMesh, outputMesh);

      // Write out the processed mesh
      std::string dirName = std::string(AXOM_BIN_DIR) + "/meshes";
      std::string fileName = "mir_clippingcase_pyramid_" + std::to_string(actualClippingCase) + "_" + bitString + ".vtk";
      // outputMesh.print();
      outputMesh.writeMeshToFile(dirName, fileName, "/");

      // Ensure that all elements have a positive volume
      axom::float64 totalVolume = 0.0;
      for (int eID = 0; eID < outputMesh.m_elems.size(); ++eID)
      {
        // Get the element's shape
        mir::Shape shapeType = (mir::Shape) outputMesh.m_shapeTypes[eID];

        // Get the element's points
        std::vector<int> elementVertices;
        for (int vID = 0; vID < outputMesh.m_bdry[eID].size(); ++vID)
        {
          elementVertices.push_back(outputMesh.m_bdry[eID][vID]);
        }
        std::vector<mir::Point2> points;
        for (int vID = 0; vID < mir::utilities::numVerts(shapeType); ++vID)
        {
          int originalVID = elementVertices[vID];
          points.push_back(outputMesh.m_vertexPositions[originalVID] );
        }

        axom::float64 volume = mir::utilities::computeShapeVolume( shapeType, points.data());

        EXPECT_TRUE(  volume > 0.0  );

        totalVolume += volume;
      }
      
      // Ensure that the total volume of the generated elements equals the original
      EXPECT_NEAR( 0.239, totalVolume, 0.00001 );
    }  
  }
    
}

//----------------------------------------------------------------------

// TEST(clipping_table_mesh_generation, triangular_prism_meshes)
// {
//   mir::Shape shape = mir::Shape::Triangular_Prism;
//   int numVerts = mir::utilities::numVerts(shape);

//   for (unsigned int actualClippingCase = 0; actualClippingCase < mir::triangularPrismClipTableVec.size(); ++actualClippingCase)
//   {
//     // Initialize the mesh
//     int numElements = 1;
//     int numVertices = 6;

//     // Create the mesh connectivity information
//     mir::CellTopologyData topology;
//     topology.m_evInds = { 0,1,2,3,4,5 };
//     topology.m_evBegins = { 0,6 };
//     topology.m_veInds = { 0,0,0,0,0,0 };
//     topology.m_veBegins = { 0,1,2,3,4,5,6 };

//     mir::VertSet  verts = mir::VertSet(numVertices);
//     mir::ElemSet  elems = mir::ElemSet(numElements);

//     // Calculate the vertex volume fractions needed to clip with the current case
//     std::vector<axom::float64> matOneVF;
//     std::vector<axom::float64> matTwoVF;
//     std::string bitString("");
//     for (unsigned int bitIndex = 0; bitIndex < (unsigned int) numVerts; ++bitIndex)
//     {
//       unsigned int shiftedBit = 1;
//       shiftedBit = shiftedBit << (numVerts - 1 - bitIndex);

//       axom::float64 matOneValue;
//       if (actualClippingCase & shiftedBit)
//       {
//         matOneValue = 1.0;
//         bitString += "1";
//       }
//       else
//       {
//         matOneValue = 0.0;
//         bitString += "0";
//       }

//       matOneVF.push_back( matOneValue );
//       matTwoVF.push_back( 1.0 - matOneValue );
//     }

//     std::vector<std::vector<axom::float64> > vertexVF = {
//       matOneVF,
//       matTwoVF
//     };

//     std::vector<std::vector<axom::float64> > elementVF = {
//       { 0.5 },
//       { 0.5 }
//     };


//     std::vector<mir::Point2> points =
//     {
//       mir::Point2::make_point( 0.5, 0.717, 0.0 ),
//       mir::Point2::make_point( 0.0, 0.0, 0.0 ),
//       mir::Point2::make_point( 1.0, 0.0, 0.0 ),
//       mir::Point2::make_point( 0.5, 0.717, -2.0 ),
//       mir::Point2::make_point( 0.0, 0.0, -2.0 ),
//       mir::Point2::make_point( 1.0, 0.0, -2.0 )
//     };

//     mir::CellMapData mapData;
//     mapData.m_elementDominantMaterials = { mir::NULL_MAT };
//     mapData.m_elementParents = { 0 };
//     mapData.m_vertexPositions = points;
//     mapData.m_shapeTypes = { mir::Shape::Triangular_Prism };

//     // Build the mesh
//     mir::MIRMesh testMesh;
//     testMesh.initializeMesh(verts, elems, 2, topology, mapData, elementVF);
//     testMesh.constructMeshVolumeFractionsVertex(vertexVF);

//     // Clip the mesh using the actualClippingCase index
//     mir::MIRMesh outputMesh;
//     mir::InterfaceReconstructor interfaceReconstructor;
//     interfaceReconstructor.computeReconstructedInterface(testMesh, outputMesh);

//     // Write out the processed mesh
//     std::string dirName = std::string(AXOM_BIN_DIR) + "/meshes";
//     std::string fileName = "mir_clippingcase_triangular_prism_" + std::to_string(actualClippingCase) + "_" + bitString + ".vtk";
//     outputMesh.writeMeshToFile(dirName, fileName, "/");

//     // Ensure that all elements have a positive volume
//     axom::float64 totalVolume = 0.0;
//     for (int eID = 0; eID < outputMesh.m_elems.size(); ++eID)
//     {
//       // Get the element's shape
//       mir::Shape shapeType = (mir::Shape) outputMesh.m_shapeTypes[eID];

//       // Get the element's points
//       std::vector<int> elementVertices;
//       for (int vID = 0; vID < outputMesh.m_bdry[eID].size(); ++vID)
//       {
//         elementVertices.push_back(outputMesh.m_bdry[eID][vID]);
//       }
//       std::vector<mir::Point2> points;
//       for (int vID = 0; vID < mir::utilities::numVerts(shapeType); ++vID)
//       {
//         int originalVID = elementVertices[vID];
//         points.push_back(outputMesh.m_vertexPositions[originalVID] );
//       }

//       axom::float64 volume = mir::utilities::computeShapeVolume( shapeType, points.data());

//       EXPECT_TRUE(  volume >= 0.0  );

//       totalVolume += volume;
//     }

//     // Ensure that the total volume of the generated elements equals the original
//     EXPECT_NEAR( 0.717, totalVolume, 0.00001 );
//   }    
// }

// //----------------------------------------------------------------------

// TEST(clipping_table_mesh_generation, hexahedron_meshes)
// {
//   mir::Shape shape = mir::Shape::Hexahedron;
//   int numVerts = mir::utilities::numVerts(shape);

//   for (unsigned int actualClippingCase = 0; actualClippingCase < mir::hexahedronClipTableVec.size(); ++actualClippingCase)
//   {
//     // Initialize the mesh
//     int numElements = 1;
//     int numVertices = 8;

//     // Create the mesh connectivity information
//     mir::CellTopologyData topology;
//     topology.m_evInds = { 0,1,2,3,4,5,6,7 };
//     topology.m_evBegins = { 0,8 };
//     topology.m_veInds = { 0,0,0,0,0,0,0,0 };
//     topology.m_veBegins = { 0,1,2,3,4,5,6,7,8 };

//     mir::VertSet  verts = mir::VertSet(numVertices);
//     mir::ElemSet  elems = mir::ElemSet(numElements);

//     // Calculate the vertex volume fractions needed to clip with the current case
//     std::vector<axom::float64> matOneVF;
//     std::vector<axom::float64> matTwoVF;
//     std::string bitString("");
//     for (unsigned int bitIndex = 0; bitIndex < (unsigned int) numVerts; ++bitIndex)
//     {
//       unsigned int shiftedBit = 1;
//       shiftedBit = shiftedBit << (numVerts - 1 - bitIndex);

//       axom::float64 matOneValue;
//       if (actualClippingCase & shiftedBit)
//       {
//         matOneValue = 1.0;
//         bitString += "1";
//       }
//       else
//       {
//         matOneValue = 0.0;
//         bitString += "0";
//       }

//       matOneVF.push_back( matOneValue );
//       matTwoVF.push_back( 1.0 - matOneValue );
//     }

//     std::vector<std::vector<axom::float64> > vertexVF = {
//       matOneVF,
//       matTwoVF
//     };

//     std::vector<std::vector<axom::float64> > elementVF = {
//       { 0.5 },
//       { 0.5 }
//     };

//     std::vector<mir::Point2> points =
//     {
//       mir::Point2::make_point( 0.0, 0.0, 0.0 ),
//       mir::Point2::make_point( 1.0, 0.0, 0.0 ),
//       mir::Point2::make_point( 1.0, 0.0, -1.0 ),
//       mir::Point2::make_point( 0.0, 0.0, -1.0 ),
//       mir::Point2::make_point( 0.0, 1.0, 0.0 ),
//       mir::Point2::make_point( 1.0, 1.0, 0.0 ),
//       mir::Point2::make_point( 1.0, 1.0, -1.0 ),
//       mir::Point2::make_point( 0.0, 1.0, -1.0 )
//     };

//     mir::CellMapData mapData;
//     mapData.m_elementDominantMaterials = { mir::NULL_MAT };
//     mapData.m_elementParents = { 0 };
//     mapData.m_vertexPositions = points;
//     mapData.m_shapeTypes = { mir::Shape::Hexahedron };

//     // Build the mesh
//     mir::MIRMesh testMesh;
//     testMesh.initializeMesh(verts, elems, 2, topology, mapData, elementVF);
//     testMesh.constructMeshVolumeFractionsVertex(vertexVF);

//     // Clip the mesh using the actualClippingCase index
//     mir::MIRMesh outputMesh;
//     mir::InterfaceReconstructor interfaceReconstructor;
//     interfaceReconstructor.computeReconstructedInterface(testMesh, outputMesh);

//     // Write out the processed mesh
//     std::string dirName = std::string(AXOM_BIN_DIR) + "/meshes";
//     std::string fileName = "mir_clippingcase_hexahedron_" + std::to_string(actualClippingCase) + "_" + bitString + ".vtk";
//     outputMesh.writeMeshToFile(dirName, fileName, "/");

//     // Ensure that all elements have a positive volume
//     axom::float64 totalVolume = 0.0;
//     for (int eID = 0; eID < outputMesh.m_elems.size(); ++eID)
//     {
//       // Get the element's shape
//       mir::Shape shapeType = (mir::Shape) outputMesh.m_shapeTypes[eID];

//       // Get the element's points
//       std::vector<int> elementVertices;
//       for (int vID = 0; vID < outputMesh.m_bdry[eID].size(); ++vID)
//       {
//         elementVertices.push_back(outputMesh.m_bdry[eID][vID]);
//       }
//       std::vector<mir::Point2> points;
//       for (int vID = 0; vID < mir::utilities::numVerts(shapeType); ++vID)
//       {
//         int originalVID = elementVertices[vID];
//         points.push_back(outputMesh.m_vertexPositions[originalVID] );
//       }

//       axom::float64 volume = mir::utilities::computeShapeVolume( shapeType, points.data());

//       EXPECT_TRUE(  volume >= 0.0  );
      
//       totalVolume += volume;
//     }

//     // Ensure that the total volume of the generated elements equals the original
//     EXPECT_NEAR( 1.0, totalVolume, 0.00001 );
//   }    
// }

//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::UnitTestLogger logger;  // create & initialize test logger,

  result = RUN_ALL_TESTS();
  return result;
}


#endif //  MIR_CELL_CLIPPER_TEST_H_

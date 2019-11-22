// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MIR_CELL_CLIPPER_TEST_H_
#define MIR_CELL_CLIPPER_TEST_H_

#include "gtest/gtest.h"

#include "axom/slic.hpp"
#include "axom/mir.hpp"

using namespace axom;

//----------------------------------------------------------------------

TEST(mir_clipping_case, triangle_case_zero)
{
  mir::Shape shape = mir::Shape::Triangle;
  std::vector<axom::float64> matOneVF = {0.0, 0.0, 0.0};
  std::vector<axom::float64> matTwoVF = {1.0, 1.0, 1.0};

  mir::CellClipper clipper;
  unsigned int clippingCase = clipper.determineClippingCase(shape, matOneVF, matTwoVF);

  EXPECT_EQ( 0 , clippingCase);
}

//----------------------------------------------------------------------

TEST(mir_clipping_case, triangle_case_three)
{
  mir::Shape shape = mir::Shape::Triangle;
  std::vector<axom::float64> matOneVF = {0.0, 1.0, 1.0};
  std::vector<axom::float64> matTwoVF = {1.0, 0.0, 0.0};

  mir::CellClipper clipper;
  unsigned int clippingCase = clipper.determineClippingCase(shape, matOneVF, matTwoVF);
  EXPECT_EQ( 3 , clippingCase);
}

//----------------------------------------------------------------------

TEST(mir_clipping_case, triangle_case_four)
{
  mir::Shape shape = mir::Shape::Triangle;
  std::vector<axom::float64> matOneVF = {1.0, 0.0, 0.0};
  std::vector<axom::float64> matTwoVF = {0.0, 1.0, 1.0};

  mir::CellClipper clipper;
  unsigned int clippingCase = clipper.determineClippingCase(shape, matOneVF, matTwoVF);

  EXPECT_EQ( 4 , clippingCase);
}

//----------------------------------------------------------------------

TEST(mir_clipping_case, triangle_case_seven)
{
  mir::Shape shape = mir::Shape::Triangle;
  std::vector<axom::float64> matOneVF = {1.0, 1.0, 1.0};
  std::vector<axom::float64> matTwoVF = {0.0, 0.0, 0.0};

  mir::CellClipper clipper;
  unsigned int clippingCase = clipper.determineClippingCase(shape, matOneVF, matTwoVF);
  
  EXPECT_EQ( 7 , clippingCase);
}

//----------------------------------------------------------------------

TEST(mir_clipping_case, quad_case_zero)
{
  mir::Shape shape = mir::Shape::Quad;
  std::vector<axom::float64> matOneVF = {0.0, 0.0, 0.0, 0.0};
  std::vector<axom::float64> matTwoVF = {1.0, 1.0, 1.0, 1.0};

  mir::CellClipper clipper;
  unsigned int clippingCase = clipper.determineClippingCase(shape, matOneVF, matTwoVF);
  EXPECT_EQ( 0 , clippingCase);
}

//----------------------------------------------------------------------

TEST(mir_clipping_case, quad_case_one)
{
  mir::Shape shape = mir::Shape::Quad;
  std::vector<axom::float64> matOneVF = {0.0, 0.0, 0.0, 1.0};
  std::vector<axom::float64> matTwoVF = {1.0, 1.0, 1.0, 0.0};

  mir::CellClipper clipper;
  unsigned int clippingCase = clipper.determineClippingCase(shape, matOneVF, matTwoVF);
  EXPECT_EQ( 1 , clippingCase);
}

//----------------------------------------------------------------------

TEST(mir_clipping_case, quad_case_two)
{
  mir::Shape shape = mir::Shape::Quad;
  std::vector<axom::float64> matOneVF = {0.0, 0.0, 1.0, 0.0};
  std::vector<axom::float64> matTwoVF = {1.0, 1.0, 0.0, 1.0};

  mir::CellClipper clipper;
  unsigned int clippingCase = clipper.determineClippingCase(shape, matOneVF, matTwoVF);
  EXPECT_EQ( 2 , clippingCase);
}

//----------------------------------------------------------------------

TEST(mir_clipping_case, quad_case_three)
{
  mir::Shape shape = mir::Shape::Quad;
  std::vector<axom::float64> matOneVF = {0.0, 0.0, 1.0, 1.0};
  std::vector<axom::float64> matTwoVF = {1.0, 1.0, 0.0, 0.0};

  mir::CellClipper clipper;
  unsigned int clippingCase = clipper.determineClippingCase(shape, matOneVF, matTwoVF);
  EXPECT_EQ( 3 , clippingCase);
}

//----------------------------------------------------------------------

TEST(mir_clipping_case, quad_case_five)
{
  mir::Shape shape = mir::Shape::Quad;
  std::vector<axom::float64> matOneVF = {0.0, 1.0, 0.0, 1.0};
  std::vector<axom::float64> matTwoVF = {1.0, 0.0, 1.0, 0.0};

  mir::CellClipper clipper;
  unsigned int clippingCase = clipper.determineClippingCase(shape, matOneVF, matTwoVF);
  EXPECT_EQ( 5 , clippingCase);
}

//----------------------------------------------------------------------

TEST(mir_clipping_case, quad_case_ten)
{
  mir::Shape shape = mir::Shape::Quad;
  std::vector<axom::float64> matOneVF = {1.0, 0.0, 1.0, 0.0};
  std::vector<axom::float64> matTwoVF = {0.0, 1.0, 0.0, 1.0};

  mir::CellClipper clipper;
  unsigned int clippingCase = clipper.determineClippingCase(shape, matOneVF, matTwoVF);
  EXPECT_EQ( 10 , clippingCase);
}

//----------------------------------------------------------------------

TEST(mir_clipping_case, quad_case_fifteen)
{
  mir::Shape shape = mir::Shape::Quad;
  std::vector<axom::float64> matOneVF = {1.0, 1.0, 1.0, 1.0};
  std::vector<axom::float64> matTwoVF = {0.0, 0.0, 0.0, 0.0};

  mir::CellClipper clipper;
  unsigned int clippingCase = clipper.determineClippingCase(shape, matOneVF, matTwoVF);
  EXPECT_EQ( 15 , clippingCase);
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

  EXPECT_EQ( 0, newElements[0][0] );
  EXPECT_EQ( 1, newElements[0][1] );
  EXPECT_EQ( 5, newElements[0][2] );
  EXPECT_EQ( 7, newElements[0][3] );

  EXPECT_EQ( 7, newElements[1][0] );
  EXPECT_EQ( 5, newElements[1][1] );
  EXPECT_EQ( 2, newElements[1][2] );
  EXPECT_EQ( 3, newElements[1][3] );

  EXPECT_EQ( 0, newVertices[0][0] );
  EXPECT_EQ( 0, newVertices[1][0] );
  EXPECT_EQ( 1, newVertices[2][0] );
  EXPECT_EQ( 1, newVertices[3][0] );
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

  EXPECT_EQ( 4, newElements[0][0] );
  EXPECT_EQ( 1, newElements[0][1] );
  EXPECT_EQ( 5, newElements[0][2] );
  
  EXPECT_EQ( 0, newElements[1][0] );
  EXPECT_EQ( 4, newElements[1][1] );
  EXPECT_EQ( 5, newElements[1][2] );
  EXPECT_EQ( 2, newElements[1][3] );

  EXPECT_EQ( 0, newElements[2][0] );
  EXPECT_EQ( 2, newElements[2][1] );
  EXPECT_EQ( 6, newElements[2][2] );
  EXPECT_EQ( 7, newElements[2][3] );

  EXPECT_EQ( 7, newElements[3][0] );
  EXPECT_EQ( 6, newElements[3][1] );
  EXPECT_EQ( 3, newElements[3][2] );

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

int main(int argc, char* argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::UnitTestLogger logger;  // create & initialize test logger,

  result = RUN_ALL_TESTS();
  return result;
}


#endif //  MIR_CELL_CLIPPER_TEST_H_

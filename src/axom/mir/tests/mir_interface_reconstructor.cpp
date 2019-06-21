// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MIR_SMOKE_H_
#define MIR_SMOKE_H_

#include "gtest/gtest.h"

#include "axom/slic.hpp"
#include "axom/mir.hpp"

using namespace axom;

TEST(mir_quad_clipping, quad_clipping_case_zero)
{
  // Initialize a 3x3 mesh with 2 materials
  mir::MeshTester meshGenerator;
  mir::MIRMesh testMesh = meshGenerator.initQuadClippingTestMesh();

  // Initialize its volume fractions to custom values to guarantee this clipping case
  std::vector<std::vector<axom::float64> > vertexVF;
  vertexVF.resize(2);
  vertexVF[0] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  vertexVF[1] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
  testMesh.constructMeshVolumeFractionsVertex(vertexVF);

  int upperLeftVertexID = 5;
  int lowerLeftVertexID = 9;
  int lowerRightVertexID = 10;
  int upperRightVertexID = 6;
  
  mir::InterfaceReconstructor reconstructor;
  unsigned int clippingCase = reconstructor.determineQuadClippingCase( testMesh, 0, 1, upperLeftVertexID, lowerLeftVertexID, lowerRightVertexID, upperRightVertexID );

  EXPECT_EQ(  clippingCase,  0  );
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


#endif //  MIR_SMOKE_H_

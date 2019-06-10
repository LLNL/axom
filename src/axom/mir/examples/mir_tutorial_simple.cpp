// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core.hpp"  // for axom macros
// #include "axom/mir.hpp"  // for Mir classes & functions
#include "axom/slam.hpp"

#include "../MIRMesh.hpp"
#include "../InterfaceReconstructor.hpp"

// namespace aliases
//namespace mir     = axom::mir;
namespace numerics = axom::numerics;
namespace slam = axom::slam;
namespace mir = axom::mir;

//--------------------------------------------------------------------------------

/*!
 * \brief Tutorial main
 */
int main( int AXOM_NOT_USED(argc), char** AXOM_NOT_USED(argv) )
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

  std::vector<axom::float64*> materialVolumeFractionsData;
  materialVolumeFractionsData.resize(numMaterials);
  axom::float64 greenVolumeFractions[] =  {1.0, 1.0, 1.0, 1.0, 0.5, 0.2, 0.2, 0.0, 0.0};
  materialVolumeFractionsData[GREEN] = greenVolumeFractions;
  axom::float64 blueVolumeFractions[] = {0.0, 0.0, 0.0, 0.0, 0.5, 0.8, 0.8, 1.0, 1.0};
  materialVolumeFractionsData[BLUE] = blueVolumeFractions;


  mir::Point2 points[numVertices];
  {
    points[0] = mir::Point2( 0.0, 3.0 );
    points[1] = mir::Point2( 1.0, 3.0 );
    points[2] = mir::Point2( 2.0, 3.0 );
    points[3] = mir::Point2( 3.0, 3.0 );

    points[4] = mir::Point2( 0.0, 2.0 );
    points[5] = mir::Point2( 1.0, 2.0 );
    points[6] = mir::Point2( 2.0, 2.0 );
    points[7] = mir::Point2( 3.0, 2.0 );

    points[8] = mir::Point2( 0.0, 1.0 );
    points[9] = mir::Point2( 1.0, 1.0 );
    points[10] = mir::Point2( 2.0, 1.0 );
    points[11] = mir::Point2( 3.0, 1.0 );

    points[12] = mir::Point2( 0.0, 0.0 );
    points[13] = mir::Point2( 1.0, 0.0 );
    points[14] = mir::Point2( 2.0, 0.0 );
    points[15] = mir::Point2( 3.0, 0.0 );
  }

  int elementParents[9] = { 0,1,2,3,4,5,6,7,8 }; // For the base mesh, the parents are always themselves

  std::vector<int> elementDominantMaterials = {NULL_MAT, NULL_MAT, NULL_MAT, NULL_MAT, NULL_MAT, NULL_MAT, NULL_MAT, NULL_MAT, NULL_MAT};

  // Build the mesh
  mir::MIRMesh testMesh;
  testMesh.InitializeMesh(evInds, evBegins, veInds, veBegins, verts, elems, numMaterials);
  testMesh.constructMeshRelations();
  testMesh.constructMeshVolumeFractionMaps(materialVolumeFractionsData);
  testMesh.constructVertexPositionMap(points);
  testMesh.constructElementParentMap(elementParents);
  testMesh.constructElementDominantMaterialMap(elementDominantMaterials);

  // Begin material interface reconstruction
  mir::InterfaceReconstructor reconstructor(&testMesh);
  mir::MIRMesh processedMesh = reconstructor.computeReconstructedInterface();
  processedMesh.print();
  processedMesh.writeMeshToFile("/Users/sterbentz3/Desktop/processedTestMesh.vtk");

  return 0;
}

//--------------------------------------------------------------------------------


/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */



#include "gtest/gtest.h"

#include "common/CommonTypes.hpp"
#include "common/Timer.hpp"

#include "quest/InOutOctree.hpp"
#include "quest/Orientation.hpp"
#include "quest/Triangle.hpp"

#include "mint/Mesh.hpp"

#include "slic/slic.hpp"

#include "quest_test_utilities.hpp"


namespace {
    const int NUM_PT_TESTS = 50000;
    const int DIM = 3;
}

typedef quest::InOutOctree<DIM> Octree3D;

typedef Octree3D::GeometricBoundingBox GeometricBoundingBox;
typedef Octree3D::SpacePt SpacePt;
typedef Octree3D::SpaceVector SpaceVector;
typedef Octree3D::GridPt GridPt;
typedef Octree3D::BlockIndex BlockIndex;


#include <cstdlib>
#include <limits>

// Uncomment the line below for true randomized points
#ifndef INOUT_OCTREE_TESTER_SHOULD_SEED
//  #define INOUT_OCTREE_TESTER_SHOULD_SEED
#endif

#ifdef INOUT_OCTREE_TESTER_SHOULD_SEED
  #include <ctime>      // for time() used by srand()
#endif




void queryOctahedronMesh(mint::Mesh*& mesh, const GeometricBoundingBox& bbox)
{
    const double bbMin = bbox.getMin()[0];
    const double bbMax = bbox.getMax()[0];

    Octree3D octree(bbox, mesh);
    octree.generateIndex();

    // Query the mesh containment
    asctoolkit::utilities::Timer timer(true);
    for(int i=0; i < NUM_PT_TESTS; ++i)
    {
        SpacePt pt;

//        switch(i)  // test a few special points and a lot of random points
//        {
//        case 0: case 1:  case 2:            // Test point at mesh vertices
//        case 3: case 4:  case 5:
//            pt = verts[i];
//            break;
//        case 6:  case 7:  case 8:  case 9:  // Test point at triangle centers
//        case 10: case 11: case 12: case 13:
//        {
//            int tIdx = (i-6)*VERTS_PER_TRI;
//            pt = getCentroid( verts[ tvRelation[tIdx]]
//                            , verts[ tvRelation[tIdx+1]]
//                            , verts[ tvRelation[tIdx+2]] );
//        }
//            break;
//        case 14: case 15: case 16: case 17:   // Test point at edge centers
//        case 18: case 19: case 20: case 21:
//        case 22: case 23: case 24: case 25:
//        {
//            int eIdx = (i-14)*VERTS_PER_EDGE;
//            pt = getCentroid( verts[ evRelation[eIdx]]
//                            , verts[ evRelation[eIdx +1]] );
//        }
//            break;
//        case 26:                 // origin
//            pt = SpacePt();
//            break;
//        case 27:                 // outside bounding box
//            pt = SpacePt(2* bbMin);
//            break;
//        case 28:                 // outside bounding box
//            pt = SpacePt(2* bbMax);
//            break;
//        default:                // random points in bounding box
            pt = quest::utilities::randomSpacePt<DIM>(bbMin, bbMax);
//            break;
//        }

        double absCoordSum = std::abs(pt[0]) + std::abs(pt[1]) + std::abs(pt[2]);

        // For the time being, we allow the within() test to fail when the
        // query point is sufficiently close to the surface
        bool expectInside = absCoordSum < 1.;
        EXPECT_TRUE( octree.within(pt) == expectInside
                || asctoolkit::utilities::isNearlyEqual(absCoordSum, 1.) )
            << "Point " << pt << " was not "
            << (expectInside? "inside" : "outside")
            << " surface of octahedron as expected."
            << " Sum of absolute values of coords was: " << absCoordSum
            << " and point is inside when this is less than 1.";

    }
    timer.stop();

    SLIC_INFO("-- querying octahedron with " << NUM_PT_TESTS
              << " points took " << timer.elapsed() << " seconds.");
    SLIC_INFO("***");
}



TEST( quest_inout_octree, octahedron_mesh)
{
    SLIC_INFO("*** This test creates a simple mesh of an octahedron and tests point containment.\n");

    // Generate the InOutOctree
    mint::Mesh* mesh = quest::utilities::make_octahedron_mesh();
    // quest::utilities::write_vtk(mesh, "octahedron.vtk");

    ///
    SpacePt ptNeg1(-1.);
    SpacePt ptPos1( 1.);
    GeometricBoundingBox bbox1(ptNeg1, ptPos1);
    SLIC_INFO("Testing InOutOctree on octahedron mesh with bounding box " << bbox1);
    queryOctahedronMesh(mesh, bbox1);

    ///
    SpacePt ptNeg2(-2.);
    SpacePt ptPos2( 2.);
    GeometricBoundingBox bbox2(ptNeg2, ptPos2);
    SLIC_INFO("Testing InOutOctree on octahedron mesh with bounding box " << bbox2);
    queryOctahedronMesh(mesh, bbox2);

    bbox2.shift( SpaceVector(0.01));
    SLIC_INFO("Testing InOutOctree on octahedron mesh with shifted bounding box " << bbox2);
    queryOctahedronMesh(mesh, bbox2);


    delete mesh;
    mesh = ATK_NULLPTR;
}


//----------------------------------------------------------------------
//----------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using asctoolkit::slic::UnitTestLogger;

int main(int argc, char * argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,

  asctoolkit::slic::setLoggingMsgLevel( asctoolkit::slic::message::Debug);


  // finalized when exiting main scope

#ifdef INOUT_OCTREE_TESTER_SHOULD_SEED
  std::srand( std::time(0) );
#else
  std::srand( 105);
#endif


  result = RUN_ALL_TESTS();

  return result;
}


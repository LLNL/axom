/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#include "gtest/gtest.h"

#include "axom/quest/geom/SpatialOctree.hpp"

#include "axom/slic/interface/slic.hpp"


TEST( quest_spatial_octree, spatial_octree_point_location)
{
  SLIC_INFO("*** This test verifies that a query point falls into "
            << " a child block.");

  static const int DIM = 3;
  typedef axom::quest::BlockData LeafNodeType;

  typedef axom::quest::SpatialOctree<DIM, LeafNodeType> OctreeType;
  typedef OctreeType::BlockIndex BlockIndex;
  typedef OctreeType::SpacePt SpacePt;
  typedef OctreeType::GeometricBoundingBox GeometricBoundingBox;


  GeometricBoundingBox bb(SpacePt(10), SpacePt(20));

  // Generate a point within the bounding box
  double alpha = 2./3.;
  SpacePt queryPt = SpacePt::lerp(bb.getMin(), bb.getMax(), alpha);
  EXPECT_TRUE( bb.contains(queryPt));

  OctreeType octree(bb);

  // Check that the point lies in a leaf of the tree
  // and that this is the root of the tree
  BlockIndex leafBlock = octree.findLeafBlock(queryPt);
  EXPECT_TRUE( octree.isLeaf(leafBlock));
  EXPECT_EQ( octree.root(), leafBlock );

  GeometricBoundingBox leafBB = octree.blockBoundingBox(leafBlock);
  EXPECT_TRUE( leafBB.contains( queryPt ));
  EXPECT_TRUE( bb.contains(leafBB));

  SLIC_INFO(
    "Query pt: " << queryPt
                 <<"\n\t" << ( leafBB.contains(queryPt) ? " was" : " was NOT" )
                 <<" contained in bounding box " << leafBB
                 <<"\n\t of octree root " << leafBlock   );

  for(int i=0 ; i< octree.maxInternalLevel() ; ++i)
  {
    octree.refineLeaf( leafBlock );
    leafBlock = octree.findLeafBlock(queryPt);
    EXPECT_TRUE( octree.isLeaf(leafBlock));

    leafBB = octree.blockBoundingBox(leafBlock);
    EXPECT_TRUE( leafBB.contains( queryPt ));
    EXPECT_TRUE( bb.contains(leafBB));

    SLIC_INFO(
      "Query pt: "
      << queryPt
      <<"\n\t" << ( leafBB.contains(queryPt) ? " was" : " was not")
      <<" contained in bounding box " << leafBB
      <<"\n\t of leaf " << leafBlock
      <<" in the octree. ");
  }


}


//----------------------------------------------------------------------
//----------------------------------------------------------------------
#include "axom/slic/core/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}

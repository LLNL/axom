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

#include "quest/Point.hpp"
#include "quest/SpatialOctree.hpp"

//------------------------------------------------------------------------------
TEST( quest_octree, topological_octree_parent_child)
{
  SLIC_INFO("*** This test exercises the parent/child relation in quest::OctreeBase");

  static const int DIM = 2;
  typedef int CoordType;
  typedef quest::LeafData LeafNodeType;

  typedef quest::OctreeBase<DIM, LeafNodeType> OctreeType;
  typedef OctreeType::GridPt GridPt;
  typedef OctreeType::BlockIndex BlockIndex;
  typedef BlockIndex::ChildIndexSet OctreeChildIndexSet;

  OctreeType octree;

  BlockIndex rootBlock = octree.root();
  EXPECT_EQ( rootBlock.pt(), GridPt());
  EXPECT_EQ( rootBlock.level(), 0);

  OctreeChildIndexSet childIndexSet;
  int numChildren  = childIndexSet.size();
  for(int i = 0; i < numChildren; ++i)
  {
      // Root's children are at level one and have coordinates that are either zero or one
      BlockIndex childBlock = octree.child( rootBlock, childIndexSet[i]);
      EXPECT_EQ( childBlock.level(), 1);
      EXPECT_EQ( childBlock.level(), rootBlock.childLevel());

      int recombineIndex = 0;
      for(int dim = 0; dim < DIM; ++dim)
      {
          CoordType coordVal = childBlock.pt()[dim];
          EXPECT_TRUE(  coordVal == 0 || coordVal == 1);

          bool expBit = (childIndexSet[i] & 1<<dim);
          EXPECT_EQ(  coordVal, expBit? 1 : 0);

          recombineIndex += (coordVal<<dim);
      }
      EXPECT_EQ(childIndexSet[i], recombineIndex);

      // Parent of child is self
      EXPECT_EQ( octree.parent( childBlock), rootBlock);
  }

}

//------------------------------------------------------------------------------
TEST( quest_octree, topological_octree_refine)
{
  SLIC_INFO("*** This test exercises the block refinement in quest::OctreeBase"
              <<"\nSpecifically, that refining the root block adds all its children to the octree.");

  static const int DIM = 3;
  typedef quest::LeafData LeafNodeType;

  typedef quest::OctreeBase<DIM, LeafNodeType> OctreeType;
  typedef OctreeType::BlockIndex BlockIndex;

  OctreeType octree;

  BlockIndex rootBlock = octree.root();
  octree.refineLeaf( rootBlock );

  int numChildren  = BlockIndex::numChildren();
  for(int i = 0; i < numChildren; ++i)
  {
      EXPECT_TRUE( octree.isLeaf( octree.child(rootBlock, i)) );
  }

}


TEST( quest_octree, spatial_octree_point_location)
{
    SLIC_INFO("*** This test verifies that a query point falls into "
            << " a child block.");

    static const int DIM = 3;
    typedef quest::LeafData  LeafNodeType;

    typedef quest::SpatialOctree<DIM, LeafNodeType> OctreeType;
    typedef OctreeType::BlockIndex BlockIndex;
    typedef OctreeType::SpacePt SpacePt;
    typedef OctreeType::GeometricBoundingBox GeometricBoundingBox;


    GeometricBoundingBox bb(SpacePt(10), SpacePt(20));

    // Generate a point within the bounding box
    double alpha = 2./3.;
    SpacePt queryPt = SpacePt::lerp(bb.getMin(), bb.getMax(), alpha);
    EXPECT_TRUE( bb.contains(queryPt));

    OctreeType octree(bb);

    // Check that the point lies in a leaf of the tree and that this is the root of the tree
    BlockIndex leafBlock = octree.findLeafBlock(queryPt);
    EXPECT_TRUE( octree.isLeaf(leafBlock));
    EXPECT_EQ( octree.root(), leafBlock );

    GeometricBoundingBox leafBB = octree.blockBoundingBox(leafBlock);
    EXPECT_TRUE( leafBB.contains( queryPt ));
    EXPECT_TRUE( bb.contains(leafBB));

    SLIC_INFO("Query pt: " << queryPt
              <<"\n\t" << ( leafBB.contains(queryPt) ? " was" : " was NOT" )
              <<" contained in bounding box " << leafBB
              <<"\n\t of octree root " << leafBlock
              );

    for(int i=0; i< octree.maxInternalLevel(); ++i)
    {
        octree.refineLeaf( leafBlock );
        leafBlock = octree.findLeafBlock(queryPt);
        EXPECT_TRUE( octree.isLeaf(leafBlock));

        leafBB = octree.blockBoundingBox(leafBlock);
        EXPECT_TRUE( leafBB.contains( queryPt ));
        EXPECT_TRUE( bb.contains(leafBB));

        SLIC_INFO("Query pt: " << queryPt
                  <<"\n\t" << ( leafBB.contains(queryPt) ? " was" : " was not")
                  <<" contained in bounding box " << leafBB
                  <<"\n\t of leaf " << leafBlock
                  <<" in the octree. ");
    }


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

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}


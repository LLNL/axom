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
#include "quest/Octree.hpp"

//------------------------------------------------------------------------------
TEST( quest_octree, octree_parent_child)
{
  static const int DIM = 2;
  typedef int CoordType;
  typedef int LeafNodeType;

  typedef quest::TopologicalOctree<DIM, LeafNodeType> OctreeType;
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
TEST( quest_octree, octree_refine)
{
  static const int DIM = 3;
  typedef int LeafNodeType;

  typedef quest::TopologicalOctree<DIM, LeafNodeType> OctreeType;
  typedef OctreeType::BlockIndex BlockIndex;


  std::cout<<"Quest::Octree -- Testing that refining the root block"
          <<" adds all its children to the octree.\n";

  OctreeType octree;

  BlockIndex rootBlock = octree.root();
  octree.refineLeaf( rootBlock );

  int numChildren  = OctreeType::NUM_CHILDREN;
  for(int i = 0; i < numChildren; ++i)
  {
      EXPECT_TRUE( octree.isLeaf( octree.child(rootBlock, i)) );
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


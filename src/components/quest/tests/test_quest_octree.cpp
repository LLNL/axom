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
TEST( quest_numeric_array, octree_parent_child)
{
  static const int DIM = 2;
  typedef int CoordType;

  typedef quest::Point<CoordType, DIM> GridPt;
  typedef quest::TopologicalOctree<DIM, int> OctreeType;
  typedef OctreeType::BlockIndex BlockIndex;

  OctreeType octree;

  BlockIndex rootIndex = octree.root();
  EXPECT_EQ( rootIndex.first, GridPt());
  EXPECT_EQ( rootIndex.second, 0);

  int numChildren  = (1<<DIM);
  for(int i = 0; i < numChildren; ++i)
  {
      // Root's children are at level one and have coordinates that are either zero or one
      BlockIndex childBlock = octree.child( rootIndex.first, rootIndex.second, i);
      EXPECT_EQ( childBlock.second, 1);

      for(int j = 0; j < DIM; ++j)
      {
          int coordVal = childBlock.first[j];
          EXPECT_TRUE(  coordVal == 0 || coordVal == 1);
      }
      // Parent of child is self
      BlockIndex parentBlock = octree.parent( childBlock.first, childBlock.second);
      EXPECT_EQ( parentBlock.first, rootIndex.first);
      EXPECT_EQ( parentBlock.second, rootIndex.second);
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


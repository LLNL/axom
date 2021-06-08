// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/spin/BVHTree.hpp"

using namespace axom;

//------------------------------------------------------------------------------
TEST(spin_bucket_tree, insert_object)
{
  // some typedef short-cuts
  using TreeType = spin::BVHTree<int, 3>;
  using BoxType = primal::BoundingBox<double, 3>;
  using PointType = primal::Point<double, 3>;

  BoxType bb(PointType::zero(), PointType::ones());

  TreeType bucketTree(3, 2);

  EXPECT_TRUE(bucketTree.empty());
  EXPECT_EQ(0, bucketTree.getNumLevels());
  EXPECT_EQ(2, bucketTree.getMaxNumLevels());
  EXPECT_EQ(0, bucketTree.getNumberOfObjects());

  for(int i = 0; i < 5; ++i)
  {
    bucketTree.insert(bb, i);

    EXPECT_FALSE(bucketTree.empty());
    EXPECT_EQ(2, bucketTree.getMaxNumLevels());
    EXPECT_EQ(1, bucketTree.getNumLevels());
    EXPECT_EQ(i + 1, bucketTree.getBucketNumObjects(0));

    bb.expand(0.5);
  }

  const int N = bucketTree.getBucketNumObjects(0);
  EXPECT_EQ(5, N);

  const int* objs = bucketTree.getBucketObjectArray(0);
  ASSERT_TRUE(objs != nullptr);

  const BoxType& bucketBox = bucketTree.getBucketBox(0);
  BoxType expected_bb(PointType::zero(), PointType::ones());
  for(int i = 0; i < 5; ++i)
  {
    const int objIdx = objs[i];

    EXPECT_EQ(expected_bb, bucketTree.getObjectBox(objIdx));
    EXPECT_TRUE(bucketBox.contains(bucketTree.getObjectBox(objIdx)));
    EXPECT_EQ(i, bucketTree.getObjectData(objIdx));

    expected_bb.expand(0.5);
  }
}

//------------------------------------------------------------------------------
#include "axom/slic/core/SimpleLogger.hpp"
using axom::slic::SimpleLogger;

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  SimpleLogger logger;  // create & initialize test logger,

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}

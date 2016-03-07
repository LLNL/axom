/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

/*!
 *******************************************************************************
 * \file quest_bucket_tree.cpp
 *
 * \date Jan 23, 2016
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */
#include "gtest/gtest.h"

// QueST includes
#include "quest/BoundingBox.hpp"
#include "quest/BVHTree.hpp"
#include "quest/Point.hpp"

//------------------------------------------------------------------------------
TEST( quest_bucket_tree, insert_object )
{
  // some typedef short-cuts
  typedef quest::BVHTree< int,3 > TreeType;
  typedef quest::BoundingBox< double,3 > BoxType;
  typedef quest::Point< double,3 > PointType;

  BoxType bb( PointType::zero(), PointType::ones() );

  TreeType bucketTree( 3,2 );

  EXPECT_TRUE( bucketTree.empty() );
  EXPECT_EQ( 0, bucketTree.getNumLevels() );
  EXPECT_EQ( 2, bucketTree.getMaxNumLevels() );
  EXPECT_EQ( 0, bucketTree.getNumberOfObjects() );

  for ( int i=0; i < 5; ++i ) {

    bucketTree.insert( bb, i );

    EXPECT_FALSE( bucketTree.empty() );
    EXPECT_EQ( 2, bucketTree.getMaxNumLevels() );
    EXPECT_EQ( 1, bucketTree.getNumLevels() );
    EXPECT_EQ( i+1, bucketTree.getBucketNumObjects( 0 ) );

    bb.expand( 0.5 );

  }

  const int N     = bucketTree.getBucketNumObjects( 0 );
  EXPECT_EQ( 5, N );

  const int* objs = bucketTree.getBucketObjectArray( 0 );
  ASSERT_TRUE( objs != ATK_NULLPTR );

  const BoxType& bucketBox = bucketTree.getBucketBox( 0 );
  BoxType expected_bb( PointType::zero(), PointType::ones() );
  for ( int i=0; i < 5; ++i ) {

     const int objIdx = objs[ i ];

     EXPECT_EQ( expected_bb, bucketTree.getObjectBox( objIdx ) );
     EXPECT_TRUE( bucketBox.contains( bucketTree.getObjectBox( objIdx ) ) );
     EXPECT_EQ( i, bucketTree.getObjectData( objIdx ) );

     expected_bb.expand( 0.5 );
  }


}

//------------------------------------------------------------------------------
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

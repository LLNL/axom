/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


/*
 *  This file tests the virtual interface of slam sets
 */

#include "gtest/gtest.h"

#include "slam/Set.hpp"
#include "slam/RangeSet.hpp"

static int const NUM_ELEMS = 5;

TEST(gtest_slam_set_virtualbase,construct)
{
  asctoolkit::slam::Set* s = new asctoolkit::slam::RangeSet(NUM_ELEMS);

  // Tests function: isValid()
  EXPECT_TRUE( s->isValid() );

  // Tests function: empty()
  EXPECT_FALSE( s->empty());

  // Tests function: size()
  EXPECT_EQ(NUM_ELEMS, s->size());

  // Tests function: at()
  typedef asctoolkit::slam::Set::IndexType IndexType;
  for(IndexType idx = 0; idx < s->size(); ++idx)
  {
    EXPECT_EQ(idx, s->at(idx));
  }

  delete s;
}


TEST(gtest_slam_set_virtualbase,equality)
{
  asctoolkit::slam::Set* s1 = new asctoolkit::slam::RangeSet(NUM_ELEMS);
  asctoolkit::slam::Set* s2 = new asctoolkit::slam::RangeSet(NUM_ELEMS);
  asctoolkit::slam::Set* s3 = new asctoolkit::slam::RangeSet(2 * NUM_ELEMS);

  // Tests function: isValid()
  EXPECT_TRUE(  s1->isValid() );
  EXPECT_TRUE(  s2->isValid() );
  EXPECT_TRUE(  s3->isValid() );


  // Tests operator==
  EXPECT_TRUE(  *s1 == *s2);
  EXPECT_TRUE(  *s2 == *s1);

  EXPECT_FALSE( *s1 == *s3);
  EXPECT_FALSE( *s2 == *s3);

  // Tests operator!=
  EXPECT_FALSE( *s1 != *s2);
  EXPECT_TRUE(  *s1 != *s3);
  EXPECT_TRUE(  *s2 != *s3);

  // Tests implicit usage of equality operator
  EXPECT_EQ( *s1, *s2);
  EXPECT_NE( *s2, *s3);

  // reclaim memory
  delete s1;
  delete s2;
  delete s3;
}

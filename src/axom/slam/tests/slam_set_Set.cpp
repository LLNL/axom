/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
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


/*
 *  This file tests the virtual interface of slam sets
 */

#include "gtest/gtest.h"

#include "axom/slam/Set.hpp"
#include "axom/slam/RangeSet.hpp"

static int const NUM_ELEMS = 5;

TEST(slam_set_virtualbase,construct)
{
  axom::slam::Set* s = new axom::slam::RangeSet(NUM_ELEMS);

  // Tests function: isValid()
  EXPECT_TRUE( s->isValid() );

  // Tests function: empty()
  EXPECT_FALSE( s->empty());

  // Tests function: size()
  EXPECT_EQ(NUM_ELEMS, s->size());

  // Tests function: at()
  typedef axom::slam::Set::IndexType IndexType;
  for(IndexType idx = 0 ; idx < s->size() ; ++idx)
  {
    EXPECT_EQ(idx, s->at(idx));
  }

  delete s;
}


TEST(slam_set_virtualbase,equality)
{
  axom::slam::Set* s1 = new axom::slam::RangeSet(NUM_ELEMS);
  axom::slam::Set* s2 = new axom::slam::RangeSet(NUM_ELEMS);
  axom::slam::Set* s3 = new axom::slam::RangeSet(2 * NUM_ELEMS);

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

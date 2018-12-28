// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)


/*
 *  This file tests the virtual interface of slam sets
 */

#include "gtest/gtest.h"

#include "axom/slam/Set.hpp"
#include "axom/slam/RangeSet.hpp"

static int const NUM_ELEMS = 5;

namespace slam = axom::slam;
using RangeSetType = slam::RangeSet<>;

TEST(slam_set_virtualbase,construct)
{
  slam::Set* s = new RangeSetType(NUM_ELEMS);

  // Tests function: isValid()
  EXPECT_TRUE( s->isValid() );

  // Tests function: empty()
  EXPECT_FALSE( s->empty());

  // Tests function: size()
  EXPECT_EQ(NUM_ELEMS, s->size());

  // Tests function: at()
  for(slam::Set::PositionType idx = 0 ; idx < s->size() ; ++idx)
  {
    EXPECT_EQ(idx, s->at(idx));
  }

  delete s;
}


TEST(slam_set_virtualbase,equality)
{
  slam::Set* s1 = new RangeSetType(NUM_ELEMS);
  slam::Set* s2 = new RangeSetType(NUM_ELEMS);
  slam::Set* s3 = new RangeSetType(2 * NUM_ELEMS);

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

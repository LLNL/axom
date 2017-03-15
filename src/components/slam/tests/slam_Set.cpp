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
 * testSet.cxx
 *
 *  Created on: Apr 23, 2015
 *      Author: weiss27
 */

#include "gtest/gtest.h"

#include "slam/Set.hpp"
#include "slam/RangeSet.hpp"

static unsigned int const NUM_ELEMS = 5;

TEST(gtest_slam_set,construct_set)
{
  axom::slam::Set* s = new axom::slam::RangeSet(NUM_ELEMS);

  delete s;

  EXPECT_TRUE( true );
}

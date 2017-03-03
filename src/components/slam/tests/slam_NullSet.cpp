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
 * \file testNullSet.cpp
 *
 * Unit tests for the NullSet class
 */

#include "gtest/gtest.h"

#include "slic/slic.hpp"

#include "slam/Set.hpp"
#include "slam/NullSet.hpp"

TEST(gtest_slam_set,construct_nullset)
{
  axom::slam::Set* s = new axom::slam::NullSet();

  EXPECT_TRUE(s->empty());

  delete s;

  EXPECT_TRUE( true );
}

TEST(gtest_slam_set,subscript_fails_nullset)
{
  SLIC_INFO("Testing subscript access on NullSet -- code is expected to assert and die.");

  typedef axom::slam::Set::PositionType SetPosition;
  axom::slam::NullSet n;

  EXPECT_EQ(n.size(), SetPosition()) << "size of null set is defined to be zero";

#ifdef AXOM_DEBUG
  // NOTE: AXOM_DEBUG is disabled in release mode, so this test will only fail in debug mode

  // add this line to avoid a warning in the output about thread safety
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH_IF_SUPPORTED(n[0],"") << "subscript operator on null set asserts";
#else
  SLIC_INFO("Skipped assertion failure check in release mode.");
#endif
}


//----------------------------------------------------------------------
//----------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

int main(int argc, char * argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}

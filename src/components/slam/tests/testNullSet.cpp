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

#include "slam/Set.hpp"
#include "slam/NullSet.hpp"

TEST(gtest_slam_set,construct_nullset)
{
  asctoolkit::slam::Set* s = new asctoolkit::slam::NullSet();

  EXPECT_TRUE(s->empty());

  delete s;

  EXPECT_TRUE( true );
}

TEST(gtest_slam_set,subscript_fails_nullset)
{
  std::cout << "\n****** Testing subscript access on NullSet -- code is expected to assert and die." << std::endl;

  typedef asctoolkit::slam::Set::PositionType SetPosition;
  asctoolkit::slam::NullSet n;

  EXPECT_EQ(n.size(), SetPosition()) << "size of null set is defined to be zero";

#ifdef ATK_DEBUG
  // NOTE: ATK_ASSSERT is disabled in release mode, so this test will only fail in debug mode

  // add this line to avoid a warning in the output about thread safety
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  ASSERT_DEATH(n[0],"") << "subscript operator on null set asserts";
#else
  std::cout << "Did not check for assertion failure since assertions are compiled out in release mode." << std::endl;
#endif
}

/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
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
 * \file slam_set_NullSet.cpp
 *
 * Unit tests for the NullSet class in Slam.
 */

#include "gtest/gtest.h"

#include "slic/slic.hpp"

#include "slam/Set.hpp"
#include "slam/NullSet.hpp"

TEST(slam_set_nullset,construct)
{
  axom::slam::NullSet ns;

  // Test function: isValid()
  EXPECT_TRUE(  ns.isValid() );

  // Test functions: empty() and size()
  EXPECT_TRUE(  ns.empty() );
  EXPECT_EQ(  0,  ns.size() );

  // Test that the parent of a nullset is itself
  EXPECT_EQ(  ns, *ns.parentSet() );

  // A Slam NullSet is not a subset of another set
  EXPECT_FALSE( ns.isSubset() );
}

TEST(slam_set_nullset,subscript_fails)
{
  SLIC_INFO("Testing subscript access on NullSet"
            <<" -- code is expected to assert and die.");

  typedef axom::slam::Set::PositionType SetPosition;
  axom::slam::NullSet n;

  EXPECT_EQ(n.size(), SetPosition())
    << "size of null set is defined to be zero";

#ifdef AXOM_DEBUG
  // NOTE: AXOM_DEBUG is disabled in release mode,
  // so this test will only fail in debug mode

  // add this line to avoid a warning in the output about thread safety
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH_IF_SUPPORTED(n[0],"")
    << "subscript operator on null set asserts";
#else
  SLIC_INFO("Skipped assertion failure check in release mode.");
#endif
}


//----------------------------------------------------------------------
//----------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}

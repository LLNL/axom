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

#include "slic/slic.hpp"

#include "quest/intersection_impl.hpp"

// #include "quest_test_utilities.hpp"


TEST( quest_intersection_impl, fuzzy_comparisons )
{
  const double eps = 0.1;

  SLIC_INFO("This test demonstrates the fuzzy comparison"
      << " operators used in quest's intersection tests"
      << " using a large epsilon (" << eps << ")" );

  // Testing 'greater than' function
  // Note: Excludes the eps window
  EXPECT_TRUE( quest::detail::isGt(1.15, 1.0, eps) );
  // --
  EXPECT_FALSE( quest::detail::isGt(1.05, 1.0, eps) );
  EXPECT_FALSE( quest::detail::isGt(1.0, 1.0, eps) );
  EXPECT_FALSE( quest::detail::isGt(0.95, 1.0, eps) );
  EXPECT_FALSE( quest::detail::isGt(0.85, 1.0, eps) );


  // Testing 'less than or equal' function
  // Note: Includes eps window
  // Note: isLeq() has opposite value of isGt()
  EXPECT_FALSE( quest::detail::isLeq(1.15, 1.0, eps) );
  // --
  EXPECT_TRUE( quest::detail::isLeq(1.05, 1.0, eps) );
  EXPECT_TRUE( quest::detail::isLeq(1.0, 1.0, eps) );
  EXPECT_TRUE( quest::detail::isLeq(0.95, 1.0, eps) );
  EXPECT_TRUE( quest::detail::isLeq(0.85, 1.0, eps) );



  // Testing 'less than' function
  // Note: Excludes eps window
  EXPECT_FALSE( quest::detail::isLt(1.15, 1.0, eps) );
  EXPECT_FALSE( quest::detail::isLt(1.05, 1.0, eps) );
  EXPECT_FALSE( quest::detail::isLt(1.0, 1.0, eps) );
  EXPECT_FALSE( quest::detail::isLt(0.95, 1.0, eps) );
  // --
  EXPECT_TRUE( quest::detail::isLt(0.85, 1.0, eps) );


  // Testing 'greater than or equal' function
  // Note: Includes eps window
  // Note: isGeq() has opposite value of isLt()
  EXPECT_TRUE( quest::detail::isGeq(1.15, 1.0, eps) );
  EXPECT_TRUE( quest::detail::isGeq(1.05, 1.0, eps) );
  EXPECT_TRUE( quest::detail::isGeq(1.0, 1.0, eps) );
  EXPECT_TRUE( quest::detail::isGeq(0.95, 1.0, eps) );
  // --
  EXPECT_FALSE( quest::detail::isGeq(0.85, 1.0, eps) );

}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using asctoolkit::slic::UnitTestLogger;

int main(int argc, char * argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  // create & initialize test logger, finalized when exiting main scope
  UnitTestLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}

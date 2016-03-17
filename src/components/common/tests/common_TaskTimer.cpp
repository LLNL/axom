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
 * \file common_TaskTimer.cpp
 *
 * \date Feb 10, 2016
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */

#include "gtest/gtest.h"

#include "common/Timer.hpp"
#include "slic/slic.hpp"

#ifdef WIN32
#include "windows.h"
void sleep(int numSeconds)
{
  Sleep( numSeconds );
}
#else
#include <unistd.h> // for sleep()
#endif

TEST(gtest_common_Timer, timer_check )
{
  asctoolkit::utilities::Timer t;

  SLIC_INFO("Checking that a newly constructed timer indicates 0 time elapsed");
  EXPECT_EQ(0., t.elapsed());

  t.start();

  sleep( 1 );

  t.stop();

  SLIC_INFO("Simple test for elapsed time in different units.");
  EXPECT_GT( t.elapsedTimeInMicroSec(), t.elapsedTimeInMilliSec() );
  EXPECT_GT( t.elapsedTimeInMilliSec(), t.elapsedTimeInSec() );
  EXPECT_EQ( t.elapsed(), t.elapsedTimeInSec() );


  SLIC_INFO("Testing that reset() indicates zero elapsed time.");
  t.reset();
  ASSERT_DOUBLE_EQ( 0., t.elapsed());
}


TEST(gtest_common_Timer, timer_check_duration )
{
  asctoolkit::utilities::Timer t;
  t.start();

  sleep( 1 );

  t.stop();
  double e = t.elapsed();
  SLIC_INFO( "Elapsed: " << e );

  EXPECT_GE( e, 1.0 );
  EXPECT_LT( e, 2.0 );
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

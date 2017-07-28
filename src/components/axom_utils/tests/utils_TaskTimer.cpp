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
 * \file utils_TaskTimer.cpp
 */

#include "gtest/gtest.h"

#include "axom_utils/Timer.hpp"

#ifdef WIN32
#include "windows.h"
void sleep(int numSeconds)
{
  int numMilliSecs = numSeconds * 1000;
  Sleep( numMilliSecs );
}
#else
#include <unistd.h> // for sleep()
#endif

TEST(gtest_utils_Timer, timer_check )
{
  axom::utilities::Timer t;

  std::cout << "Checking that a newly constructed timer indicates 0 time elapsed" << std::endl;
  EXPECT_EQ(0., t.elapsed());

  t.start();

  sleep( 1 );

  t.stop();

  std::cout << "Simple test for elapsed time in different units." << std::endl;
  EXPECT_GT( t.elapsedTimeInMicroSec(), t.elapsedTimeInMilliSec() );
  EXPECT_GT( t.elapsedTimeInMilliSec(), t.elapsedTimeInSec() );
  EXPECT_EQ( t.elapsed(), t.elapsedTimeInSec() );


  std::cout <<  "Testing that reset() indicates zero elapsed time." << std::endl;
  t.reset();
  ASSERT_DOUBLE_EQ( 0., t.elapsed());
}


TEST(gtest_utils_Timer, timer_check_duration )
{
  axom::utilities::Timer t;
  t.start();

  sleep( 1 );

  t.stop();
  double e = t.elapsed();
  std::cout << "Elapsed: " << e << std::endl;

  EXPECT_GE( e, 1.0 );
  EXPECT_LT( e, 2.0 );
}


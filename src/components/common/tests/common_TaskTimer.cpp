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

#include "common/TaskTimer.hpp"
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

TEST(gtest_common_TaskTimer, timer_check )
{
  asctoolkit::utilities::TaskTimer t;
  t.start();

  sleep( 2 );

  t.stop();
  double e = t.elapsed();
  SLIC_INFO( "Elapsed: " << e );

  EXPECT_GE( e, 2.0 );
  EXPECT_LT( e, 3.0 );
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

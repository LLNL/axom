// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core/utilities/Timer.hpp"

#ifdef WIN32
  #include "windows.h"
void sleep(int numSeconds)
{
  int numMilliSecs = numSeconds * 1000;
  Sleep(numMilliSecs);
}
#else
  #include <unistd.h>  // for sleep()
#endif

TEST(utils_Timer, timer_check)
{
  axom::utilities::Timer t;

  std::cout << "Checking that new timer indicates 0 time elapsed" << std::endl;
  EXPECT_EQ(0., t.elapsed());

  t.start();

  sleep(1);

  t.stop();

  std::cout << "Simple test for elapsed time in different units." << std::endl;
  EXPECT_GT(t.elapsedTimeInMicroSec(), t.elapsedTimeInMilliSec());
  EXPECT_GT(t.elapsedTimeInMilliSec(), t.elapsedTimeInSec());
  EXPECT_EQ(t.elapsed(), t.elapsedTimeInSec());

  std::cout << "Testing that reset() indicates 0 elapsed time." << std::endl;
  t.reset();
  ASSERT_DOUBLE_EQ(0., t.elapsed());
}

TEST(utils_Timer, timer_check_duration)
{
  axom::utilities::Timer t;
  t.start();

  sleep(1);

  t.stop();
  double e = t.elapsed();
  std::cout << "Elapsed: " << e << std::endl;

  EXPECT_GE(e, 1.0);
  EXPECT_LT(e, 2.0);
}

TEST(utils_Timer, timer_check_sum)
{
  const int N = 3;
  axom::utilities::Timer t1(false);
  axom::utilities::Timer t2(false);
  for(int n = 0; n < N; ++n)
  {
    t2.start();
    t1.start();
    sleep(1);
    t1.stop();
    sleep(1);
    t2.stop();
  }

  std::cout << "t1 measured: " << t1.summed() << "s in " << t1.periodCount() << " periods" << std::endl;
  std::cout << "t2 measured: " << t2.summed() << "s in " << t2.periodCount() << " periods" << std::endl;

  EXPECT_EQ(t1.periodCount(), N);
  EXPECT_EQ(t2.periodCount(), N);

  // Estimated inaccuracy per start-stop period.
  // from the C++ code, function call, etc.
  const double tol = 1.4e-4;

  EXPECT_GE(t1.summed()/N, 1 - tol);
  EXPECT_LE(t1.summed()/N, 1 + tol);
  EXPECT_GE(t2.summed()/N, 2 - tol);
  EXPECT_LE(t2.summed()/N, 2 + tol);
}

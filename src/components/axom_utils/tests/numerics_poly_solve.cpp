/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
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

// Axom includes
#include "axom_utils/poly_solve.hpp"

// Google Test include
#include "gtest/gtest.h"

int count_mismatches(double *standard, double *test, int n)
{
  int mcount = 0;

  for (int i = 0; i < n; ++i) {
    mcount += (standard[i] != test[i]);
  }

  return mcount;
}


TEST( numerics_poly_solve, poly_solve_linear )
{
  double t1[2];
  double r1[1];
  double s1[1];
  int n = 1;
  int rc;

  {
    SCOPED_TRACE("Line 1 through origin.");
    t1[0] = 1; t1[1] = 0;
    r1[0] = 0;
    s1[0] = 0;
    rc = axom::numerics::solve_linear(t1, r1, n);
    EXPECT_EQ(0, rc);
    EXPECT_EQ(n, 1);
    EXPECT_EQ(0, count_mismatches(s1, r1, 1));
  }

  {
    SCOPED_TRACE("Line 2 through origin.");
    t1[0] = 18; t1[1] = 0;
    r1[0] = 0;
    s1[0] = 0;
    rc = axom::numerics::solve_linear(t1, r1, n);
    EXPECT_EQ(0, rc);
    EXPECT_EQ(n, 1);
    EXPECT_EQ(0, count_mismatches(s1, r1, 1));
  }

  {
    SCOPED_TRACE("Off origin 1");
    t1[1] = 0.5; t1[0] = -1;
    r1[0] = 0;
    s1[0] = 2;
    rc = axom::numerics::solve_linear(t1, r1, n);
    EXPECT_EQ(0, rc);
    EXPECT_EQ(n, 1);
    EXPECT_EQ(0, count_mismatches(s1, r1, 1));
  }

  {
    SCOPED_TRACE("Off origin 2");
    t1[0] = 0.5; t1[1] = -1;
    r1[0] = 0;
    s1[0] = 2;
    rc = axom::numerics::solve_linear(t1, r1, n);
    EXPECT_EQ(0, rc);
    EXPECT_EQ(n, 1);
    EXPECT_EQ(0, count_mismatches(s1, r1, 1));
  }

  {
    // Not sure exactly how this should work.
    SCOPED_TRACE("X-axis");
    t1[0] = 0; t1[1] = 0;
    r1[0] = 0;
    s1[0] = 0;
    rc = axom::numerics::solve_linear(t1, r1, n);
    EXPECT_EQ(0, rc);
    EXPECT_EQ(n, 1);
    EXPECT_EQ(0, count_mismatches(s1, r1, 1));
  }
}

TEST( numerics_poly_solve, poly_solve_quadratic )
{
}

TEST( numerics_poly_solve, poly_solve_cubic )
{
}


int main(int argc, char * argv[])
{
  ::testing::InitGoogleTest(&argc, argv);

  int result = RUN_ALL_TESTS();
  return result;
}

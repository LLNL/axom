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
  double coeff[2];
  double roots[1];
  double expected[1];
  int n = 1;
  int rc;

  {
    SCOPED_TRACE("Line 1 through origin.");
    coeff[0] = 0; coeff[1] = 1;
    roots[0] = 0;
    expected[0] = 0;
    rc = axom::numerics::solve_linear(coeff, roots, n);
    EXPECT_EQ(rc, 0);
    EXPECT_EQ(n, 1);
    EXPECT_EQ(count_mismatches(expected, roots, 1), 0);
  }

  {
    SCOPED_TRACE("Line 2 through origin.");
    coeff[0] = 0; coeff[1] = 18;
    roots[0] = 0;
    expected[0] = 0;
    rc = axom::numerics::solve_linear(coeff, roots, n);
    EXPECT_EQ(rc, 0);
    EXPECT_EQ(n, 1);
    EXPECT_EQ(count_mismatches(expected, roots, 1), 0);
  }

  {
    SCOPED_TRACE("Off origin 1");
    coeff[0] = -1; coeff[1] = 0.5;
    roots[0] = 0;
    expected[0] = 2;
    rc = axom::numerics::solve_linear(coeff, roots, n);
    EXPECT_EQ(rc, 0);
    EXPECT_EQ(n, 1);
    EXPECT_EQ(count_mismatches(expected, roots, 1), 0);
  }

  {
    SCOPED_TRACE("Off origin 2");
    coeff[0] = 0.5; coeff[1] = -1;
    roots[0] = 0;
    expected[0] = 0.5;
    rc = axom::numerics::solve_linear(coeff, roots, n);
    EXPECT_EQ(rc, 0);
    EXPECT_EQ(n, 1);
    EXPECT_EQ(count_mismatches(expected, roots, 1), 0);
  }

  {
    // Not sure exactly how this should work.
    SCOPED_TRACE("X-axis");
    coeff[0] = 0; coeff[1] = 0;
    roots[0] = 0;
    expected[0] = 0;
    rc = axom::numerics::solve_linear(coeff, roots, n);
    EXPECT_EQ(rc, 0);
    EXPECT_EQ(n, -1);
    EXPECT_EQ(count_mismatches(expected, roots, 1), 0);
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

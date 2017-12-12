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
#include "axom_utils/Utilities.hpp" // for isNearlyEqual()

// Google Test include
#include "gtest/gtest.h"

int count_mismatches(double *standard, double *test, int n)
{
  int mcount = 0;

  for (int i = 0; i < n; ++i) {
    if (! axom::utilities::isNearlyEqual(standard[i], test[i])) {
      ++mcount;
    }
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
    SCOPED_TRACE("X-axis");
    coeff[0] = 0; coeff[1] = 0;
    roots[0] = 0;
    expected[0] = 0;
    rc = axom::numerics::solve_linear(coeff, roots, n);
    // rc == 0 because there are real solutions
    EXPECT_EQ(rc, 0);
    // n == -1 because there are infinitely many solutions
    EXPECT_EQ(n, -1);
    EXPECT_EQ(count_mismatches(expected, roots, 1), 0);
  }
}

TEST( numerics_poly_solve, poly_solve_quadratic )
{
  double coeff[3];
  double roots[2];
  double expected[2];
  int n = 2;
  int rc;

  {
    // y = (x + 2.3)(x + 2.3)
    SCOPED_TRACE("Double root at x = -2.3");
    coeff[0] = 5.29; coeff[1] = 4.6; coeff[2] = 1;
    roots[0] = 0; roots[1] = 0;
    expected[0] = -2.3; expected[1] = -2.3;
    rc = axom::numerics::solve_quadratic(coeff, roots, n);
    EXPECT_EQ(rc, 0);
    EXPECT_EQ(n, 1);
    EXPECT_EQ(count_mismatches(expected, roots, 2), 0);
  }

  {
    // y = (-x + 1.5)(x - 1.5)
    SCOPED_TRACE("Double root at x = 1.5 (opening down)");
    coeff[0] = -2.25; coeff[1] = 3; coeff[2] = -1;
    roots[0] = 0; roots[1] = 0;
    expected[0] = 1.5; expected[1] = 1.5;
    rc = axom::numerics::solve_quadratic(coeff, roots, n);
    EXPECT_EQ(rc, 0);
    EXPECT_EQ(n, 1);
    EXPECT_EQ(count_mismatches(expected, roots, 2), 0);
  }

  {
    // y = 3.2(x + 0.7)(x - 2)
    SCOPED_TRACE("Roots at -0.7 and 2");
    coeff[0] = -4.48; coeff[1] = -4.16; coeff[2] = 3.2;
    roots[0] = 0; roots[1] = 0;
    expected[0] = 2; expected[1] = -0.7;
    rc = axom::numerics::solve_quadratic(coeff, roots, n);
    EXPECT_EQ(rc, 0);
    EXPECT_EQ(n, 2);
    EXPECT_EQ(count_mismatches(expected, roots, 2), 0);
  }

  {
    // y = 0.1x^2 + 0.2x + 6
    SCOPED_TRACE("No real roots (opening up)");
    coeff[0] = 6; coeff[1] = 0.2; coeff[2] = 0.1;
    roots[0] = 0; roots[1] = 0;
    expected[0] = 0; expected[1] = 0;
    rc = axom::numerics::solve_quadratic(coeff, roots, n);
    EXPECT_EQ(rc, -1);
    EXPECT_EQ(n, 0);
    EXPECT_EQ(count_mismatches(expected, roots, 2), 0);
  }

  {
    // y = -5x^2 + 0.2x -20
    SCOPED_TRACE("No real roots (opening up)");
    coeff[0] = 6; coeff[1] = 0.2; coeff[2] = 0.1;
    roots[0] = 0; roots[1] = 0;
    expected[0] = 0; expected[1] = 0;
    rc = axom::numerics::solve_quadratic(coeff, roots, n);
    EXPECT_EQ(rc, -1);
    EXPECT_EQ(n, 0);
    EXPECT_EQ(count_mismatches(expected, roots, 2), 0);
  }
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

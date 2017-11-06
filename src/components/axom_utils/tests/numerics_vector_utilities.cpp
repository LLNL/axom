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

#include "gtest/gtest.h"

#include "axom_utils/vector_utilities.hpp"
#include "axom_utils/Utilities.hpp"

namespace numerics = axom::numerics;
namespace utilities = axom::utilities;

TEST( numerics_vector_utilities, dot_product_test )
{
  const int dim = 10;
  const double EPS = 1E-4;

  double * u = new double[dim]();  // 0 initialization
  double * v = new double[dim]();  // 0 initialization

  u[0] = 1.;
  v[1] = 1.;

  EXPECT_NEAR(numerics::dot_product< double >(u, v, dim), 0., EPS);

  u[1] = 1.;
  v[0] = 1.;

  EXPECT_NEAR(numerics::dot_product< double >(u, v, dim), 2., EPS);

  u[5] = 10.;
  v[5] = -10.;

  EXPECT_NEAR(numerics::dot_product< double >(u, v, dim), -98., EPS);

  // free up allocated memory
  delete [] u;
  delete [] v;
}

TEST( numerics_vector_utilities, make_orthogonal_test )
{
  const int dim = 3;
  const double ONE_THIRD = 0.3333333333;
  const double EPS = 1E-4;

  double * u = new double[dim]();  // 0 initialization
  double * v = new double[dim]();  // 0 initialization

  u[0] = 1.;
  u[1] = 1.;
  u[2] = -1.;

  v[0] = 1.;
  v[1] = 1.;
  v[2] = 1.;

  numerics::make_orthogonal< double >(u, v, dim);

  EXPECT_NEAR(u[0], 1. - ONE_THIRD, EPS);
  EXPECT_NEAR(u[1], 1. - ONE_THIRD, EPS);
  EXPECT_NEAR(u[2], -1. - ONE_THIRD, EPS);

  numerics::make_orthogonal(u, u, dim);

  EXPECT_NEAR(u[0], 0., EPS);
  EXPECT_NEAR(u[1], 0., EPS);
  EXPECT_NEAR(u[2], 0., EPS);

  // free up allocated memory
  delete [] u;
  delete [] v;
}

TEST( numerics_vector_utilities, orthonormalize_test )
{
  const int dim = 3;
  const int size = 2;
  const double EPS = 1E-4;
  const double ONE_OVER_SQRT_TWO = 0.70710678;

  double * basis = new double[dim*dim]();  // 0 initialization

  double * u = basis;
  double * v = basis + dim;

  u[0] = 1.;
  u[1] = 1.;
  v[0] = 2.;

  EXPECT_TRUE(numerics::orthonormalize< double >(basis, size, dim));

  EXPECT_NEAR(u[0], ONE_OVER_SQRT_TWO, EPS);
  EXPECT_NEAR(u[1], ONE_OVER_SQRT_TWO, EPS);
  EXPECT_NEAR(u[2], 0., EPS);

  EXPECT_NEAR(v[0], ONE_OVER_SQRT_TWO, EPS);
  EXPECT_NEAR(v[1], -ONE_OVER_SQRT_TWO, EPS);
  EXPECT_NEAR(v[2], 0., EPS);

  // doesn't span
  EXPECT_FALSE(numerics::orthonormalize< double >(basis, dim, dim));

  // free up allocated memory
  delete [] basis;
}

TEST( numerics_vector_utilities, normalize_test )
{
  const int dim = 3;
  const double EPS = 1E-4;
  const double ONE_OVER_SQRT_THREE = 0.57735026;
  const double ONE_OVER_SQRT_TEN = 0.316227766;

  double * u = new double[dim]();  // 0 initialization
  double * v = new double[dim]();  // 0 initialization

  u[0] = 1.;
  u[1] = 1.;
  u[2] = 1.;

  v[0] = 0.;
  v[1] = 3.;
  v[2] = 1.;

  EXPECT_TRUE(numerics::normalize< double >(u, dim));

  EXPECT_NEAR(u[0], ONE_OVER_SQRT_THREE, EPS);
  EXPECT_NEAR(u[1], ONE_OVER_SQRT_THREE, EPS);
  EXPECT_NEAR(u[2], ONE_OVER_SQRT_THREE, EPS);

  EXPECT_TRUE(numerics::normalize< double >(v, dim));

  EXPECT_NEAR(v[0], 0., EPS);
  EXPECT_NEAR(v[1], 3.*ONE_OVER_SQRT_TEN, EPS);
  EXPECT_NEAR(v[2], ONE_OVER_SQRT_TEN, EPS);

  // now make u = (0, 0, 0)
  numerics::make_orthogonal(u, u, dim);

  // normalize should fail
  EXPECT_FALSE(numerics::normalize< double >(u, dim));
}

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

#include "axom_utils/vector_utilities.hpp"
#include "axom_utils/Utilities.hpp"

namespace numerics = axom::numerics;
namespace utilities = axom::utilities;

TEST( numerics_vector_utilities, dot_product_test )
{
  int dim = 10;

  double *u = new double[dim];
  double *v = new double[dim];

  u[0] = 1.;
  v[1] = 1.;

  EXPECT_TRUE(utilities::isNearlyEqual< double >(
    numerics::dot_product< double >(u, v, dim), 0.));

  u[1] = 1.;
  v[0] = 1.;

  EXPECT_TRUE(utilities::isNearlyEqual< double >(
    numerics::dot_product< double >(u, v, dim), 2.));

  u[5] = 10.;
  v[5] = -10.;

  EXPECT_TRUE(utilities::isNearlyEqual< double >(
    numerics::dot_product< double >(u, v, dim), -98.));

  // free up allocated memory
  delete [] u;
  delete [] v;
}


TEST( numerics_vector_utilities, make_orthogonal_test )
{

  int dim = 3;

  double *u = new double[dim];
  double *v = new double[dim];

  const double ONE_THIRD = 0.3333333333;

  u[0] = 1.;
  u[1] = 1.;
  u[2] = -1.;

  v[0] = 1.;
  v[1] = 1.;
  v[2] = 1.;

  numerics::make_orthogonal< double >(u, v, dim);

  EXPECT_TRUE(utilities::isNearlyEqual(u[0], 1. - ONE_THIRD, 1e-4));
  EXPECT_TRUE(utilities::isNearlyEqual(u[1], 1. - ONE_THIRD, 1e-4));
  EXPECT_TRUE(utilities::isNearlyEqual(u[2], -1. - ONE_THIRD, 1e-4));

  numerics::make_orthogonal(u, u, dim);

  EXPECT_TRUE(utilities::isNearlyEqual(u[0], 0., 1e-4));
  EXPECT_TRUE(utilities::isNearlyEqual(u[1], 0., 1e-4));
  EXPECT_TRUE(utilities::isNearlyEqual(u[2], 0., 1e-4));

  // free up allocated memory
  delete [] u;
  delete [] v;
}

TEST( numerics_vector_utilities, orthonormalize_test ) {

  int dim = 3;

  int size = 2;

  double *basis = new double[dim*dim];

  double *u = basis;
  double *v = basis + dim;

  const double ONE_OVER_SQRT_TWO = 0.70710678;

  u[0] = 1.;
  u[1] = 1.;
  v[0] = 2.;
  v[1] = 0.;

  EXPECT_TRUE(numerics::orthonormalize< double >(basis, size, dim));

  EXPECT_TRUE(utilities::isNearlyEqual(u[0], ONE_OVER_SQRT_TWO, 1e-4));
  EXPECT_TRUE(utilities::isNearlyEqual(u[1], ONE_OVER_SQRT_TWO, 1e-4));
  EXPECT_TRUE(utilities::isNearlyEqual(u[2], 0., 1e-4));

  EXPECT_TRUE(utilities::isNearlyEqual(v[0], ONE_OVER_SQRT_TWO, 1e-4));
  EXPECT_TRUE(utilities::isNearlyEqual(v[1], -ONE_OVER_SQRT_TWO, 1e-4));
  EXPECT_TRUE(utilities::isNearlyEqual(v[2], 0., 1e-4));

  // doesn't span
  EXPECT_FALSE(numerics::orthonormalize< double >(basis, dim, dim));

  // free up allocated memory
  delete [] basis;
}

TEST( numerics_vector_utilities, normalize_test ) {

  int dim = 3;

  double *u = new double[dim];
  double *v = new double[dim];

  u[0] = 1.;
  u[1] = 1.;
  u[2] = 1.;

  v[0] = 0.;
  v[1] = 3.;
  v[2] = 1.;

  const double ONE_OVER_SQRT_THREE = 0.57735026;
  const double ONE_OVER_SQRT_TEN = 0.316227766;

  EXPECT_TRUE(numerics::normalize< double >(u, dim));

  EXPECT_TRUE(utilities::isNearlyEqual(u[0], ONE_OVER_SQRT_THREE, 1e-4));
  EXPECT_TRUE(utilities::isNearlyEqual(u[1], ONE_OVER_SQRT_THREE, 1e-4));
  EXPECT_TRUE(utilities::isNearlyEqual(u[2], ONE_OVER_SQRT_THREE, 1e-4));

  EXPECT_TRUE(numerics::normalize< double >(v, dim));

  EXPECT_TRUE(utilities::isNearlyEqual(v[0], 0., 1e-4));
  EXPECT_TRUE(utilities::isNearlyEqual(v[1], 3.*ONE_OVER_SQRT_TEN, 1e-4));
  EXPECT_TRUE(utilities::isNearlyEqual(v[2], ONE_OVER_SQRT_TEN, 1e-4));

  // now make u = (0, 0, 0)
  numerics::make_orthogonal(u, u, dim);

  // normalize should fail
  EXPECT_FALSE(numerics::normalize< double >(u, dim));
}

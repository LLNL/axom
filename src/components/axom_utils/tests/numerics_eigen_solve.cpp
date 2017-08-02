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

#include "axom_utils/Matrix.hpp"
#include "axom_utils/Utilities.hpp"
#include "axom_utils/eigen_solve.hpp"

namespace numerics = axom::numerics;
namespace utilities = axom::utilities;

TEST( numerics_eigen_solve, eigen_solve_with_diagonal_matrix )
{
  const int depth = 1000;
  const int N = 10;

  double* u = new double[N*N];
  double* lambdas = new double[N];

  numerics::Matrix< double > A(N, N);

  // make a diagonal matrix with descending values on the diagonal
  for (int i = 0; i < N; i++)
    A(i,i) = N - i;

  EXPECT_EQ(1, numerics::eigen_solve(A, N, u, lambdas, depth));

  // eigenvalues should be N, N - 1, ... etc. vectors should be standard basis
  // vectors in standard order
  for (int i = 0; i < N; i++) {
    EXPECT_TRUE(utilities::isNearlyEqual(lambdas[i], static_cast<double>(N - i)));
    for (int j = 0; j < N; j++) {
      if (j == i) {
        EXPECT_TRUE(utilities::isNearlyEqual(u[i*N + j]*u[i*N + j], 1.));
      } else {
        EXPECT_TRUE(utilities::isNearlyEqual(u[i*N + j], 0.));
      }
    }
  }

  delete [] u;
  delete [] lambdas;
}

TEST( numerics_eigen_solve, eigen_solve_with_partial_diagonal )
{
  const int depth = 1000;
  const int N = 10;

  double* u = new double[N*N];
  double* lambdas = new double[N];

  numerics::Matrix< double > A(N, N);

  // same thing as above test, but zero out the diagonal halfway through
  for (int i = 0; i < N/2; i++)
    A(i,i) = N - i;

  EXPECT_EQ(1, numerics::eigen_solve(A, N, u, lambdas, depth));

  // now eigenvals should be N, N - 1, ... N/2 + 1, and 0's for the rest
  // vecs are same as before for first half. After that any basis for the
  // remaining subspace would be valid; we'll verify the last half of
  // eigenvecs are orthonormal

  for (int i = 0; i < N/2; i++) {
    EXPECT_TRUE(utilities::isNearlyEqual(lambdas[i], static_cast<double>(N - i)));
    for (int j = 0; j < N; j++) {
      if (j == i) {
        EXPECT_TRUE(utilities::isNearlyEqual(u[i*N + j]*u[i*N + j], 1.));
      } else {
        EXPECT_TRUE(utilities::isNearlyEqual(u[i*N + j], 0.));
      }
    }
  }

  for (int i = N/2; i < N; i++) {
    EXPECT_TRUE(utilities::isNearlyEqual(lambdas[i], 0.));
    for (int j = N/2; j < N; j++) {
      double dot = numerics::dot_product< double >(u + i*N, u + j*N, N);
      if (j == i) {  // check it's unit norm
        EXPECT_TRUE(utilities::isNearlyEqual(dot, 1.));
      } else {  // check they're orthogonal
        EXPECT_TRUE(utilities::isNearlyEqual(dot, 0.));
      }
    }
  }

  delete [] u;
  delete [] lambdas;
}

TEST( numerics_eigen_solve, eigen_solve_with_two_by_two )
{
  const int depth = 1000;
  const int N = 2;

  double* u = new double[N*N];
  double* lambdas = new double[N];

  numerics::Matrix< double > A(N, N);

  A(0,0) = 3.;
  A(0,1) = 1.;
  A(1,0) = 1.;
  A(1,1) = 3.;

  int rc = numerics::eigen_solve(A, N, u, lambdas, depth);
  EXPECT_EQ(1, rc);

  // check lambdas are correct
  EXPECT_TRUE(utilities::isNearlyEqual(lambdas[0], 4.));
  EXPECT_TRUE(utilities::isNearlyEqual(lambdas[1], 2.));

  // check it has correct eigen vecs; we're doing this weird squaring thing
  // because they could be in one of two orders; namely (1/sqrt(2), 1/sqrt(2))
  // followed by (1/sqrt(2), -1/sqrt(2)) or some rearrangement or multiplied
  // by -1
  EXPECT_TRUE(utilities::isNearlyEqual(u[0]*u[0], 0.5, 1.0E-3));
  EXPECT_TRUE(utilities::isNearlyEqual(u[1]*u[1], 0.5, 1.0E-3));
  EXPECT_TRUE(utilities::isNearlyEqual(u[2]*u[2], 0.5, 1.0E-3));
  EXPECT_TRUE(utilities::isNearlyEqual(u[3]*u[3], 0.5, 1.0E-3));

  delete [] u;
  delete [] lambdas;
}

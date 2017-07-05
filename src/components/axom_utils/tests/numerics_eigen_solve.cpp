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

  for (int i = 0; i < N; i++) {
    A(i,i) = N - i;
  }

  int rc = numerics::eigen_solve(A, N, depth, u, lambdas);
  EXPECT_EQ(0, rc);

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

TEST( numerics_eigen_solve, eigen_solve_with_two_by_two ) {
  const int depth = 1000;
  const int N = 2;

  double* u = new double[N*N];
  double* lambdas = new double[N];

  numerics::Matrix< double > A(N, N);

  A(0,0) = 3.;
  A(0,1) = 1.;
  A(1,0) = 1.;
  A(1,1) = 3.;

  int rc = numerics::eigen_solve(A, N, depth, u, lambdas);
  EXPECT_EQ(0, rc);

  for (int i = 0; i < N; i++) {
    std::cout << "Find " << i << "th eigenval = " << lambdas[i] << std::endl;
  }

  // check lambdas are correct
  EXPECT_TRUE(utilities::isNearlyEqual(lambdas[0], 4.));
  EXPECT_TRUE(utilities::isNearlyEqual(lambdas[1], 2.));

  // check it has correct eigen vecs approximately
  EXPECT_TRUE(utilities::isNearlyEqual(u[0], 0.7071, 1.0E-3));
  EXPECT_TRUE(utilities::isNearlyEqual(u[1], 0.7071, 1.0E-3));
  EXPECT_TRUE(utilities::isNearlyEqual(u[2], 0.7071, 1.0E-3));
  EXPECT_TRUE(utilities::isNearlyEqual(u[3], -0.7071, 1.0E-3));

  delete [] u;
  delete [] lambdas;
}

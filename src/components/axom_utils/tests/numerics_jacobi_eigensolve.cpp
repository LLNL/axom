/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
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
#include "axom/Types.hpp"                    // for AXOM_NULLPTR

#include "axom_utils/Matrix.hpp"             // for numerics::Matrix
#include "axom_utils/Utilities.hpp"          // for random_real()
#include "axom_utils/jacobi_eigensolve.hpp"  // for jacobi_eigensolve()
#include "axom_utils/matvecops.hpp"          // for matrix operators

// gtest includes
#include "gtest/gtest.h"

// namespace aliases
namespace numerics  = axom::numerics;
namespace utilities = axom::utilities;

TEST( numerics_jacobi_eigensolve, random_symmetric_matrix )
{
  using IndexType = typename numerics::Matrix< double >::IndexType;

  // 40x40 test. Generate random symmetric matrix A
  int N = 40;
  numerics::Matrix<double> A_test(N,N);
  for (IndexType i=0; i < N; ++i)
  {
    for (IndexType j=0; j <= i; ++j)
    {
      // A_test(i,j) = (double) rand()/ (double) RAND_MAX;
      A_test(i,j) = utilities::random_real<double> (0.0, 1.0);
      A_test(j,i) = A_test(i,j);
    }
  }

  double *lambdas = new double[N];
  numerics::Matrix<double> V(N,N);

  numerics::jacobi_eigensolve<double>( A_test, V, lambdas);

  // compute A*V - V*lambda
  numerics::Matrix<double> tmp(N,N);
  numerics::Matrix<double> test(N,N);
  numerics::matrix_multiply(A_test, V, tmp);

  for (IndexType i=0; i < N; ++i)
  {
    for (IndexType j=0; j < N; ++j)
    {
      test(i,j) = tmp(i,j) - V(i,j)*lambdas[j];
    }
  }

  double p1norm    = numerics::matrix_norm( test, numerics::P1_NORM );
  double inftynorm = numerics::matrix_norm( test, numerics::INF_NORM );
  double frobnorm  = numerics::matrix_norm( test, numerics::FROBENIUS_NORM );

  std::cout << "p1norm: "    << p1norm    << std::endl;
  std::cout << "inftynorm: " << inftynorm << std::endl;
  std::cout << "frobnorm: "  << frobnorm  << std::endl;
}

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  int result = RUN_ALL_TESTS();
  return result;
}

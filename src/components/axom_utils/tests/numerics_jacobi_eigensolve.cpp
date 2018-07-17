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
#include "axom_utils/Utilities.hpp"          // for random_real()/isNearlyEqual
#include "axom_utils/jacobi_eigensolve.hpp"  // for jacobi_eigensolve()
#include "axom_utils/matvecops.hpp"          // for matrix operators

// gtest includes
#include "gtest/gtest.h"

// C/C++ includes
#include <sstream> // for std::ostringstream

// namespace aliases
namespace numerics  = axom::numerics;
namespace utilities = axom::utilities;

TEST( numerics_jacobi_eigensolve, random_symmetric_matrix )
{
  using IndexType = typename numerics::Matrix< double >::IndexType;
  constexpr unsigned int seed = 123456789;
  constexpr double TOL        = 1.e-12;
  constexpr int MAX_SIZE      = 41;

  // 40x40 test. Generate random symmetric matrix A
  for ( int k=2 ; k < MAX_SIZE ; ++k  )
  {
    const int N = k;

    std::ostringstream oss;
    oss << "checking [" << N << " x " << N << "] matrix";
    SCOPED_TRACE( oss.str() );

    numerics::Matrix<double> A_test(N,N);
    for (IndexType i=0 ; i < N ; ++i)
    {
      for (IndexType j=0 ; j <= i ; ++j)
      {
        A_test(i,j) = utilities::random_real<double> (0.0, 1.0, seed);
        A_test(j,i) = A_test(i,j);
      }
    }

    double* lambdas = new double[N];
    numerics::Matrix<double> V(N,N);

    int numIterations = 0;
    numerics::jacobi_eigensolve( A_test, V, lambdas,
                                 numerics::JACOBI_DEFAULT_MAX_ITERATIONS,
                                 &numIterations,
                                 numerics::JACOBI_DEFAULT_TOLERANCE );

    EXPECT_TRUE( numIterations > 0 );
    EXPECT_TRUE( numIterations < 16 );

    // compute A*V - V*lambda
    numerics::Matrix<double> tmp(N,N);
    numerics::Matrix<double> test(N,N);
    numerics::matrix_multiply(A_test, V, tmp);

    for (IndexType i=0 ; i < N ; ++i)
    {
      for (IndexType j=0 ; j < N ; ++j)
      {
        test(i,j) = tmp(i,j) - V(i,j)*lambdas[j];
        EXPECT_TRUE( utilities::isNearlyEqual( test(i,j), 0.0, TOL ) );
      }
    }

    double p1norm    = numerics::matrix_norm( test, numerics::P1_NORM );
    double inftynorm = numerics::matrix_norm( test, numerics::INF_NORM );
    double frobnorm  = numerics::matrix_norm( test, numerics::FROBENIUS_NORM );

    EXPECT_TRUE( utilities::isNearlyEqual( p1norm, 0.0, TOL ) );
    EXPECT_TRUE( utilities::isNearlyEqual( inftynorm, 0.0, TOL ) );
    EXPECT_TRUE( utilities::isNearlyEqual( frobnorm, 0.0, TOL ) );

    delete [] lambdas;
    lambdas = AXOM_NULLPTR;

  } // END for all matrix sizes

}

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  int result = RUN_ALL_TESTS();
  return result;
}

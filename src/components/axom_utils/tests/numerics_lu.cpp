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

#include "gtest/gtest.h"

#include "axom_utils/numerics/LU.hpp"
#include "axom_utils/numerics/Matrix.hpp"
#include "axom_utils/numerics/matvecops.hpp"

namespace numerics = axom::numerics;

TEST( numerics_lu, lu_decompose_non_square_matrix )
{
  int pivots[3];
  const int M = 2;
  const int N = 3;

  EXPECT_TRUE( M != N );
  numerics::Matrix< double > A = numerics::Matrix< double >::ones(M,N);
  int rc = numerics::lu_decompose( A, pivots );
  EXPECT_EQ( numerics::LU_NONSQUARE_MATRIX, rc );
}

//------------------------------------------------------------------------------
TEST( numerics_lu, lu_decompose_singular_matrix )
{
  int pivots[3];
  const int N = 3;

  numerics::Matrix< double > A = numerics::Matrix< double >::zeros( N, N );
  int rc = numerics::lu_decompose( A, pivots );
  EXPECT_EQ( numerics::LU_SINGULAR_MATRIX, rc );
}

//------------------------------------------------------------------------------
TEST( numerics_lu, lu_decompose )
{
  int pivots[3];
  const int N = 3;

  // initial matrix
  numerics::Matrix< double > A( N,N );
  A( 0,0 ) = 1; A( 0,1 ) = (-2); A( 0,2 ) = 3;
  A( 1,0 ) = 2; A( 1,1 ) = (-5); A( 1,2 ) = 12;
  A( 2,0 ) = 0; A( 2,1 ) =    2; A( 2,2 ) = (-10);

  // decompose matrix
  numerics::Matrix< double > A_decomposed = A;
  int rc = numerics::lu_decompose( A_decomposed, pivots );
  EXPECT_EQ( numerics::LU_SUCCESS, rc );

  // extract lower/upper triangular matrix
  numerics::Matrix< double > L = numerics::lower_triangular( A_decomposed );
  numerics::Matrix< double > U = numerics::upper_triangular( A_decomposed );

  // construct permutation matrix on lhs
  numerics::Matrix< double > P = numerics::Matrix< double >::identity( N );
  for ( int i=0 ; i < N ; ++i )
  {
    if ( pivots[i] != i )
    {
      P.swapRows( i, pivots[i] );
    }
  }

  numerics::Matrix< double > PA( N, N );
  numerics::matrix_multiply( P, A, PA );

  numerics::Matrix< double > LU( N, N );
  numerics::matrix_multiply( L, U, LU );

  for ( int i=0 ; i < N ; ++i )
  {
    for ( int j=0 ; j < N ; ++j )
    {
      EXPECT_EQ( PA(i,j), LU(i,j) );
    } // END for all colummns
  } // END for all rows

}
//------------------------------------------------------------------------------
TEST( numerics_lu, lu_solve )
{
  int pivots[3];
  const int N = 3;

  // initial matrix
  numerics::Matrix< double > A( N,N );
  A( 0,0 ) = 1; A( 0,1 ) = 2; A( 0,2 ) = 4;
  A( 1,0 ) = 3; A( 1,1 ) = 8; A( 1,2 ) = 14;
  A( 2,0 ) = 2; A( 2,1 ) = 6; A( 2,2 ) = 13;

  double b[3] = { 3, 13, 4 };
  double x[3];
  numerics::Matrix< double > A_decomposed = A;
  int rc = numerics::lu_decompose( A_decomposed, pivots );
  EXPECT_EQ( numerics::LU_SUCCESS, rc );

  rc = numerics::lu_solve( A_decomposed, pivots, b, x );
  EXPECT_EQ( numerics::LU_SUCCESS, rc );
  EXPECT_DOUBLE_EQ( 3,  x[0] );
  EXPECT_DOUBLE_EQ( 4,  x[1] );
  EXPECT_DOUBLE_EQ( -2, x[2] );

  double rhs[ 3 ];
  numerics::matrix_vector_multiply( A, x, rhs );

  using IndexType = numerics::Matrix< double >::IndexType;
  for ( IndexType i=0 ; i < N ; ++i )
  {
    EXPECT_DOUBLE_EQ( b[i], rhs[ i ] );
  }

}

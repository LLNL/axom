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

#include "axom_utils/Matrix.hpp"
#include "axom_utils/linear_solve.hpp"

namespace numerics = axom::numerics;

TEST( numerics_linear_solve, linear_solve_with_identity_matrix )
{
  const int N = 25;

  for ( int i=1 ; i < N ; ++i )
  {

    numerics::Matrix< double > A = numerics::Matrix< double >::identity( i );

    // form right-hand side
    double* b = new double[ i ];
    for ( int j=0 ; j < i ; ++j )
    {
      b[ j ] = j;
    }

    // allocate solution vector
    double* x = new double[ i ];

    // solve
    int rc = numerics::linear_solve( A, b, x );
    EXPECT_EQ( 0, rc );

    // check answer
    for ( int j=0 ; j < i ; ++j )
    {
      EXPECT_DOUBLE_EQ( j, x[j] );
    }

    delete [] b;
    delete [] x;
  }

}

//------------------------------------------------------------------------------
TEST( numerics_linear_solve, linear_solve2x2 )
{
  const int N = 2;
  numerics::Matrix< double > A( N,N );
  A( 0,0 )=2; A( 0,1 )=1;
  A( 1,0 )=3; A( 1,1 )=4;

  double b[ 2 ] = { 1, 14 };
  double x[ 2 ];

  int rc = numerics::linear_solve( A, b, x );
  EXPECT_EQ( 0, rc );

  EXPECT_DOUBLE_EQ( -2, x[0] );
  EXPECT_DOUBLE_EQ(  5, x[1] );
}

//------------------------------------------------------------------------------
TEST( numerics_linear_solve, linear_solve3x3 )
{
  const int N = 3;
  numerics::Matrix< double > A( N,N );
  A( 0,0 )=4; A( 0,1 )=5; A( 0,2 )=-2;
  A( 1,0 )=7; A( 1,1 )=-1; A( 1,2 )= 2;
  A( 2,0 )=3; A( 2,1 )=1; A( 2,2 )= 4;

  double b[ 3 ] = { -14, 42, 28 };
  double x[ 3 ];

  int rc = numerics::linear_solve( A, b, x );
  EXPECT_EQ( 0, rc );

  EXPECT_DOUBLE_EQ(  4, x[0] );
  EXPECT_DOUBLE_EQ( -4, x[1] );
  EXPECT_DOUBLE_EQ(  5, x[2] );
}

//------------------------------------------------------------------------------
TEST( numerics_linear_solve, linear_solve_singular )
{
  const int N = 3;
  numerics::Matrix< double > A  = numerics::Matrix< double >::zeros( N,N );
  double b[ 3 ] = { -14, 42, 28 };
  double x[ 3 ];

  int rc = numerics::linear_solve( A, b, x );
  EXPECT_NE( 0, rc );
}

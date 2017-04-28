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

// C/C++ includes
#include <sstream> // for std::ostringstream

namespace numerics = axom::numerics;

//-----------------------------------------------------------------------------
// HELPER ROUTINES
//-----------------------------------------------------------------------------
namespace {

void testConstAccess( const numerics::Matrix< double >& A)
{
  const int MROWS = 10;
  const int NCOLS = 10;
  EXPECT_EQ( MROWS, A.getNumRows() );
  EXPECT_EQ( NCOLS, A.getNumColumns() );

  for ( int irow=0; irow < MROWS; ++irow ) {
     for ( int jcol=0; jcol < NCOLS; ++jcol ) {
        double expval = static_cast< double >( irow*NCOLS+jcol );
        EXPECT_EQ( expval, A(irow,jcol) );
     } // END for all columns
  } // END for all rows

}

//------------------------------------------------------------------------------
void testCopyConstructor( numerics::Matrix< double > A )
{
  const int MROWS = 10;
  const int NCOLS = 10;
  EXPECT_EQ( MROWS, A.getNumRows() );
  EXPECT_EQ( NCOLS, A.getNumColumns() );

  for ( int irow=0; irow < MROWS; ++irow ) {
     for ( int jcol=0; jcol < NCOLS; ++jcol ) {
        double expval = static_cast< double >( irow*NCOLS+jcol );
        EXPECT_EQ( expval, A(irow,jcol) );
     } // END for all columns
  } // END for all rows

}

} /* end unnamed namespace */

//-----------------------------------------------------------------------------
// TEST ROUTINES
//------------------------------------------------------------------------------
TEST( numerics_matrix, basic_constructor )
{
  const int MROWS = 5;
  const int NCOLS = 10;

  numerics::Matrix< double > A( MROWS, NCOLS );
  EXPECT_EQ( MROWS, A.getNumRows() );
  EXPECT_EQ( NCOLS, A.getNumColumns() );
  EXPECT_FALSE( A.isSquare() );

  for ( int irow=0; irow < MROWS; ++irow ) {
     for ( int jcol=0; jcol < NCOLS; ++jcol ) {
        EXPECT_DOUBLE_EQ( 0.0, A(irow,jcol) );
     }
  }
}

//------------------------------------------------------------------------------
TEST( numerics_matrix, basic_constructor_with_value )
{
  const int MROWS = 5;
  const int NCOLS = 10;
  const double FILL_VAL = 2.5;

  numerics::Matrix< double > A( MROWS, NCOLS, FILL_VAL );
  EXPECT_EQ( MROWS, A.getNumRows() );
  EXPECT_EQ( NCOLS, A.getNumColumns() );
  EXPECT_FALSE( A.isSquare() );

  for ( int irow=0; irow < MROWS; ++irow ) {
    for ( int jcol=0; jcol < NCOLS; ++jcol ) {
       EXPECT_DOUBLE_EQ( FILL_VAL, A(irow,jcol) );
     }
  }

}

//------------------------------------------------------------------------------
TEST( numerics_matrix, array_constructor )
{
  const int MROWS = 2;
  const int NCOLS = 2;
  int data[ 4 ]   = { 1, 2, 3, 4 };

  numerics::Matrix< int > A( MROWS, NCOLS, data );
  EXPECT_EQ( MROWS, A.getNumRows() );
  EXPECT_EQ( NCOLS, A.getNumColumns() );

  int idx = 0;
  for ( int j=0; j < NCOLS; ++j ) {
     for ( int i=0; i < MROWS; ++i ) {
        ++idx;
        EXPECT_EQ( idx, A(i,j) );
     }
  }

}

//------------------------------------------------------------------------------
TEST( numerics_matrix, is_square )
{
  const int MROWS = 10;
  const int NCOLS = 10;

  numerics::Matrix< double > A( MROWS, NCOLS );
  EXPECT_TRUE( A.isSquare() );
}

//------------------------------------------------------------------------------
TEST( numerics_matrix, random_access_operators )
{
  const int MROWS = 10;
  const int NCOLS = 10;

  numerics::Matrix< double > A( MROWS, NCOLS );

  for ( int irow=0; irow < MROWS; ++irow ) {
     for ( int jcol=0; jcol < NCOLS; ++jcol ) {

        double val     = static_cast< double >( irow*NCOLS+jcol );
        A( irow,jcol ) = val;
        EXPECT_EQ( val, A(irow,jcol) );
     } // END for all columns
  } // END for all rows

  testConstAccess( A );

}

//------------------------------------------------------------------------------
TEST( numerics_matrix, copy_constructor )
{
  const int MROWS = 10;
  const int NCOLS = 10;

  numerics::Matrix< double > A( MROWS, NCOLS );

  for ( int irow=0; irow < MROWS; ++irow ) {
     for ( int jcol=0; jcol < NCOLS; ++jcol ) {

        double val     = static_cast< double >( irow*NCOLS+jcol );
        A( irow,jcol ) = val;
        EXPECT_EQ( val, A(irow,jcol) );
     } // END for all columns
  } // END for all rows

  testCopyConstructor( A );
}

//------------------------------------------------------------------------------
TEST( numerics_matrix, assignment )
{
  const int MROWS = 3;
  const int NCOLS = 3;

  numerics::Matrix< int > A( MROWS, NCOLS );
  A.fill( 3 );

  numerics::Matrix< int > B( 2, 2 );
  B.fill( 1 );

  B = A;

  EXPECT_EQ( MROWS, B.getNumRows() );
  EXPECT_EQ( NCOLS, B.getNumColumns() );

  for ( int i=0; i < MROWS; ++i ) {
     for ( int j=0; j < NCOLS; ++j ) {
         EXPECT_EQ( A(i,j), B(i,j) );
     }
  }
}

//------------------------------------------------------------------------------
TEST( numerics_matrix, getColumn )
{
  const int N=3;
  numerics::Matrix< int > M = numerics::Matrix< int >::identity( N );
  EXPECT_EQ( N, M.getNumRows() );
  EXPECT_EQ( N, M.getNumColumns() );

  int expected[] = {
      1, 0, 0,
      0, 1, 0,
      0, 0, 1
  };

  for ( int j=0; j < N; ++j ) {
     const int* column = M.getColumn( j );
     for ( int i=0; i < N; ++i ) {
        EXPECT_EQ( expected[ j*N+i ], column[ i ] );
     } // END for all i
  } // END for all j

}

//------------------------------------------------------------------------------
TEST( numerics_matrix, getDiagonal )
{
  const int N=3;
  numerics::Matrix< int > M = numerics::Matrix< int >::identity( N );
  EXPECT_EQ( N, M.getNumRows() );
  EXPECT_EQ( N, M.getNumColumns() );
  EXPECT_EQ( N, M.getDiagonalSize() );

  int* diagonal = new int[ N ];
  M.getDiagonal( diagonal );

  for ( int i=0; i < N; ++i ) {
     EXPECT_EQ( 1, diagonal[ i ] );
  }

  delete [] diagonal;
}

//------------------------------------------------------------------------------
TEST( numerics_matrix, fillDiagonal )
{
  const int N=3;
  numerics::Matrix< int > M = numerics::Matrix< int >::identity( N );
  EXPECT_EQ( N, M.getNumRows() );
  EXPECT_EQ( N, M.getNumColumns() );
  EXPECT_EQ( N, M.getDiagonalSize() );

  M.fillDiagonal( 3 );

  int* diagonal = new int[ N ];
   M.getDiagonal( diagonal );

  for ( int i=0; i < N; ++i ) {
    EXPECT_EQ( 3, diagonal[ i ] );
  }

  delete [] diagonal;
}

//------------------------------------------------------------------------------
TEST( numerics_matrix, fillRow )
{
  const int FILL_VAL   = 3;
  const int TARGET_ROW = 1;

  const int N=3;
  numerics::Matrix< int > M = numerics::Matrix< int >::identity( N );
  EXPECT_EQ( N, M.getNumRows() );
  EXPECT_EQ( N, M.getNumColumns() );

  M.fillRow( TARGET_ROW, FILL_VAL );

  for ( int i=0; i < N; ++i ) {
     EXPECT_EQ( FILL_VAL, M( TARGET_ROW, i ) );
  }

}

//------------------------------------------------------------------------------
TEST( numerics_matrix, fillColumn )
{
  const int FILL_VAL   = 3;
  const int TARGET_COL = 1;

  const int N=3;
  numerics::Matrix< int > M = numerics::Matrix< int >::identity( N );
  EXPECT_EQ( N, M.getNumRows() );
  EXPECT_EQ( N, M.getNumColumns() );

  M.fillColumn( TARGET_COL, FILL_VAL );

  for ( int i=0; i < N; ++i ) {
     EXPECT_EQ( FILL_VAL, M( i, TARGET_COL ) );
  }

}

//------------------------------------------------------------------------------
TEST( numerics_matrix, fill )
{
  const int FILL_VAL = 3;

  const int N=3;
  numerics::Matrix< int > M = numerics::Matrix< int >::identity( N );
  EXPECT_EQ( N, M.getNumRows() );
  EXPECT_EQ( N, M.getNumColumns() );

  M.fill( FILL_VAL );

  for ( int i=0; i < N; ++i ) {
     for ( int j=0; j < N; ++j ) {
        EXPECT_EQ( FILL_VAL, M(i,j) );
     }
  }

}

//------------------------------------------------------------------------------
TEST( numerics_matrix, swapRows )
{
  const int M=2;
  const int N=3;
  numerics::Matrix< int > A = numerics::Matrix< int >::zeros(M,N);

  EXPECT_EQ( 2, M );
  int FILL_VAL[2] = { 3, 9 };

  A.fillRow( 0, FILL_VAL[0] );
  A.fillRow( 1, FILL_VAL[1] );

  A.swapRows( 0, 1 );

  typedef typename numerics::Matrix< int >::IndexType IndexType;

  for ( IndexType i=0; i < N; ++i ) {
     int* column = A.getColumn( i );
     for ( IndexType j=0; j < M; ++j ) {
        EXPECT_EQ( FILL_VAL[ (j+1) % M ], column[ j ] );
     }
  }

}

//------------------------------------------------------------------------------
TEST( numerics_matrix, swapColumns )
{
  const int M=3;
  const int N=4;

  typedef typename numerics::Matrix< int >::IndexType IndexType;

  // setup a test matrix
  numerics::Matrix< int > A( M,N );
  for ( IndexType i=0; i < N; ++i ) {
     A.fillColumn( i, i+1 );
  }

  int first_column = 0;
  int last_column  = N-1;
  A.swapColumns( first_column, last_column );

  int* first_column_data = A.getColumn( first_column );
  int* last_column_data  = A.getColumn( last_column );
  for ( IndexType i=0; i < M; ++i ) {
     EXPECT_EQ( last_column+1, first_column_data[ i ] );
     EXPECT_EQ( first_column+1, last_column_data[ i ] );
  }

}

//------------------------------------------------------------------------------
TEST( numerics_matrix, addition )
{
  const int N=3;
  numerics::Matrix< int > I1 = numerics::Matrix< int >::identity( N );
  numerics::Matrix< int > I2 = numerics::Matrix< int >::identity( N );

  numerics::Matrix< int > A = I1 + I2;
  const int nrows = A.getNumRows();
  const int ncols = A.getNumColumns();

  EXPECT_EQ( N, nrows );
  EXPECT_EQ( N, ncols );

  for ( int i=0; i < nrows; ++i ) {
     for ( int j=0; j < ncols; ++j ) {
        int expected = (i==j)? 2 : 0;
        EXPECT_EQ( expected, A(i,j) );
     } // END for all columns
  } // END for all rows

}

//------------------------------------------------------------------------------
TEST( numerics_matrix, subtraction )
{
  const int N=3;
  numerics::Matrix< int > I1 = numerics::Matrix< int >::identity( N );
  numerics::Matrix< int > I2 = numerics::Matrix< int >::identity( N );

  numerics::Matrix< int > A = I1 - I2;
  const int nrows = A.getNumRows();
  const int ncols = A.getNumColumns();

  EXPECT_EQ( N, nrows );
  EXPECT_EQ( N, ncols );

  for ( int i=0; i < nrows; ++i ) {
     for ( int j=0; j < ncols; ++j ) {
        int expected = 0;
        EXPECT_EQ( expected, A(i,j) );
     } // END for all columns
  } // END for all rows

}

//------------------------------------------------------------------------------
TEST( numerics_matrix, matrix_scalar_multiplication )
{
  const int N=3;
  numerics::Matrix< int > I = numerics::Matrix< int >::identity( N );

  const int scalar = 2;

  numerics::Matrix< int > A = scalar*I;
  const int nrows = A.getNumRows();
  const int ncols = A.getNumColumns();

  EXPECT_EQ( N, nrows );
  EXPECT_EQ( N, ncols );

  for ( int i=0; i < nrows; ++i ) {
     for ( int j=0; j < ncols; ++j ) {
        int expected = (i==j)? scalar : 0;
        EXPECT_EQ( expected, A(i,j) );
     }
  }
}

//------------------------------------------------------------------------------
TEST( numerics_matrix, matrix_vector_multiplication )
{
  const int N=2;
  numerics::Matrix< int > A( N,N );
  A( 0,0 ) = 1; A( 0,1 ) = 2;
  A( 1,0 ) = 3; A( 1,1 ) = 4;

  int v[2] = { 1, 2 };
  int expected[2] = { 5, 11 };

  numerics::Matrix< int > result = A*v;
  EXPECT_EQ( N, result.getNumRows() );
  EXPECT_EQ( 1, result.getNumColumns() );

  for ( int i=0; i < N; ++i ) {
     EXPECT_EQ( expected[i], result(i,0) );
  }

}

//------------------------------------------------------------------------------
TEST( numerics_matrix, matrix_matrix_multiplication )
{
  // Setup test matrix A
  numerics::Matrix< int > A( 2, 3 );
  A( 0,0 ) = 1; A( 0,1 ) = 2; A( 0,2 ) = 3;
  A( 1,0 ) = 4; A( 1,1 ) = 5; A( 1,2 ) = 6;

  // Setup test matrix B
  numerics::Matrix< int > B( 3, 2 );
  B( 0,0 ) = 1; B( 0,1 ) = 2;
  B( 1,0 ) = 3; B( 1,1 ) = 4;
  B( 2,0 ) = 5; B( 2,1 ) = 6;

  numerics::Matrix< int > C = A * B;
  EXPECT_EQ( 2, C.getNumRows() );
  EXPECT_EQ( 2, C.getNumColumns() );

  EXPECT_EQ( 22, C(0,0) );
  EXPECT_EQ( 28, C(0,1) );
  EXPECT_EQ( 49, C(1,0) );
  EXPECT_EQ( 64, C(1,1) );

  numerics::Matrix< int > I2 = numerics::Matrix< int >::identity( 2 );
  numerics::Matrix< int > C2 = C*I2;
  EXPECT_EQ( 22, C2(0,0) );
  EXPECT_EQ( 28, C2(0,1) );
  EXPECT_EQ( 49, C2(1,0) );
  EXPECT_EQ( 64, C2(1,1) );
}

//------------------------------------------------------------------------------
TEST( numerics_matrix, transpose )
{
  // Setup matrix: (1 2 3)
  //               (4 5 6)
  numerics::Matrix< int > A( 2, 3 );
  A( 0,0 ) = 1; A( 0,1 ) = 2; A( 0,2 ) = 3;
  A( 1,0 ) = 4; A( 1,1 ) = 5; A( 1,2 ) = 6;

  // Matrix transpose
  numerics::Matrix< int > M = numerics::transpose( A );
  EXPECT_EQ( 3, M.getNumRows() );
  EXPECT_EQ( 2, M.getNumColumns() );

  EXPECT_EQ( 1, M(0,0) );
  EXPECT_EQ( 2, M(1,0) );
  EXPECT_EQ( 3, M(2,0) );

  EXPECT_EQ( 4, M(0,1) );
  EXPECT_EQ( 5, M(1,1) );
  EXPECT_EQ( 6, M(2,1) );
}

//------------------------------------------------------------------------------
TEST( numerics_matrix, output_stream )
{
  numerics::Matrix< int > A( 2, 3 );
  A( 0,0 ) = 1; A( 0,1 ) = 2; A( 0,2 ) = 3;
  A( 1,0 ) = 4; A( 1,1 ) = 5; A( 1,2 ) = 6;

  std::ostringstream expected_stream;
  expected_stream << "[ 1 2 3 ]\n" << "[ 4 5 6 ]\n";

  std::ostringstream actual_stream;
  actual_stream << A;

  EXPECT_EQ( expected_stream.str(), actual_stream.str() );
}

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
#include "axom_utils/Determinants.hpp"

namespace numerics = axom::numerics;

TEST( numerics_determinants, determinant_of_In )
{
  const int N = 25;

  for ( int i=2 ; i < N ; ++i )
  {
    numerics::Matrix< double > In = numerics::Matrix< double >::identity( i );
    double det = numerics::determinant( In );
    EXPECT_DOUBLE_EQ( 1.0, det );
  }

}

//------------------------------------------------------------------------------
TEST( numerics_determinants, determinant5x5 )
{

  const int N = 5;
  numerics::Matrix< double > A( N,N );
  A( 0,0 )=1; A( 0,1 )=2; A( 0,2 )=4; A( 0,3 )=3; A( 0,4 )=0;
  A( 1,0 )=2; A( 1,1 )=1; A( 1,2 )=-1; A( 1,3 )=1; A( 1,4 )=3;
  A( 2,0 )=4; A( 2,1 )=-1; A( 2,2 )=-2; A( 2,3 )=5; A( 2,4 )=1;
  A( 3,0 )=7; A( 3,1 )=3; A( 3,2 )=6; A( 3,3 )=2; A( 3,4 )=1;
  A( 4,0 )=1; A( 4,1 )=0; A( 4,2 )=-1; A( 4,3 )=1; A( 4,4 )=1;

  double computed_det = numerics::determinant( A );
  double expected_det = -34.0;
  EXPECT_DOUBLE_EQ( expected_det, computed_det );
}

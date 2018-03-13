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

#include "mint/Extent.hpp"
using namespace axom;

TEST( mint_extent, basic )
{
  int ext[6]={ 0,5,0,5,0,5 };
  mint::Extent< int > extent(3,ext);

  EXPECT_EQ( 3, extent.getDimension() );

  int numNodes = 1;
  int numCells = 1;
  for ( int i=0 ; i < 3 ; ++i )
  {
    EXPECT_EQ(  extent.min( i ),  0  );
    EXPECT_EQ(  extent.max( i ),  5  );
    EXPECT_EQ(  extent.size( i ), 6 );

    numNodes *= extent.size( i );
    numCells *= (extent.size( i )-1);
  }

  EXPECT_EQ(  numNodes, extent.getNumNodes() );
  EXPECT_EQ(  numCells, extent.getNumCells() );

  const int imin = extent.min(0);
  const int imax = extent.max(0);
  const int jmin = extent.min(1);
  const int jmax = extent.max(1);
  const int kmin = extent.min(2);
  const int kmax = extent.max(2);

  int count = 0;
  for ( int k=kmin ; k <= kmax ; ++k )
  {
    for ( int j=jmin ; j <= jmax ; ++j )
    {
      for ( int i=imin ; i <= imax ; ++i )
      {

        const int idx = extent.getLinearIndex( i,j,k );
        EXPECT_EQ( count, idx );

        int ii = -1;
        int jj = -1;
        int kk = -1;
        extent.getGridIndex( idx, ii, jj, kk );

        EXPECT_EQ(i,ii);
        EXPECT_EQ(j,jj);
        EXPECT_EQ(k,kk);

        ++count;

      }    // END for all i
    }   // END for all j
  }  // END for all k

}

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
#include "mint/config.hpp"
using namespace axom::mint;

TEST( mint_mesh_extent, basic )
{
  int64 ext[6]={ 0,5,0,5,0,5 };
  Extent extent(3,ext);

  EXPECT_EQ( 3, extent.getDimension() );

  IndexType numNodes = 1;
  IndexType numCells = 1;
  for ( int i=0 ; i < 3 ; ++i )
  {
    EXPECT_EQ( extent.min( i ),  0  );
    EXPECT_EQ( extent.max( i ),  5  );
    EXPECT_EQ( extent.size( i ), 6 );

    numNodes *= extent.size( i );
    numCells *= (extent.size( i )-1);
  }

  EXPECT_EQ( numNodes, extent.getNumNodes() );
  EXPECT_EQ( numCells, extent.getNumCells() );

  const int64 imin = extent.min(0);
  const int64 imax = extent.max(0);
  const int64 jmin = extent.min(1);
  const int64 jmax = extent.max(1);
  const int64 kmin = extent.min(2);
  const int64 kmax = extent.max(2);

  IndexType count = 0;
  for ( int64 k=kmin ; k <= kmax ; ++k )
  {
    for ( int64 j=jmin ; j <= jmax ; ++j )
    {
      for ( int64 i=imin ; i <= imax ; ++i )
      {

        const IndexType idx = extent.getLinearIndex( i,j,k );
        EXPECT_EQ( count, idx );

        IndexType ii = -1;
        IndexType jj = -1;
        IndexType kk = -1;
        extent.getGridIndex( idx, ii, jj, kk );

        EXPECT_EQ(i,ii);
        EXPECT_EQ(j,jj);
        EXPECT_EQ(k,kk);

        ++count;

      }    // END for all i
    }   // END for all j
  }  // END for all k

}

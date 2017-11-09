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
#include "mint/DataTypes.hpp"
using namespace axom::mint;

TEST( mint_extent, basic )
{
  globalIndex ext[6]={ 0,5,0,5,0,5 };
  Extent extent(3,ext);

  EXPECT_EQ( 3, extent.getDimension() );

  localIndex numNodes = 1;
  localIndex numCells = 1;
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

  const globalIndex imin = extent.min(0);
  const globalIndex imax = extent.max(0);
  const globalIndex jmin = extent.min(1);
  const globalIndex jmax = extent.max(1);
  const globalIndex kmin = extent.min(2);
  const globalIndex kmax = extent.max(2);

  localIndex count = 0;
  for ( globalIndex k=kmin ; k <= kmax ; ++k )
  {
    for ( globalIndex j=jmin ; j <= jmax ; ++j )
    {
      for ( globalIndex i=imin ; i <= imax ; ++i )
      {

        const localIndex idx = extent.getLinearIndex( i,j,k );
        EXPECT_EQ( count, idx );

        localIndex ii = -1;
        localIndex jj = -1;
        localIndex kk = -1;
        extent.getGridIndex( idx, ii, jj, kk );

        EXPECT_EQ(i,ii);
        EXPECT_EQ(j,jj);
        EXPECT_EQ(k,kk);

        ++count;

      }    // END for all i
    }   // END for all j
  }  // END for all k

}

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

#include "axom/mint/config.hpp"
#include "axom/mint/mesh/Extent.hpp"

// aliases
namespace mint    = axom::mint;
using int64       = mint::int64;
using IndexType   = mint::IndexType;


//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------
namespace
{

void check_constructor( int ndims, const int64* ext )
{
  EXPECT_TRUE( (ndims >= 1) && (ndims <= 3) );
  EXPECT_TRUE( (ext != nullptr) );

  mint::Extent extent( ndims, ext );
  EXPECT_EQ( extent.getDimension(), ndims );

  int64 expected_num_nodes = 1;
  int64 expected_num_cells = 1;
  int64 expected_jp = 0;
  int64 expected_kp = 0;

  for ( int i=0 ; i < ndims ; ++i )
  {
    const int64 ilo = ext[ i*2   ];
    const int64 ihi = ext[ i*2+1 ];
    EXPECT_EQ( extent.min( i ), ilo );
    EXPECT_EQ( extent.max( i ), ihi );

    const int64 expected_size = ihi-ilo+1;
    EXPECT_EQ( extent.size( i ), expected_size );

    expected_num_nodes *= expected_size;
    expected_num_cells *= ( expected_size-1 );

    if ( i== 1 )
    {
      expected_jp = expected_size;
    }

    if ( i==2 )
    {
      expected_kp =  expected_jp * expected_size;
    }

  } // END for all dimensions

  EXPECT_EQ( extent.getNumNodes(), expected_num_nodes );
  EXPECT_EQ( extent.getNumCells(), expected_num_cells );
  EXPECT_EQ( extent.jp(), expected_jp );
  EXPECT_EQ( extent.kp(), expected_kp );
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------

TEST( mint_mesh_extent, constructor )
{
  constexpr int NDIMS = 3;
  int64 extent1[]     = {  0,4,  0,4,  0,4 };
  int64 extent2[]     = { -2,2, -2,2, -2,2 };

  for ( int idim=1 ; idim <= NDIMS ; ++idim )
  {
    check_constructor( idim, extent1 );
    check_constructor( idim, extent2 );

    mint::Extent ext1( idim, extent1 );
    mint::Extent ext2( idim, extent2 );
    EXPECT_EQ( ext1.getNumNodes(), ext2.getNumNodes() );
    EXPECT_EQ( ext1.getNumCells(), ext2.getNumCells() );

  }

}

//------------------------------------------------------------------------------
TEST( mint_mesh_extent, local_to_global_bijective_mapping )
{
  constexpr int NDIMS = 3;
  int64 extent[]      = { -2,2, -2,2, -2,2 };

  for ( int idim=1 ; idim <= NDIMS ; ++idim )
  {
    mint::Extent ext( idim, extent );

    switch ( idim )
    {
    case 1:
    {

      const int64 imin = ext.min( mint::I_DIRECTION );
      const int64 imax = ext.max( mint::I_DIRECTION );
      for ( int64 i=imin ; i <= imax ; ++i )
      {
        int64 gijk[ 1 ]  = { i  };
        IndexType lijk[ 1 ]  = { -1 };
        int64 comp_gijk[ 1 ] = { -1 };

        ext.shiftToLocal( gijk, lijk );
        EXPECT_TRUE( lijk[ 0 ] >= 0 && lijk[ 0 ] < ext.size( 0 ) );

        ext.shiftToGlobal( lijk, comp_gijk );
        EXPECT_EQ( gijk[ 0 ], comp_gijk[ 0 ] );

      } // END for all i

    }   // END 1-D
    break;
    case 2:
    {

      const int64 imin = ext.min( mint::I_DIRECTION );
      const int64 imax = ext.max( mint::I_DIRECTION );
      const int64 jmin = ext.min( mint::J_DIRECTION );
      const int64 jmax = ext.max( mint::J_DIRECTION );

      for ( int64 j=jmin ; j <= jmax ; ++j )
      {
        for ( int64 i=imin ; i <= imax ; ++i )
        {
          int64 gijk[ 2 ]  = {  i,  j };
          IndexType lijk[ 2 ]  = { -1, -1 };
          int64 comp_gijk[ 2 ] = { -1, -1 };

          ext.shiftToLocal( gijk, lijk );
          EXPECT_TRUE( lijk[ 0 ] >= 0 && lijk[ 0 ] < ext.size( 0 ) );
          EXPECT_TRUE( lijk[ 1 ] >= 0 && lijk[ 1 ] < ext.size( 1 ) );

          ext.shiftToGlobal( lijk, comp_gijk );
          EXPECT_EQ( gijk[ 0 ], comp_gijk[ 0 ] );
          EXPECT_EQ( gijk[ 1 ], comp_gijk[ 1 ] );

        } // END for all i
      } // END for all j

    }   // END 2-D
    break;
    default:
      EXPECT_EQ( idim, 3 );
      {

        const int64 imin = ext.min( mint::I_DIRECTION );
        const int64 imax = ext.max( mint::I_DIRECTION );
        const int64 jmin = ext.min( mint::J_DIRECTION );
        const int64 jmax = ext.max( mint::J_DIRECTION );
        const int64 kmin = ext.min( mint::K_DIRECTION );
        const int64 kmax = ext.max( mint::K_DIRECTION );

        for ( int64 k=kmin ; k <= kmax ; ++k )
        {
          for ( int64 j=jmin ; j <= jmax ; ++j )
          {
            for ( int64 i=imin ; i <= imax ; ++i )
            {
              int64 gijk[ 3 ]  = {  i,  j,  k };
              IndexType lijk[ 3 ]  = { -1, -1, -1 };
              int64 comp_gijk[ 3 ] = { -1, -1, -1 };

              ext.shiftToLocal( gijk, lijk );
              EXPECT_TRUE( lijk[ 0 ] >= 0 && lijk[ 0 ] < ext.size( 0 ) );
              EXPECT_TRUE( lijk[ 1 ] >= 0 && lijk[ 1 ] < ext.size( 1 ) );
              EXPECT_TRUE( lijk[ 2 ] >= 0 && lijk[ 2 ] < ext.size( 2 ) );

              ext.shiftToGlobal( lijk, comp_gijk );
              EXPECT_EQ( gijk[ 0 ], comp_gijk[ 0 ] );
              EXPECT_EQ( gijk[ 1 ], comp_gijk[ 1 ] );
              EXPECT_EQ( gijk[ 2 ], comp_gijk[ 2 ] );

            } // END for all i
          } // END for all j
        } // END for all k

      } // END 3-D

    } // END switch

  } // END for all NDIMS
}

//------------------------------------------------------------------------------
TEST( mint_mesh_extent, nodal_linear_to_ijk_bijective_map )
{
  constexpr int NDIMS = 3;
  int64 extent[]      = {  0,4,  0,4,  0,4 };

  for ( int idim=2 ; idim <= NDIMS ; ++idim )
  {
    mint::Extent ext( idim, extent );

    switch ( idim )
    {
    case 2:
    {

      const IndexType Ni = ext.size( mint::I_DIRECTION );
      const IndexType Nj = ext.size( mint::J_DIRECTION );

      for ( IndexType j=0 ; j < Nj ; ++j )
      {
        for ( IndexType i=0 ; i < Ni ; ++i )
        {
          const IndexType idx = ext.getLinearIndex( i, j );
          EXPECT_TRUE( (idx >= 0) && (idx < ext.getNumNodes()) );

          IndexType ii = -1;
          IndexType jj = -1;
          ext.getGridIndex( idx, ii, jj );
          EXPECT_EQ( ii, i );
          EXPECT_EQ( jj, j );

        } // END for all i
      } // END for all j

    }   // END 2-D
    break;
    default:
      EXPECT_EQ( idim, 3 );
      {

        const IndexType Ni = ext.size( mint::I_DIRECTION );
        const IndexType Nj = ext.size( mint::J_DIRECTION );
        const IndexType Nk = ext.size( mint::K_DIRECTION );

        for ( IndexType k=0 ; k < Nk ; ++k )
        {
          for ( IndexType j=0 ; j < Nj ; ++j )
          {
            for ( IndexType i=0 ; i < Ni ; ++i )
            {
              const IndexType idx = ext.getLinearIndex( i, j, k );
              EXPECT_TRUE( (idx >=0) && (idx < ext.getNumNodes() ) );

              IndexType ii = -1;
              IndexType jj = -1;
              IndexType kk = -1;
              ext.getGridIndex( idx, ii, jj, kk );
              EXPECT_EQ( ii, i );
              EXPECT_EQ( jj, j );
              EXPECT_EQ( kk, k );
            } // END for all i
          } // END for all j
        } // END for all k

      } // END 3-D
    } // END switch

  } // END for all NDIMS

}

//------------------------------------------------------------------------------
TEST( mint_mesh_extent, cell_linear_to_ijk_bijective_map )
{

  constexpr int NDIMS = 3;
  int64 extent[]      = {  0,4,  0,4,  0,4 };

  for ( int idim=2 ; idim <= NDIMS ; ++idim )
  {
    mint::Extent ext( idim, extent );

    switch ( idim )
    {
    case 2:
    {

      const IndexType Ni = ext.size( mint::I_DIRECTION ) - 1;
      const IndexType Nj = ext.size( mint::J_DIRECTION ) - 1;

      for ( IndexType j=0 ; j < Nj ; ++j )
      {
        for ( IndexType i=0 ; i < Ni ; ++i )
        {
          const IndexType idx = ext.getCellLinearIndex( i, j );
          EXPECT_TRUE( (idx >= 0) && (idx < ext.getNumNodes()) );

          IndexType ii = -1;
          IndexType jj = -1;
          ext.getCellGridIndex( idx, ii, jj );
          EXPECT_EQ( ii, i );
          EXPECT_EQ( jj, j );

        } // END for all i
      } // END for all j

    } // END 2-D
    break;
    default:
      EXPECT_EQ( idim, 3 );
      {

        const IndexType Ni = ext.size( mint::I_DIRECTION ) - 1;
        const IndexType Nj = ext.size( mint::J_DIRECTION ) - 1;
        const IndexType Nk = ext.size( mint::K_DIRECTION ) - 1;

        for ( IndexType k=0 ; k < Nk ; ++k )
        {
          for ( IndexType j=0 ; j < Nj ; ++j )
          {
            for ( IndexType i=0 ; i < Ni ; ++i )
            {
              const IndexType idx = ext.getCellLinearIndex( i, j, k );
              EXPECT_TRUE( (idx >=0) && (idx < ext.getNumCells() ) );

              IndexType ii = -1;
              IndexType jj = -1;
              IndexType kk = -1;
              ext.getCellGridIndex( idx, ii, jj, kk );
              EXPECT_EQ( ii, i );
              EXPECT_EQ( jj, j );
              EXPECT_EQ( kk, k );
            } // END for all i
          } // END for all j
        } // END for all k

      } // END 3-D
    } // END switch

  } // END for all NDIMS

}

//------------------------------------------------------------------------------
TEST( mint_mesh_DeathTest, invalid_construction )
{
  const char* IGNORE_OUTPUT = ".*";
  int64 extent[]      = {  0,4,  0,4,  0,4 };

  EXPECT_DEATH_IF_SUPPORTED( mint::Extent(0,extent), IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( mint::Extent(9,extent), IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( mint::Extent(3,nullptr), IGNORE_OUTPUT );
}

//------------------------------------------------------------------------------
#include "axom/slic/core/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}

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
#include "axom/mint/config.hpp"               // for compile-time definitions

// Mint includes
#include "axom/mint/mesh/blueprint.hpp"       // for blueprint functions
#include "axom/mint/mesh/CellTypes.hpp"       // for CellType enum
#include "axom/mint/mesh/RectilinearMesh.hpp" // for RectilinearMesh
#include "axom/mint/mesh/ParticleMesh.hpp"    // for ParticleMesh
#include "StructuredMesh_helpers.hpp"         // for StructuredMesh test helpers

// Slic includes
#include "axom/slic/interface/slic.hpp"       // for slic macros

// Sidre includes
#ifdef MINT_USE_SIDRE
#include "axom/sidre/core/sidre.hpp"
namespace sidre = axom::sidre;
#endif

#include "gtest/gtest.h"                      // for gtest macros

// C/C++ includes
#include <cmath>                              // for exp
using namespace axom::mint;

// globals
const char* IGNORE_OUTPUT = ".*";

//------------------------------------------------------------------------------
//  HELPER METHODS
//------------------------------------------------------------------------------
namespace
{

void exponential_distribution( double origin, IndexType N, double* x )
{
  EXPECT_TRUE( x != nullptr );

  constexpr double beta = 0.05;
  const double expbeta  = exp( beta );
  const double invf     = 1 / ( expbeta - 1.0 );

  x[ 0 ] = origin;
  for ( int i=1 ; i < N ; ++i )
  {
    const double prev = x[ i-1 ];
    const double dx   = ( exp( i*beta ) - 1.0 ) * invf;
    x[ i ] = prev + dx;
  }

}

//------------------------------------------------------------------------------
void check_coordinate( const double* x, const double* expected, IndexType N )
{
  EXPECT_TRUE( x != nullptr );
  EXPECT_TRUE( expected != nullptr );
  EXPECT_TRUE( N >= 0 );

  for ( IndexType i=0 ; i < N ; ++i )
  {
    EXPECT_DOUBLE_EQ( x[ i ], expected[ i ] );
  }
}

//------------------------------------------------------------------------------
void check_fill_coords( RectilinearMesh* m )
{
  EXPECT_TRUE( m != nullptr );

  const int ndims = m->getDimension();
  for ( int idim=0 ; idim < ndims ; ++idim )
  {
    const IndexType N = m->getNodeExtent( idim );
    double* x = m->getCoordinateArray( idim );
    exponential_distribution( 42.0, N, x );
  }
}

//------------------------------------------------------------------------------

} // END anonymous namespace

//------------------------------------------------------------------------------
//  UNIT TESTS
//------------------------------------------------------------------------------
TEST( mint_mesh_rectilinear_mesh_DeathTest, invalid_construction )
{
  const IndexType N[]  = {   5,   5,   5 };
  double x[ 5 ];

  // check 1st native constructor
  EXPECT_DEATH_IF_SUPPORTED( RectilinearMesh( 42, N ), IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( RectilinearMesh(2,AXOM_NULLPTR), IGNORE_OUTPUT );

  // check 2nd native constructor
  EXPECT_DEATH_IF_SUPPORTED( RectilinearMesh( -1,N[1],N[2] ), IGNORE_OUTPUT );

  // check external constructor
  EXPECT_DEATH_IF_SUPPORTED( RectilinearMesh(AXOM_NULLPTR,x), IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( RectilinearMesh(N,AXOM_NULLPTR), IGNORE_OUTPUT );

#ifdef MINT_USE_SIDRE

  sidre::DataStore ds;
  sidre::Group* root          = ds.getRoot();
  sidre::Group* valid_group   = root->createGroup( "mesh" );
  sidre::Group* particle_mesh = root->createGroup( "particle_mesh" );
  ParticleMesh( 3, 10, particle_mesh );

  // check pull constructor
  EXPECT_DEATH_IF_SUPPORTED( RectilinearMesh(nullptr,""), IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( RectilinearMesh(root,""), IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( RectilinearMesh(particle_mesh,""), IGNORE_OUTPUT );

  // check 1st push constructor
  EXPECT_DEATH_IF_SUPPORTED( RectilinearMesh(42, N, valid_group),
                             IGNORE_OUTPUT );
  EXPECT_EQ( valid_group->getNumGroups(), 0 );
  EXPECT_EQ( valid_group->getNumViews(), 0 );

  EXPECT_DEATH_IF_SUPPORTED( RectilinearMesh(3,nullptr,valid_group),
                             IGNORE_OUTPUT );
  EXPECT_EQ( valid_group->getNumGroups(), 0 );
  EXPECT_EQ( valid_group->getNumViews(), 0 );

  EXPECT_DEATH_IF_SUPPORTED( RectilinearMesh(3,N,AXOM_NULLPTR),
                             IGNORE_OUTPUT );
  EXPECT_EQ( valid_group->getNumGroups(), 0 );
  EXPECT_EQ( valid_group->getNumViews(), 0 );

  // check 2nd push constructor
  EXPECT_DEATH_IF_SUPPORTED( RectilinearMesh( nullptr, N[0] ),
                             IGNORE_OUTPUT );

  EXPECT_DEATH_IF_SUPPORTED( RectilinearMesh( valid_group, -1 ),
                             IGNORE_OUTPUT );

#endif

}

//------------------------------------------------------------------------------
TEST( mint_mesh_rectilinear_mesh, native_constructor )
{
  constexpr int NDIMS    = 3;
  const IndexType N[]    = {   5,   5,   5 };

  for ( int idim=1 ; idim <= NDIMS ; ++idim )
  {
    RectilinearMesh* m = new RectilinearMesh( idim, N );
    internal::check_constructor( m, STRUCTURED_RECTILINEAR_MESH, idim, N );
    EXPECT_FALSE( m->isExternal() );
    EXPECT_FALSE( m->hasSidreGroup() );
    check_fill_coords( m );
    internal::check_create_fields( m );
    delete m;

    switch ( idim )
    {
    case 1:
      m = new RectilinearMesh( N[ 0 ] );
      break;
    case 2:
      m = new RectilinearMesh( N[ 0 ], N[ 1 ] );
      break;
    default:
      EXPECT_EQ( idim, 3 );
      m = new RectilinearMesh( N[ 0 ], N[ 1 ], N[ 2 ] );
    }

    internal::check_constructor( m, STRUCTURED_RECTILINEAR_MESH, idim, N );
    check_fill_coords( m );
    internal::check_create_fields( m );

    delete m;
  } // END for all dimensions

}

//------------------------------------------------------------------------------
TEST( mint_mesh_rectilinear_mesh, external_costructor )
{
  constexpr int NDIMS        = 3;
  const IndexType N[]        = {   5,   5,   5 };
  constexpr double MAGIC_VAL = 42.0;

  double* X = new double[ N[ 0 ] ];
  double* Y = new double[ N[ 1 ] ];
  double* Z = new double[ N[ 2 ] ];
  exponential_distribution( MAGIC_VAL, N[ 0 ], X );
  exponential_distribution( MAGIC_VAL, N[ 1 ], Y );
  exponential_distribution( MAGIC_VAL, N[ 2 ], Z );

  double* x = new double[ N[ 0 ] ];
  double* y = new double[ N[ 1 ] ];
  double* z = new double[ N[ 2 ] ];
  exponential_distribution( MAGIC_VAL, N[ 0 ], x );
  exponential_distribution( MAGIC_VAL, N[ 1 ], y );
  exponential_distribution( MAGIC_VAL, N[ 2 ], z );

  for ( int idim=1 ; idim <= NDIMS ; ++idim )
  {
    RectilinearMesh* m = nullptr;

    switch ( idim )
    {
    case 1:
      m = new RectilinearMesh( N, x );
      EXPECT_EQ( m->getCoordinateArray(X_COORDINATE), x );
      check_coordinate( m->getCoordinateArray(X_COORDINATE), x, N[0] );
      break;
    case 2:
      m = new RectilinearMesh( N, x, y );
      EXPECT_EQ( m->getCoordinateArray(X_COORDINATE), x );
      EXPECT_EQ( m->getCoordinateArray(Y_COORDINATE), y );
      check_coordinate( m->getCoordinateArray(X_COORDINATE), x, N[0] );
      check_coordinate( m->getCoordinateArray(Y_COORDINATE), y, N[1] );
      break;
    default:
      m = new RectilinearMesh( N, x, y, z );
      EXPECT_EQ( m->getCoordinateArray(X_COORDINATE), x );
      EXPECT_EQ( m->getCoordinateArray(Y_COORDINATE), y );
      EXPECT_EQ( m->getCoordinateArray(Z_COORDINATE), z );
      check_coordinate( m->getCoordinateArray(X_COORDINATE), x, N[0] );
      check_coordinate( m->getCoordinateArray(Y_COORDINATE), y, N[1] );
      check_coordinate( m->getCoordinateArray(Z_COORDINATE), z, N[2] );

    } // END switch

    EXPECT_TRUE( m != AXOM_NULLPTR );
    internal::check_constructor( m, STRUCTURED_RECTILINEAR_MESH, idim, N );
    EXPECT_FALSE( m->hasSidreGroup() );
    EXPECT_TRUE( m->isExternal() );

    // deallocate
    delete m;
    m = nullptr;

    // ensure coordinates are not changed after mesh gets deleted
    EXPECT_TRUE( x != nullptr );
    EXPECT_TRUE( y != nullptr );
    EXPECT_TRUE( z != nullptr );
    check_coordinate( x, X, N[ 0 ] );
    check_coordinate( y, Y, N[ 1 ] );
    check_coordinate( z, Z, N[ 2 ] );

  } // END for all dimensions

  delete [] x;
  delete [] X;
  delete [] y;
  delete [] Y;
  delete [] z;
  delete [] Z;
}

//------------------------------------------------------------------------------
#ifdef MINT_USE_SIDRE

TEST( mint_mesh_rectilinear_mesh, sidre_constructor )
{
  constexpr int NDIMS        = 3;
  const IndexType N[]        = {   5,   5,   5 };
  constexpr double MAGIC_VAL = 42.0;

  double* expected_coords = new double[ N[ 0 ] ];
  exponential_distribution( MAGIC_VAL, N[0], expected_coords );

  for ( int idim=1 ; idim <= NDIMS ; ++idim )
  {
    // STEP 0: create a data-store with two groups
    sidre::DataStore ds;
    sidre::Group* root  = ds.getRoot( );
    sidre::Group* m1grp = root->createGroup( "m1" );
    sidre::Group* m2grp = root->createGroup( "m2" );

    // STEP 1: populate meshes in Sidre using the 2 sidre constructors
    // BEGIN SCOPE
    {
      RectilinearMesh* m1 = new RectilinearMesh( idim, N, m1grp );
      EXPECT_TRUE( m1->hasSidreGroup() );
      EXPECT_FALSE( m1->isExternal() );
      internal::check_constructor( m1, STRUCTURED_RECTILINEAR_MESH, idim, N );
      check_fill_coords( m1 );
      internal::check_create_fields( m1 );

      delete m1;
      m1 = nullptr;

      RectilinearMesh* m2 = nullptr;
      switch ( idim )
      {
      case 1:
        m2 = new RectilinearMesh( m2grp, N[ I_DIRECTION ] );
        break;
      case 2:
        m2 = new RectilinearMesh( m2grp, N[ I_DIRECTION ], N[ J_DIRECTION ] );
        break;
      default:
        EXPECT_EQ( idim, 3 );
        m2 = new RectilinearMesh( m2grp, N[ I_DIRECTION ],
                                  N[ J_DIRECTION ],
                                  N[ K_DIRECTION ]  );
      } // END switch

      internal::check_constructor( m2, STRUCTURED_RECTILINEAR_MESH, idim, N );
      EXPECT_TRUE( m2->hasSidreGroup() );
      EXPECT_FALSE( m2->isExternal() );check_fill_coords( m2 );
      internal::check_create_fields( m2 );

      delete m2;
      m2 = nullptr;
    }
    // END SCOPE

    // STEP 2: pull the 2 meshes from sidre and check correctness
    // BEGIN SCOPE
    {
      // check m1
      RectilinearMesh* m = new RectilinearMesh( m1grp );
      internal::check_constructor( m, STRUCTURED_RECTILINEAR_MESH, idim, N );
      EXPECT_TRUE( m->hasSidreGroup() );
      EXPECT_FALSE( m->isExternal() );

      for ( int ii=0 ; ii < idim ; ++ii )
      {
        check_coordinate( m->getCoordinateArray( ii ), expected_coords,
                          m->getNodeExtent( ii ) );
      }

      internal::check_fields( m, true );
      delete m;

      // check m2
      m = new RectilinearMesh( m2grp );
      internal::check_constructor( m, STRUCTURED_RECTILINEAR_MESH, idim, N );
      EXPECT_TRUE( m->hasSidreGroup() );
      EXPECT_FALSE( m->isExternal() );

      for ( int ii=0 ; ii < idim ; ++ii )
      {
        check_coordinate( m->getCoordinateArray( ii ), expected_coords,
                          m->getNodeExtent( ii ) );
      }

      internal::check_fields( m, true );
      delete m;
    }
    // ENDE SCOPE

    // STEP 3: ensure data is persistent in Sidre
    EXPECT_TRUE( blueprint::isValidRootGroup( m1grp ) );
    EXPECT_TRUE( blueprint::isValidRootGroup( m2grp ) );

  } // END for all dimensions

  delete [] expected_coords;
  expected_coords = nullptr;
}

#endif /* ENDIF MINT_USE_SIDRE */

//------------------------------------------------------------------------------
TEST( mint_mesh_rectilinear_mesh, get_node )
{
  // initialize coordinates
  const IndexType N[] = { 5,    5,    5 };
  double x[ 5 ];
  double y[ 5 ];
  double z[ 5 ];
  exponential_distribution( 42.0, 5, x );
  exponential_distribution( 0.0, 5, y );
  exponential_distribution( 0.0, 5, z );

  RectilinearMesh m( N, x, y, z );
  internal::check_constructor( &m, STRUCTURED_RECTILINEAR_MESH, 3, N );

  const IndexType kp = m.nodeKp( );
  const IndexType jp = m.nodeJp( );
  const double* xx   = m.getCoordinateArray( X_COORDINATE );
  const double* yy   = m.getCoordinateArray( Y_COORDINATE );
  const double* zz   = m.getCoordinateArray( Z_COORDINATE );
  EXPECT_EQ( x, xx );
  EXPECT_EQ( y, yy );
  EXPECT_EQ( z, zz );

  for ( IndexType k=0 ; k < N[ K_DIRECTION ] ; ++k )
  {
    const IndexType k_offset = k*kp;
    for ( IndexType j=0 ; j < N[ J_DIRECTION ] ; ++j )
    {
      const IndexType j_offset = j*jp;
      for ( IndexType i=0 ; i < N[ I_DIRECTION ] ; ++i )
      {
        const IndexType nodeID = i + j_offset + k_offset;

        double node[ 3 ];
        m.getNode( nodeID, node );

        EXPECT_DOUBLE_EQ( xx[ i ], node[ X_COORDINATE ] );
        EXPECT_DOUBLE_EQ( yy[ j ], node[ Y_COORDINATE ] );
        EXPECT_DOUBLE_EQ( zz[ k ], node[ Z_COORDINATE ] );

      }  // END for all i
    }  // END for all j
  }  // END for all k

}

//------------------------------------------------------------------------------
#include "axom/slic/core/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

int main ( int argc, char* argv[] )
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,

  // finalized when exiting main scope
  result = RUN_ALL_TESTS();
  return result;
}

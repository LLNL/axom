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
#include "axom/mint/config.hpp"          // for compile-time type definitions

// Mint includes
#include "axom/mint/mesh/blueprint.hpp"       // for blueprint functions
#include "axom/mint/mesh/CellTypes.hpp"       // for CellTypes enum definition
#include "axom/mint/mesh/CurvilinearMesh.hpp" // for CurivilinearMesh
#include "axom/mint/mesh/ParticleMesh.hpp"    // for ParticleMesh
#include "axom/mint/mesh/internal/MeshHelpers.hpp" // for internal::dim
#include "StructuredMesh_helpers.hpp" // for StructuredMesh test helpers

// Slic includes
#include "axom/slic/interface/slic.hpp"            // for slic macros

// Sidre includes
#ifdef AXOM_MINT_USE_SIDRE
#include "axom/sidre/core/sidre.hpp"
namespace sidre = axom::sidre;
#endif

#include "gtest/gtest.h"           // for gtest macros

using namespace axom::mint;

// globals
const char* IGNORE_OUTPUT = ".*";

//------------------------------------------------------------------------------
//  HELPER METHODS
//------------------------------------------------------------------------------
namespace
{

//------------------------------------------------------------------------------
void set_coordinates( IndexType N, double* x, double* y=nullptr,
                      double* z=nullptr )
{
  const IndexType ndims = internal::dim( x, y, z );
  double* const coords[ 3 ] = { x, y, z };

  double factor = internal::PI;
  for ( int dim = 0 ; dim < ndims ; ++dim )
  {
    SLIC_ASSERT( coords[ dim ] != nullptr );
    for ( IndexType i = 0 ; i < N ; ++i )
    {
      coords[ dim ][ i ] = factor * i;
    }

    factor *= internal::PI;
  }
}

//------------------------------------------------------------------------------
void set_coordinates( CurvilinearMesh* m )
{
  const int ndims = m->getDimension();
  const IndexType numNodes = m->getNumberOfNodes();
  double* x = m->getCoordinateArray( X_COORDINATE );

  if ( ndims == 1 )
  {
    set_coordinates( numNodes, x );
  }
  else if ( ndims == 2 )
  {
    double* y = m->getCoordinateArray( Y_COORDINATE );
    set_coordinates( numNodes, x, y );
  }
  else
  {
    double* y = m->getCoordinateArray( Y_COORDINATE );
    double* z = m->getCoordinateArray( Z_COORDINATE );
    set_coordinates( numNodes, x, y, z );
  }
}

//------------------------------------------------------------------------------
void check_coordinates( IndexType N, int ndims, const double* x,
                        const double* y=nullptr,
                        const double* z=nullptr  )
{
  const double* const coords[ 3 ] = { x, y, z };

  double factor = internal::PI;
  for ( int dim = 0 ; dim < ndims ; ++dim )
  {
    ASSERT_NE( coords[ dim ], nullptr );
    for ( IndexType i = 0 ; i < N ; ++i )
    {
      EXPECT_DOUBLE_EQ( coords[ dim ][ i ], factor * i );
    }

    factor *= internal::PI;
  }
}

//------------------------------------------------------------------------------
void check_coordinates( CurvilinearMesh* m )
{
  EXPECT_TRUE( m != nullptr );
  const int ndims = m->getDimension();
  const IndexType numNodes = m->getNumberOfNodes();
  double* x = m->getCoordinateArray( X_COORDINATE );

  if ( ndims == 1 )
  {
    check_coordinates( numNodes, ndims, x );
  }
  else if ( ndims == 2 )
  {
    double* y = m->getCoordinateArray( Y_COORDINATE );
    check_coordinates( numNodes, ndims, x, y );
  }
  else
  {
    double* y = m->getCoordinateArray( Y_COORDINATE );
    double* z = m->getCoordinateArray( Z_COORDINATE );
    check_coordinates( numNodes, ndims, x, y, z );
  }
}

} // END namespace

//------------------------------------------------------------------------------
//  UNIT TESTS
//------------------------------------------------------------------------------
TEST( mint_mesh_curvilinear_mesh_DeathTest, invalid_construction )
{
  const IndexType N[]  = {   5,   5,   5 };
  double x[]           = { 0, 1, 2, 3, 4, 5 };

  // check 1st native constructor
  EXPECT_DEATH_IF_SUPPORTED( CurvilinearMesh(42,N), IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( CurvilinearMesh(3,nullptr), IGNORE_OUTPUT );

  // check 2nd native constructor
  EXPECT_DEATH_IF_SUPPORTED( CurvilinearMesh( -1,N[1], N[2]), IGNORE_OUTPUT );

  // check external constructor
  EXPECT_DEATH_IF_SUPPORTED( CurvilinearMesh(nullptr, x), IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( CurvilinearMesh(N,nullptr), IGNORE_OUTPUT );

#ifdef AXOM_MINT_USE_SIDRE

  sidre::DataStore ds;
  sidre::Group* root          = ds.getRoot();
  sidre::Group* valid_group   = root->createGroup( "mesh" );
  sidre::Group* particle_mesh = root->createGroup( "particle_mesh" );
  ParticleMesh( 3, 10, particle_mesh );

  // check pull constructor
  EXPECT_DEATH_IF_SUPPORTED( CurvilinearMesh(nullptr,""), IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( CurvilinearMesh(root,""), IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( CurvilinearMesh(particle_mesh,""), IGNORE_OUTPUT );

  // check 1st push constructor
  EXPECT_DEATH_IF_SUPPORTED( CurvilinearMesh( 42, N, valid_group ),
                             IGNORE_OUTPUT );
  EXPECT_EQ( valid_group->getNumGroups(), 0 );
  EXPECT_EQ( valid_group->getNumViews(), 0 );

  EXPECT_DEATH_IF_SUPPORTED( CurvilinearMesh( 3, nullptr, valid_group ),
                             IGNORE_OUTPUT );
  EXPECT_EQ( valid_group->getNumGroups(), 0 );
  EXPECT_EQ( valid_group->getNumViews(), 0 );

  EXPECT_DEATH_IF_SUPPORTED( CurvilinearMesh( 3, N, nullptr ),
                             IGNORE_OUTPUT );
  EXPECT_EQ( valid_group->getNumGroups(), 0 );
  EXPECT_EQ( valid_group->getNumViews(), 0 );

  // check 2nd push constructor
  EXPECT_DEATH_IF_SUPPORTED( CurvilinearMesh( nullptr, N[0] ),
                             IGNORE_OUTPUT );

  EXPECT_DEATH_IF_SUPPORTED( CurvilinearMesh( valid_group, -1 ),
                             IGNORE_OUTPUT );

#endif

}

//------------------------------------------------------------------------------
TEST( mint_mesh_curvilinear_mesh, native_constructor )
{
  constexpr int NDIMS  = 3;
  const IndexType N[]  = { 5, 6, 7 };
  const int64 extent[] = { 0, 4, 10, 15, 7, 13 };

  for ( int idim=1 ; idim <= NDIMS ; ++idim )
  {
    CurvilinearMesh* m = new CurvilinearMesh( idim, N );
    internal::check_constructor( m, STRUCTURED_CURVILINEAR_MESH, idim, N );
    EXPECT_FALSE( m->isExternal() );
    EXPECT_FALSE( m->hasSidreGroup() );
    m->setExtent( idim, extent );
    internal::check_node_extent( m, extent );
    set_coordinates( m );
    internal::check_create_fields( m );
    delete m;

    switch( idim )
    {
    case 1:
      m = new CurvilinearMesh( N[ 0 ] );
      break;
    case 2:
      m= new CurvilinearMesh( N[ 0 ], N[ 1 ] );
      break;
    default:
      EXPECT_EQ( idim, 3 );
      m = new CurvilinearMesh( N[ 0 ], N[ 1 ], N[ 2 ] );
    } // END switch

    internal::check_constructor( m, STRUCTURED_CURVILINEAR_MESH, idim, N );
    EXPECT_FALSE( m->isExternal( ) );
    EXPECT_FALSE( m->hasSidreGroup() );
    m->setExtent( idim, extent );
    internal::check_node_extent( m, extent );
    set_coordinates( m );
    internal::check_create_fields( m );
    delete m;
  } // END for all dimensions

}

//------------------------------------------------------------------------------
TEST( mint_mesh_curvilinear_mesh, external_constructor )
{
  constexpr int NDIMS = 3;
  const IndexType N[] = { 5, 6, 7 };
  const int64 extent[] = { 0, 4, 10, 15, 7, 13 };
  const IndexType maxNumNodes = N[ 0 ] * N[ 1 ] * N[ 2 ];

  double* x = new double[ maxNumNodes ];
  double* y = new double[ maxNumNodes ];
  double* z = new double[ maxNumNodes ];
  set_coordinates( maxNumNodes, x, y, z );

  IndexType curNumNodes = 1;
  IndexType curNumCells = 1;
  for ( int idim=1 ; idim <= NDIMS ; ++idim )
  {
    curNumNodes *= N[ idim - 1 ];
    curNumCells *= N[ idim - 1 ] - 1;
    CurvilinearMesh* m = nullptr;

    switch( idim )
    {
    case 1:
    {
      m = new CurvilinearMesh( N, x );
      EXPECT_EQ( x, m->getCoordinateArray( X_COORDINATE ) );
      check_coordinates( curNumNodes, idim,
                         m->getCoordinateArray( X_COORDINATE ) );
    }     // END 1D
    break;
    case 2:
    {
      m = new CurvilinearMesh( N, x, y );
      EXPECT_EQ( x, m->getCoordinateArray( X_COORDINATE ) );
      EXPECT_EQ( y, m->getCoordinateArray( Y_COORDINATE ) );
      check_coordinates( curNumNodes, idim,
                         m->getCoordinateArray( X_COORDINATE ),
                         m->getCoordinateArray( Y_COORDINATE )   );
    }     // END 2D
    break;
    default:
    {
      m = new CurvilinearMesh( N, x, y, z );
      EXPECT_EQ( x, m->getCoordinateArray( X_COORDINATE ) );
      EXPECT_EQ( y, m->getCoordinateArray( Y_COORDINATE ) );
      EXPECT_EQ( z, m->getCoordinateArray( Z_COORDINATE ) );
      check_coordinates( curNumNodes, idim,
                         m->getCoordinateArray( X_COORDINATE ),
                         m->getCoordinateArray( Y_COORDINATE ),
                         m->getCoordinateArray( Z_COORDINATE ) );
    }     // END 3D
    } // END switch

    EXPECT_TRUE( m != nullptr );
    internal::check_constructor( m, STRUCTURED_CURVILINEAR_MESH, idim, N );
    m->setExtent( idim, extent );
    internal::check_node_extent( m, extent );

    EXPECT_FALSE( m->hasSidreGroup() );
    EXPECT_TRUE( m->isExternal() );
    EXPECT_EQ( m->getNumberOfNodes(), curNumNodes );
    EXPECT_EQ( m->getNumberOfCells(), curNumCells );

    delete m;
    m = nullptr;

    // ensure array buffers are persistent
    check_coordinates( curNumNodes, idim, x, y, z );
  } // END for all dimensions

  delete[] x;
  delete[] y;
  delete[] z;
}

//------------------------------------------------------------------------------
#ifdef AXOM_MINT_USE_SIDRE

TEST( mint_mesh_curvilinear_mesh, sidre_constructor )
{
  constexpr int NDIMS  = 3;
  const IndexType N[]  = { 5, 6, 7 };
  const int64 extent[] = { 0, 4, 10, 15, 7, 13 };

  IndexType numNodes = 1;
  for ( int idim=1 ; idim <= NDIMS ; ++idim )
  {
    numNodes *= N[ idim - 1 ];

    // STEP 0: create a data-store with two (empty) groups
    sidre::DataStore ds;
    sidre::Group* root  = ds.getRoot();
    sidre::Group* m1grp = root->createGroup( "m1" );
    sidre::Group* m2grp = root->createGroup( "m2" );

    // STEP 1: populate meshes in the 2 Sidre groups using the 2 constructors.
    // BEGIN SCOPE
    {
      CurvilinearMesh* m = new CurvilinearMesh( idim, N, m1grp );
      EXPECT_TRUE( m->hasSidreGroup() );
      EXPECT_FALSE( m->isExternal() );
      internal::check_constructor( m, STRUCTURED_CURVILINEAR_MESH, idim, N );
      m->setExtent( idim, extent );
      internal::check_node_extent( m, extent );
      set_coordinates( m );
      internal::check_create_fields( m );
      delete m;

      switch ( idim )
      {
      case 1:
        m = new CurvilinearMesh( m2grp, N[ I_DIRECTION ] );
        break;
      case 2:
        m = new CurvilinearMesh( m2grp, N[ I_DIRECTION ],
                                 N[ J_DIRECTION ] );
        break;
      default:
        EXPECT_EQ( idim, 3 );
        m = new CurvilinearMesh( m2grp, N[ I_DIRECTION ],
                                 N[ J_DIRECTION ],
                                 N[ K_DIRECTION ] );
      } // END switch

      EXPECT_TRUE( m->hasSidreGroup() );
      EXPECT_FALSE( m->isExternal() );
      internal::check_constructor( m, STRUCTURED_CURVILINEAR_MESH, idim, N );
      m->setExtent( idim, extent );
      internal::check_node_extent( m, extent );
      set_coordinates( m );
      internal::check_create_fields( m );
      delete m;
    }
    // END SCOPE

    // STEP 2: pull the meshes from sidre in to new instances.
    // BEGIN SCOPE
    {
      // check m1
      CurvilinearMesh* m = new CurvilinearMesh( m1grp );
      EXPECT_TRUE( m->hasSidreGroup() );
      EXPECT_FALSE( m->isExternal() );
      internal::check_constructor( m, STRUCTURED_CURVILINEAR_MESH, idim, N );
      internal::check_fields( m, true );
      EXPECT_EQ( idim, m->getDimension() );
      EXPECT_EQ( numNodes, m->getNumberOfNodes() );
      check_coordinates( m );
      delete m;

      // check m2
      m = new CurvilinearMesh( m2grp );
      EXPECT_TRUE( m->hasSidreGroup() );
      EXPECT_FALSE( m->isExternal() );
      internal::check_constructor( m, STRUCTURED_CURVILINEAR_MESH, idim, N );
      internal::check_fields( m, true );
      EXPECT_EQ( idim, m->getDimension() );
      EXPECT_EQ( numNodes, m->getNumberOfNodes() );
      check_coordinates( m );
      delete m;
    }
    // END SCOPE

    EXPECT_TRUE( blueprint::isValidRootGroup( m1grp ) );
    EXPECT_TRUE( blueprint::isValidRootGroup( m2grp ) );

  } // END for all dimensions

}

#endif /* AXOM_MINT_USE_SIDRE */

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

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
#include "axom/mint/config.hpp"          // for compile-time definitions

// Mint includes
#include "axom/mint/mesh/blueprint.hpp"       // for blueprint functions
#include "axom/mint/mesh/CellTypes.hpp"       // for CellType enum
#include "axom/mint/mesh/RectilinearMesh.hpp" // for RectilinearMesh
#include "axom/mint/mesh/ParticleMesh.hpp"    // for ParticleMesh

// Slic includes
#include "axom/slic/interface/slic.hpp"            // for slic macros

// Sidre includes
#ifdef MINT_USE_SIDRE
#include "axom/sidre/core/sidre.hpp"
namespace sidre = axom::sidre;
#endif

#include "gtest/gtest.h"            // for gtest macros

// C/C++ includes
#include <cmath>  // for exp
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
void check_create_field( RectilinearMesh* m,
                         int association,
                         const std::string& name,
                         int numComponents=1 )
{
  EXPECT_TRUE( m != AXOM_NULLPTR );
  EXPECT_FALSE( m->hasField( name, association ) );

  double* f = m->createField< double >( name, association, numComponents );
  EXPECT_TRUE( m->hasField( name, association ) );

  IndexType expected_num_tuples;
  if ( association == NODE_CENTERED )
  {
    expected_num_tuples = m->getNumberOfNodes();
  }
  else if ( association == CELL_CENTERED )
  {
    expected_num_tuples = m->getNumberOfCells();
  }
  else if ( association == FACE_CENTERED )
  {
    expected_num_tuples = m->getNumberOfFaces();
  }
  else
  {
    expected_num_tuples = m->getNumberOfEdges();
  }

  const Field* field = m->getFieldData( association )->getField( name );
  EXPECT_TRUE( field != nullptr );
  EXPECT_EQ( f, Field::getDataPtr< double >( field ) );
  EXPECT_EQ( numComponents, field->getNumComponents() );
  EXPECT_EQ( expected_num_tuples, field->getNumTuples() );

  for ( IndexType i = 0; i < expected_num_tuples * numComponents; ++i )
  {
    f[ i ] = 3.14 * i;
  }
}

//------------------------------------------------------------------------------
void check_field_values( const Field* field )
{
  const IndexType num_values = field->getNumTuples() * field->getNumComponents();
  const double* values = Field::getDataPtr< double >( field );
  for ( IndexType i = 0; i < num_values; ++i )
  {
    EXPECT_EQ( values[ i ], 3.14 * i );
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
void check_constructor( const RectilinearMesh* m,
                        int expected_dimension,
                        const IndexType* expected_node_dimensions )
{
  EXPECT_TRUE( m != nullptr );
  EXPECT_TRUE( (expected_dimension >= 1) && (expected_dimension <= 3) );

  const int mesh_dimension = m->getDimension();

  EXPECT_EQ( mesh_dimension, expected_dimension );
  EXPECT_EQ( m->getMeshType(), STRUCTURED_RECTILINEAR_MESH );
  EXPECT_FALSE( m->hasExplicitConnectivity( ) );
  EXPECT_FALSE( m->hasMixedCellTypes( ) );
  EXPECT_TRUE( m->hasExplicitCoordinates( ) );

  CellType expected_cell_type =
    ( mesh_dimension==3 ) ? HEX : ( ( mesh_dimension==2 ) ? QUAD : SEGMENT );
  EXPECT_EQ( m->getCellType(), expected_cell_type );

  for ( int i=0 ; i < mesh_dimension ; ++i )
  {
    EXPECT_EQ( m->getNodeExtent( i ),expected_node_dimensions[ i ] );
    EXPECT_TRUE( m->getCoordinateArray( i ) != AXOM_NULLPTR );
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
    RectilinearMesh m1( idim, N );
    check_constructor( &m1, idim, N );
    check_fill_coords( &m1 );
    check_create_field( &m1, NODE_CENTERED, "n1" );
    check_create_field( &m1, CELL_CENTERED, "c1", 2 );
    check_create_field( &m1, FACE_CENTERED, "f1", 3 );
    check_create_field( &m1, EDGE_CENTERED, "e1", 4 );

    EXPECT_FALSE( m1.isExternal( ) );
    EXPECT_FALSE( m1.hasSidreGroup( ) );

    RectilinearMesh* m2 = nullptr;
    switch ( idim )
    {
    case 1:
      m2 = new RectilinearMesh( N[ 0 ] );
      break;
    case 2:
      m2 = new RectilinearMesh( N[ 0 ], N[ 1 ] );
      break;
    default:
      EXPECT_EQ( idim, 3 );
      m2 = new RectilinearMesh( N[ 0 ], N[ 1 ], N[ 2 ] );
    }

    check_constructor( m2, idim, N );
    check_fill_coords( m2 );
    check_create_field( m2, NODE_CENTERED, "n1" );
    check_create_field( m2, CELL_CENTERED, "c1", 3 );
    check_create_field( m2, FACE_CENTERED, "f1", 3 );
    check_create_field( m2, EDGE_CENTERED, "e1", 4 );

    delete m2;
    m2 = nullptr;
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
    check_constructor( m, idim, N );
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
      check_constructor( m1, idim, N );
      check_fill_coords( m1 );
      check_create_field( m1, NODE_CENTERED, "n1" );
      check_create_field( m1, CELL_CENTERED, "c1", 2 );
      check_create_field( m1, FACE_CENTERED, "f1", 3 );
      check_create_field( m1, EDGE_CENTERED, "e1", 4 );

      EXPECT_TRUE( m1->hasSidreGroup() );
      EXPECT_FALSE( m1->isExternal() );

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

      check_constructor( m2, idim, N );
      check_fill_coords( m2 );
      check_create_field( m2, NODE_CENTERED, "n1" );
      check_create_field( m2, CELL_CENTERED, "c1", 2 );
      check_create_field( m2, FACE_CENTERED, "f1", 3 );
      check_create_field( m2, EDGE_CENTERED, "e1", 4 );

      EXPECT_TRUE( m2->hasSidreGroup() );
      EXPECT_FALSE( m2->isExternal() );

      delete m2;
      m2 = nullptr;
    }
    // END SCOPE

    // STEP 2: pull the 2 meshes from sidre and check correctness
    // BEGIN SCOPE
    {
      const FieldData* fd = nullptr;
      const Field* field  = nullptr;

      // check m1
      RectilinearMesh* m1 = new RectilinearMesh( m1grp );
      check_constructor( m1, idim, N );
      EXPECT_TRUE( m1->hasSidreGroup() );
      EXPECT_FALSE( m1->isExternal() );

      for ( int ii=0 ; ii < idim ; ++ii )
      {
        check_coordinate( m1->getCoordinateArray( ii ), expected_coords,
                          m1->getNodeExtent( ii ) );
      }

      // check node-centered field on m1
      fd    = m1->getFieldData( NODE_CENTERED );
      field = fd->getField( "n1" );
      EXPECT_EQ( field->getNumTuples(), m1->getNumberOfNodes() );
      EXPECT_EQ( field->getNumComponents(), 1 );
      EXPECT_TRUE( field->isInSidre() );
      EXPECT_FALSE( field->isExternal() );
      check_field_values( field );

      // check cell-centered field on m1
      fd    = m1->getFieldData( CELL_CENTERED );
      field = fd->getField( "c1" );
      EXPECT_EQ( field->getNumTuples(), m1->getNumberOfCells() );
      EXPECT_EQ( field->getNumComponents(), 2 );
      EXPECT_TRUE( field->isInSidre() );
      EXPECT_FALSE( field->isExternal() );
      check_field_values( field );

      // check face-centered field on m1
      fd    = m1->getFieldData( FACE_CENTERED );
      field = fd->getField( "f1" );
      EXPECT_EQ( field->getNumTuples(), m1->getNumberOfFaces() );
      EXPECT_EQ( field->getNumComponents(), 3 );
      EXPECT_TRUE( field->isInSidre() );
      EXPECT_FALSE( field->isExternal() );
      check_field_values( field );

      // check edge-centered field on m1
      fd    = m1->getFieldData( EDGE_CENTERED );
      field = fd->getField( "e1" );
      EXPECT_EQ( field->getNumTuples(), m1->getNumberOfEdges() );
      EXPECT_EQ( field->getNumComponents(), 4 );
      EXPECT_TRUE( field->isInSidre() );
      EXPECT_FALSE( field->isExternal() );
      check_field_values( field );

      delete m1;
      m1 = nullptr;

      // check m2
      RectilinearMesh* m2 = new RectilinearMesh( m2grp );
      check_constructor( m2, idim, N );
      EXPECT_TRUE( m2->hasSidreGroup() );
      EXPECT_FALSE( m2->isExternal() );

      for ( int ii=0 ; ii < idim ; ++ii )
      {
        check_coordinate( m2->getCoordinateArray( ii ), expected_coords,
                          m2->getNodeExtent( ii ) );
      }

      // check node-centered field on m2
      fd    = m2->getFieldData( NODE_CENTERED );
      field = fd->getField( "n1" );
      EXPECT_EQ( field->getNumTuples(), m2->getNumberOfNodes() );
      EXPECT_EQ( field->getNumComponents(), 1 );
      EXPECT_TRUE( field->isInSidre() );
      EXPECT_FALSE( field->isExternal() );
      check_field_values( field );

      // check cell-centered field on m2
      fd    = m2->getFieldData( CELL_CENTERED );
      field = fd->getField( "c1" );
      EXPECT_EQ( field->getNumTuples(), m2->getNumberOfCells() );
      EXPECT_EQ( field->getNumComponents(), 2 );
      EXPECT_TRUE( field->isInSidre() );
      EXPECT_FALSE( field->isExternal() );
      check_field_values( field );

      // check face-centered field on m2
      fd    = m2->getFieldData( FACE_CENTERED );
      field = fd->getField( "f1" );
      EXPECT_EQ( field->getNumTuples(), m2->getNumberOfFaces() );
      EXPECT_EQ( field->getNumComponents(), 3 );
      EXPECT_TRUE( field->isInSidre() );
      EXPECT_FALSE( field->isExternal() );
      check_field_values( field );

      // check edge-centered field on m2
      fd    = m2->getFieldData( EDGE_CENTERED );
      field = fd->getField( "e1" );
      EXPECT_EQ( field->getNumTuples(), m2->getNumberOfEdges() );
      EXPECT_EQ( field->getNumComponents(), 4 );
      EXPECT_TRUE( field->isInSidre() );
      EXPECT_FALSE( field->isExternal() );
      check_field_values( field );

      delete m2;
      m2 = nullptr;
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
  const IndexType N[] = {    5,    5,    5 };
  double x[ 5 ];
  double y[ 5 ];
  double z[ 5 ];
  exponential_distribution( 42.0, 5, x );
  exponential_distribution( 0.0,  5, y );
  exponential_distribution( 0.0,  5, z );

  RectilinearMesh m( N, x, y, z );
  check_constructor( &m, 3, N );

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

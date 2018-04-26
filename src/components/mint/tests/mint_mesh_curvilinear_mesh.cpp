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
#include "mint/config.hpp"          // for compile-time type definitions

// Mint includes
#include "mint/blueprint.hpp"       // for blueprint functions
#include "mint/CellTypes.hpp"       // for CellTypes enum definition
#include "mint/CurvilinearMesh.hpp" // for mint::CurivilinearMesh
#include "mint/ParticleMesh.hpp"    // for ParticleMesh
#include "mint/UniformMesh.hpp"     // for UniformMesh

// Slic includes
#include "slic/slic.hpp"            // for slic macros

// Sidre includes
#ifdef MINT_USE_SIDRE
#include "sidre/sidre.hpp"
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
inline int dim( const double* AXOM_NOT_USED(x),
                const double* y,
                const double* z )
{
  return ( ( z != AXOM_NULLPTR ) ? 3 : ( (y != AXOM_NULLPTR ) ? 2 : 1 ) );
}

//------------------------------------------------------------------------------
template < typename T >
void delete_array( T* ptr )
{
  if ( ptr != AXOM_NULLPTR )
  {
    delete [ ] ptr;
    ptr = AXOM_NULLPTR;
  }
}

//------------------------------------------------------------------------------
void check_coordinates( const UniformMesh* m,
                        double* x,
                        double* y=AXOM_NULLPTR,
                        double* z=AXOM_NULLPTR  )
{
  EXPECT_TRUE( m != AXOM_NULLPTR );
  EXPECT_TRUE( x != AXOM_NULLPTR );

  const int mesh_dimension = m->getDimension();
  EXPECT_EQ( mesh_dimension, dim(x,y,z) );

  switch( mesh_dimension )
  {
  case 1:
    {
      const IndexType Ni = m->getNumberOfNodesAlongDim( I_DIRECTION );
      EXPECT_EQ( Ni, m->getNumberOfNodes() );
      EXPECT_EQ( Ni, m->getNumberOfNodes() );
      for ( IndexType i=0; i < Ni; ++i )
      {
        EXPECT_DOUBLE_EQ( x[ i ],  m->evaluateCoordinate( i, I_DIRECTION ) );
      } // END for all i
    } // END 1D
    break;
  case 2:
    {
      EXPECT_TRUE( y != AXOM_NULLPTR );

      const IndexType Ni = m->getNumberOfNodesAlongDim( I_DIRECTION );
      const IndexType Nj = m->getNumberOfNodesAlongDim( J_DIRECTION );
      for ( IndexType j=0; j < Nj; ++j )
      {
        for ( IndexType i=0; i < Ni; ++i )
        {

          const IndexType nodeIdx = m->getLinearIndex( i,j );
          const double xx = m->evaluateCoordinate( i, I_DIRECTION );
          const double yy = m->evaluateCoordinate( j, J_DIRECTION );
          EXPECT_DOUBLE_EQ( x[ nodeIdx ], xx );
          EXPECT_DOUBLE_EQ( y[ nodeIdx ], yy );

        } // END for all i
      } // END for all j

    } // END 2D
    break;
  default:
    {
      EXPECT_TRUE( mesh_dimension==3 );
      EXPECT_TRUE( y != AXOM_NULLPTR );
      EXPECT_TRUE( z != AXOM_NULLPTR );

      const IndexType Ni = m->getNumberOfNodesAlongDim( I_DIRECTION );
      const IndexType Nj = m->getNumberOfNodesAlongDim( J_DIRECTION );
      const IndexType Nk = m->getNumberOfNodesAlongDim( K_DIRECTION );

      for ( IndexType k=0; k < Nk; ++k )
      {
        for ( IndexType j=0; j < Nj; ++j )
        {
          for ( IndexType i=0; i < Ni; ++i )
          {
            const IndexType nodeIdx = m->getLinearIndex( i, j, k );
            const double xx = m->evaluateCoordinate( i, I_DIRECTION );
            const double yy = m->evaluateCoordinate( j, J_DIRECTION );
            const double zz = m->evaluateCoordinate( k, K_DIRECTION );
            EXPECT_DOUBLE_EQ( x[ nodeIdx ], xx );
            EXPECT_DOUBLE_EQ( y[ nodeIdx ], yy );
            EXPECT_DOUBLE_EQ( z[ nodeIdx ], zz );

          } // END for all i
        } // END for all j
      } // END for all k
    } // END 3D

  } // END switch
}

//------------------------------------------------------------------------------
void get_coordinates( const UniformMesh* m,
                      double* x,
                      double* y=AXOM_NULLPTR ,
                      double* z=AXOM_NULLPTR    )
{
  EXPECT_TRUE( m != AXOM_NULLPTR );
  EXPECT_TRUE( x != AXOM_NULLPTR );

  const int mesh_dimension = m->getDimension();

  switch( mesh_dimension )
  {
  case 1:
    {
      const IndexType Ni = m->getNumberOfNodesAlongDim( I_DIRECTION );
      EXPECT_EQ( Ni, m->getNumberOfNodes() );
      for ( IndexType i=0; i < Ni; ++i )
      {
        x[ i ] = m->evaluateCoordinate( i, I_DIRECTION );
      } // END for all i

    } // END 1D
    break;
  case 2:
    {
      EXPECT_TRUE( y != AXOM_NULLPTR );
      const IndexType Ni = m->getNumberOfNodesAlongDim( I_DIRECTION );
      const IndexType Nj = m->getNumberOfNodesAlongDim( J_DIRECTION );
      for ( IndexType j=0; j < Nj; ++j )
      {
        for ( IndexType i=0; i < Ni; ++i )
        {
          const IndexType nodeIdx = m->getLinearIndex( i,j );
          x[ nodeIdx ] = m->evaluateCoordinate( i, I_DIRECTION );
          y[ nodeIdx ] = m->evaluateCoordinate( j, J_DIRECTION );
        } // END for all i
      } // END for all j

    } // END 2D
    break;
  default:
    {
      EXPECT_TRUE( y != AXOM_NULLPTR );
      EXPECT_TRUE( z != AXOM_NULLPTR );
      const IndexType Ni = m->getNumberOfNodesAlongDim( I_DIRECTION );
      const IndexType Nj = m->getNumberOfNodesAlongDim( J_DIRECTION );
      const IndexType Nk = m->getNumberOfNodesAlongDim( K_DIRECTION );

      for ( IndexType k=0; k < Nk; ++k )
      {
        for ( IndexType j=0; j < Nj; ++j )
        {
          for ( IndexType i=0; i < Ni; ++i )
          {
            const IndexType nodeIdx = m->getLinearIndex( i, j, k );
            x[ nodeIdx ] = m->evaluateCoordinate( i, I_DIRECTION );
            y[ nodeIdx ] = m->evaluateCoordinate( j, J_DIRECTION );
            z[ nodeIdx ] = m->evaluateCoordinate( k, K_DIRECTION );

          } // END for all i
        } // END for all j
      } // END for all k

    } // END 3D

  } // END switch
}

//------------------------------------------------------------------------------
void check_create_field( CurvilinearMesh* m,
                         int association,
                         const std::string& name,
                         int numComponents=1 )
{
  EXPECT_TRUE( m != AXOM_NULLPTR );
  EXPECT_TRUE( (association==NODE_CENTERED) ||
               (association==CELL_CENTERED) );

  EXPECT_FALSE( m->hasField( name, association ) );

  double* f = m->createField< double >( name, association, numComponents );
  EXPECT_TRUE( f != AXOM_NULLPTR );
  EXPECT_TRUE( m->hasField( name, association ) );

  IndexType expected_num_tuples = ( association==NODE_CENTERED ) ?
        m->getNumberOfNodes() : m->getNumberOfCells() ;

  const Field* field = m->getFieldData( association )->getField( name );
  EXPECT_TRUE( field != AXOM_NULLPTR );
  EXPECT_EQ( f, Field::getDataPtr< double >( field ) );
  EXPECT_EQ( numComponents, field->getNumComponents() );
  EXPECT_EQ( expected_num_tuples, field->getNumTuples() );
}

//------------------------------------------------------------------------------
void check_fill_coords( CurvilinearMesh* m, double MAGIC_VAL=42.0 )
{
  EXPECT_TRUE( m != AXOM_NULLPTR );

  const IndexType numNodes = m->getNumberOfNodes();

  const int mesh_dimension = m->getDimension();
  switch( mesh_dimension )
  {
  case 1:
    {
      double* x = m->getCoordinateArray( X_COORDINATE );
      EXPECT_TRUE( x != AXOM_NULLPTR );

      for ( IndexType i=0; i < numNodes; ++i )
      {
        x[ i ] = MAGIC_VAL;
      } // END for

    } // END if 1D
    break;
  case 2:
    {
      double* x = m->getCoordinateArray( X_COORDINATE );
      double* y = m->getCoordinateArray( Y_COORDINATE );
      EXPECT_TRUE( x != AXOM_NULLPTR );
      EXPECT_TRUE( y != AXOM_NULLPTR );

      for ( IndexType i=0; i < numNodes; ++i )
      {
        x[ i ] = y[ i ] = MAGIC_VAL;
      } // END for

    } // END if 2D
    break;
  default:
    {
      EXPECT_TRUE( mesh_dimension==3 );

      double* x = m->getCoordinateArray( X_COORDINATE );
      double* y = m->getCoordinateArray( Y_COORDINATE );
      double* z = m->getCoordinateArray( Z_COORDINATE );
      EXPECT_TRUE( x != AXOM_NULLPTR );
      EXPECT_TRUE( y != AXOM_NULLPTR );
      EXPECT_TRUE( z != AXOM_NULLPTR );

      for ( IndexType i=0; i < numNodes; ++i )
      {
        x[ i ] = y[ i ] = z[ i ] = MAGIC_VAL;
      } // END for

    } // END 3D

  } // END switch

}

//------------------------------------------------------------------------------
void check_constructor( const CurvilinearMesh* m,
                        int expected_dimension,
                        const IndexType* expected_node_dimensions,
                        const int64* expected_extent
                        )
{
  EXPECT_TRUE( m != AXOM_NULLPTR );
  EXPECT_TRUE( (expected_dimension >= 1) && (expected_dimension <= 3) );

  const int mesh_dimension = m->getDimension();

  EXPECT_EQ( mesh_dimension, expected_dimension );
  EXPECT_EQ( m->getMeshType(), STRUCTURED_CURVILINEAR_MESH );
  EXPECT_FALSE( m->hasExplicitConnectivity() );
  EXPECT_FALSE( m->hasMixedCellTypes() );
  EXPECT_TRUE( m->hasExplicitCoordinates() );

  CellType expected_cell_type =
      ( mesh_dimension==3 ) ? HEX : ( ( mesh_dimension==2 ) ? QUAD : SEGMENT );
  EXPECT_EQ( m->getCellType(), expected_cell_type );

  const Extent* ext = m->getExtent();
  EXPECT_TRUE( ext != AXOM_NULLPTR );

  for ( int i=0; i < mesh_dimension; ++i )
  {
    const int offset = i*2;
    EXPECT_EQ( ext->min( i ), expected_extent[ offset ] );
    EXPECT_EQ( ext->max( i ), expected_extent[ offset+1 ] );
    EXPECT_EQ( m->getNumberOfNodesAlongDim( i ),expected_node_dimensions[ i ] );

    EXPECT_TRUE( m->getCoordinateArray( i ) != AXOM_NULLPTR );
  }

  EXPECT_EQ( m->getNumberOfNodes(), ext->getNumNodes() );
  EXPECT_EQ( m->getNumberOfCells(), ext->getNumCells() );

}

} // END namespace

//------------------------------------------------------------------------------
//  UNIT TESTS
//------------------------------------------------------------------------------
TEST( mint_mesh_curvilinear_mesh_DeathTest, invalid_construction )
{
  const int64 ext[ ]   = { 0,4, 0,4, 0,4 };
  const IndexType N[ ] = {   5,   5,   5 };
  double x[ ]          = { 0, 1, 2, 3, 4, 5 };

  // check 1st native constructor
  EXPECT_DEATH_IF_SUPPORTED( CurvilinearMesh(42,ext), IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( CurvilinearMesh(3,AXOM_NULLPTR), IGNORE_OUTPUT );

  // check 2nd native constructor
  EXPECT_DEATH_IF_SUPPORTED( CurvilinearMesh( -1,N[1], N[2]), IGNORE_OUTPUT );

  // check external constructor
  EXPECT_DEATH_IF_SUPPORTED( CurvilinearMesh(AXOM_NULLPTR, x), IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( CurvilinearMesh(ext,AXOM_NULLPTR), IGNORE_OUTPUT );

#ifdef MINT_USE_SIDRE

  sidre::DataStore ds;
  sidre::Group* root          = ds.getRoot();
  sidre::Group* valid_group   = root->createGroup( "mesh" );
  sidre::Group* particle_mesh = root->createGroup( "particle_mesh" );
  ParticleMesh( 3, 10, particle_mesh );

  // check pull constructor
  EXPECT_DEATH_IF_SUPPORTED( CurvilinearMesh(AXOM_NULLPTR,""), IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( CurvilinearMesh(root,""), IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( CurvilinearMesh(particle_mesh,""), IGNORE_OUTPUT );

  // check 1st push constructor
  EXPECT_DEATH_IF_SUPPORTED( CurvilinearMesh( 42, ext, valid_group ),
                             IGNORE_OUTPUT );
  EXPECT_EQ( valid_group->getNumGroups(), 0 );
  EXPECT_EQ( valid_group->getNumViews(), 0 );

  EXPECT_DEATH_IF_SUPPORTED( CurvilinearMesh( 3, AXOM_NULLPTR, valid_group ),
                             IGNORE_OUTPUT );
  EXPECT_EQ( valid_group->getNumGroups(), 0 );
  EXPECT_EQ( valid_group->getNumViews(), 0 );

  EXPECT_DEATH_IF_SUPPORTED( CurvilinearMesh( 3, ext, AXOM_NULLPTR ),
                             IGNORE_OUTPUT );
  EXPECT_EQ( valid_group->getNumGroups(), 0 );
  EXPECT_EQ( valid_group->getNumViews(), 0 );

  // check 2nd push constructor
  EXPECT_DEATH_IF_SUPPORTED( CurvilinearMesh( AXOM_NULLPTR, N[0] ),
                             IGNORE_OUTPUT );

  EXPECT_DEATH_IF_SUPPORTED( CurvilinearMesh( valid_group, -1 ),
                             IGNORE_OUTPUT );

#endif

}

//------------------------------------------------------------------------------
TEST( mint_mesh_curvilinear_mesh, native_constructor )
{
  constexpr int NDIMS  = 3;
  constexpr double MAGIC_NUM = 42.0;
  const int64 ext[ ]   = { 0,4, 0,4, 0,4 };
  const IndexType N[ ] = {   5,   5,   5 };

  for ( int idim=1; idim <= NDIMS; ++idim )
  {

    CurvilinearMesh m1( idim, ext );
    check_constructor( &m1, idim, N, ext );
    EXPECT_FALSE( m1.isExternal() );
    EXPECT_FALSE( m1.hasSidreGroup() );

    check_fill_coords( &m1, MAGIC_NUM );
    check_create_field( &m1, NODE_CENTERED, "n1" );
    check_create_field( &m1, CELL_CENTERED, "c1", 3 );

    CurvilinearMesh *m2 = AXOM_NULLPTR;
    switch( idim )
    {
    case 1:
      m2 = new CurvilinearMesh( N[ 0 ] );
      break;
    case 2:
      m2= new CurvilinearMesh( N[ 0 ], N[ 1 ] );
      break;
    default:
      EXPECT_EQ( idim, 3 );
      m2 = new CurvilinearMesh( N[ 0 ], N[ 1 ], N[ 2 ] );
    } // END switch

    check_constructor( m2, idim, N, ext );
    EXPECT_FALSE( m2->isExternal( ) );
    EXPECT_FALSE( m2->hasSidreGroup() );

    check_fill_coords( m2, MAGIC_NUM );
    check_create_field( m2, NODE_CENTERED, "n1" );
    check_create_field( m2, CELL_CENTERED, "c1", 3 );

    delete m2;
    m2 = AXOM_NULLPTR;
  } // END for all dimensions

}

//------------------------------------------------------------------------------
TEST( mint_mesh_curvilinear_mesh, external_constructor )
{
  constexpr int NDIMS = 3;
  const int64  ext[]  = { 0,10, 0,10, 0,10  };
  const IndexType N[] = {   11,   11,   11  };
  const double lo[]   = { 0.0, 0.0, 0.0 };
  const double hi[]   = { 5.0, 5.0, 5.0 };

  for ( int idim=1; idim <= NDIMS; ++idim )
  {
    UniformMesh um( idim, ext, lo, hi );
    const IndexType numNodes = um.getNumberOfNodes();

    double* x = AXOM_NULLPTR;
    double* y = AXOM_NULLPTR;
    double* z = AXOM_NULLPTR;

    CurvilinearMesh* m = AXOM_NULLPTR;

    switch( idim )
    {
    case 1:
      {
        x = new double[ numNodes ];
        get_coordinates( &um, x );

        m = new CurvilinearMesh( ext, x );
        EXPECT_EQ( x, m->getCoordinateArray( X_COORDINATE ) );
        check_coordinates( &um, m->getCoordinateArray( X_COORDINATE ) );
      } // END 1D
      break;
    case 2:
      {
        x = new double[ numNodes ];
        y = new double[ numNodes ];
        get_coordinates( &um, x, y );

        m = new CurvilinearMesh( ext, x, y );
        EXPECT_EQ( x, m->getCoordinateArray( X_COORDINATE ) );
        EXPECT_EQ( y, m->getCoordinateArray( Y_COORDINATE ) );
        check_coordinates( &um,
                           m->getCoordinateArray( X_COORDINATE ),
                           m->getCoordinateArray( Y_COORDINATE )   );
      } // END 2D
      break;
    default:
      {
        x = new double[ numNodes ];
        y = new double[ numNodes ];
        z = new double[ numNodes ];
        get_coordinates( &um, x, y, z );

        m = new CurvilinearMesh( ext, x, y, z );
        EXPECT_EQ( x, m->getCoordinateArray( X_COORDINATE ) );
        EXPECT_EQ( y, m->getCoordinateArray( Y_COORDINATE ) );
        EXPECT_EQ( z, m->getCoordinateArray( Z_COORDINATE ) );
        check_coordinates( &um,
                           m->getCoordinateArray( X_COORDINATE ),
                           m->getCoordinateArray( Y_COORDINATE ),
                           m->getCoordinateArray( Z_COORDINATE ) );

      } // END 3D

    } // END switch

    EXPECT_TRUE( m != AXOM_NULLPTR );
    check_constructor( m, idim, N, ext );

    EXPECT_FALSE( m->hasSidreGroup() );
    EXPECT_TRUE( m->isExternal() );
    EXPECT_EQ( m->getNumberOfNodes(), um.getNumberOfNodes() );
    EXPECT_EQ( m->getNumberOfCells(), um.getNumberOfCells() );

    delete m;
    m = AXOM_NULLPTR;

    // ensure array buffers are persistent
    check_coordinates( &um, x, y, z );

    delete_array( x );
    delete_array( y );
    delete_array( z );
  } // END for all dimensions

}

//------------------------------------------------------------------------------
#ifdef MINT_USE_SIDRE

TEST( mint_mesh_curvilinear_mesh, sidre_constructor )
{
  constexpr int NDIMS  = 3;
  constexpr double MAGIC_NUM = 42.0;
  const int64 ext[ ]   = { 0,4, 0,4, 0,4 };
  const IndexType N[ ] = {   5,   5,   5 };

  for ( int idim=1; idim <= NDIMS; ++idim )
  {
    // STEP 0: create a data-store with two (empty) groups
    sidre::DataStore ds;
    sidre::Group* root  = ds.getRoot();
    sidre::Group* m1grp = root->createGroup( "m1" );
    sidre::Group* m2grp = root->createGroup( "m2" );

    // STEP 1: populate meshes in the 2 Sidre groups using the 2 construstors.
    // BEGIN SCOPE
    {
      CurvilinearMesh* m1 = new CurvilinearMesh( idim, ext, m1grp );
      EXPECT_TRUE( m1->hasSidreGroup() );
      EXPECT_FALSE( m1->isExternal() );
      check_constructor( m1, idim, N, ext );
      check_fill_coords( m1, MAGIC_NUM );
      check_create_field( m1, NODE_CENTERED, "n1" );
      check_create_field( m1, CELL_CENTERED, "c1", 3 );

      EXPECT_FALSE( m1->isExternal() );
      EXPECT_TRUE( m1->hasSidreGroup() );

      CurvilinearMesh* m2 = AXOM_NULLPTR;
      switch ( idim )
      {
      case 1:
        m2 = new CurvilinearMesh( m2grp, N[ I_DIRECTION ] );
        break;
      case 2:
        m2 = new CurvilinearMesh( m2grp, N[ I_DIRECTION ],
                                         N[ J_DIRECTION ] );
        break;
      default:
        EXPECT_EQ( idim, 3 );
        m2 = new CurvilinearMesh( m2grp, N[ I_DIRECTION ],
                                         N[ J_DIRECTION ],
                                         N[ K_DIRECTION ] );

      } // END switch

      EXPECT_TRUE( m2->hasSidreGroup() );
      EXPECT_FALSE( m2->isExternal() );
      check_constructor( m2, idim, N, ext );
      check_fill_coords( m2, MAGIC_NUM );
      check_create_field( m2, NODE_CENTERED, "n1", 3 );
      check_create_field( m2, CELL_CENTERED, "c1" );

      delete m1;
      m1 = AXOM_NULLPTR;

      delete m2;
      m2 = AXOM_NULLPTR;
    }
    // END SCOPE

    // STEP 2: pull the meshes from sidre in to new instances.
    // BEGIN SCOPE
    {
      const FieldData* fd = AXOM_NULLPTR;
      const Field* field  = AXOM_NULLPTR;

      // check m1
      CurvilinearMesh* m1 = new CurvilinearMesh( m1grp );
      EXPECT_TRUE( m1->hasSidreGroup() );
      EXPECT_FALSE( m1->isExternal() );
      check_constructor( m1, idim, N, ext );
      EXPECT_TRUE( m1->hasField( "n1", NODE_CENTERED ) );
      EXPECT_TRUE( m1->hasField( "c1", CELL_CENTERED ) );

      // check node-centered field on m1
      fd    = m1->getFieldData( NODE_CENTERED );
      field = fd->getField( "n1" );
      EXPECT_EQ( field->getNumTuples(), m1->getNumberOfNodes() );
      EXPECT_EQ( field->getNumComponents(), 1 );
      EXPECT_TRUE( field->isInSidre() );
      EXPECT_FALSE( field->isExternal() );

      // check cell-centered field on m1
      fd    = m1->getFieldData( CELL_CENTERED );
      field = fd->getField( "c1" );
      EXPECT_EQ( field->getNumTuples(), m1->getNumberOfCells() );
      EXPECT_EQ( field->getNumComponents(), 3 );
      EXPECT_TRUE( field->isInSidre() );
      EXPECT_FALSE( field->isExternal() );

      // check m2
      CurvilinearMesh* m2 = new CurvilinearMesh( m2grp );
      EXPECT_TRUE( m2->hasSidreGroup() );
      EXPECT_FALSE( m2->isExternal() );
      check_constructor( m2, idim, N, ext );
      EXPECT_TRUE( m2->hasField( "n1", NODE_CENTERED ) );
      EXPECT_TRUE( m2->hasField( "c1", CELL_CENTERED ) );

      // check node-centered field on m2
      fd    = m2->getFieldData( NODE_CENTERED );
      field = fd->getField( "n1" );
      EXPECT_EQ( field->getNumTuples(), m2->getNumberOfNodes() );
      EXPECT_EQ( field->getNumComponents(), 3 );
      EXPECT_TRUE( field->isInSidre() );
      EXPECT_FALSE( field->isExternal() );

      // check cell-centered field on m2
      fd    = m2->getFieldData( CELL_CENTERED );
      field = fd->getField( "c1" );
      EXPECT_EQ( field->getNumTuples(), m2->getNumberOfCells() );
      EXPECT_EQ( field->getNumComponents(), 1 );
      EXPECT_TRUE( field->isInSidre() );
      EXPECT_FALSE( field->isExternal() );

      // check nodes
      EXPECT_EQ( m1->getNumberOfNodes(), m2->getNumberOfNodes() );
      EXPECT_EQ( m1->getNumberOfCells(), m2->getNumberOfCells() );

      const IndexType numNodes = m1->getNumberOfNodes();
      const int ndims          = m1->getDimension();
      EXPECT_EQ( idim, ndims );

      for ( int i=0; i < ndims; ++i )
      {
        const double *x1 = m1->getCoordinateArray( i );
        EXPECT_TRUE( x1 != AXOM_NULLPTR );

        const double* x2 = m2->getCoordinateArray( i );
        EXPECT_TRUE( x2 != AXOM_NULLPTR );

        for ( IndexType inode=0; inode < numNodes; ++inode )
        {
          EXPECT_DOUBLE_EQ( x1[ inode ], MAGIC_NUM );
          EXPECT_DOUBLE_EQ( x2[ inode ], MAGIC_NUM );
        } // END for all nodes

      } // END for all ndims

      // clean up
      delete m1;
      m1 = AXOM_NULLPTR;

      delete m2;
      m2 = AXOM_NULLPTR;
    }
    // END SCOPE

    EXPECT_TRUE( blueprint::validRootGroup( m1grp ) );
    EXPECT_TRUE( blueprint::validRootGroup( m2grp ) );

  } // END for all dimensions

}

#endif /* MINT_USE_SIDRE */

//------------------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
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

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

#include "mint/blueprint.hpp"       // for blueprint functions
#include "mint/CellTypes.hpp"       // for CellTypes enum definition
#include "mint/ParticleMesh.hpp"    // for ParticleMesh
#include "mint/UniformMesh.hpp"     // for UniformMesh

#include "slic/slic.hpp"            // for slic macros

// Sidre includes
#ifdef MINT_USE_SIDRE
#include "sidre/sidre.hpp"
namespace sidre = axom::sidre;
#endif

#include "gtest/gtest.h"            // for gtest macros

using namespace axom::mint;

// globals
const char* IGNORE_OUTPUT = ".*";

//------------------------------------------------------------------------------
//  HELPER METHODS
//------------------------------------------------------------------------------
namespace
{

//------------------------------------------------------------------------------
#ifdef MINT_USE_SIDRE

void check_sidre_group( sidre::Group* root_group,
                        int expected_dimension,
                        const double* expected_origin,
                        const double* expected_spacing,
                        const int64* expected_extent )
{
  EXPECT_TRUE( blueprint::validRootGroup( root_group ) );

  const sidre::Group* topology = blueprint::getTopologyGroup( root_group );
  EXPECT_TRUE( blueprint::validTopologyGroup( topology ) );
  EXPECT_EQ( topology->getParent()->getNumGroups(), 1 );
  EXPECT_EQ( topology->getParent()->getNumViews(), 0 );

  const sidre::Group* coordset =
    blueprint::getCoordsetGroup( root_group, topology );
  EXPECT_TRUE( blueprint::validCoordsetGroup( coordset ) );
  EXPECT_EQ( coordset->getParent()->getNumGroups(), 1 );
  EXPECT_EQ( coordset->getParent()->getNumViews(), 0 );

  int mesh_type = UNDEFINED_MESH;
  int dimension = -1;
  blueprint::getMeshTypeAndDimension( mesh_type, dimension, root_group );
  EXPECT_EQ( mesh_type, STRUCTURED_UNIFORM_MESH );
  EXPECT_EQ( dimension, expected_dimension );

  double mesh_origin[ 3 ];
  double mesh_spacing[ 3 ];
  int64 mesh_extent[ 6 ];
  blueprint::getUniformMesh( dimension, coordset, topology,
                             mesh_origin,
                             mesh_spacing,
                             mesh_extent );

  for ( int i=0 ; i < dimension ; ++i )
  {
    EXPECT_DOUBLE_EQ( expected_origin[ i ], mesh_origin[ i ] );
    EXPECT_DOUBLE_EQ( expected_spacing[ i ], mesh_spacing[ i ] );
    EXPECT_EQ( expected_extent[ i*2 ], mesh_extent[ i*2 ] );
    EXPECT_EQ( expected_extent[ i*2+1 ], mesh_extent[ i*2+1 ] );
  }

}

#endif

//------------------------------------------------------------------------------
void check_create_field( UniformMesh* m,
                         int association,
                         const std::string& name,
                         int numComponents )
{
  EXPECT_TRUE( m != AXOM_NULLPTR );
  EXPECT_TRUE( (association==NODE_CENTERED) ||
               (association==CELL_CENTERED)    );

  const Field* field = AXOM_NULLPTR;

  double* f = m->createField< double >( name, association, numComponents );
  EXPECT_TRUE( f != AXOM_NULLPTR );
  EXPECT_TRUE( m->hasField( name, association ) );

  field = m->getFieldData( association )->getField( name );

  IndexType expected_num_tuples = ( association==NODE_CENTERED ) ?
                                  m->getNumberOfNodes() : m->getNumberOfCells();

  EXPECT_EQ( expected_num_tuples, field->getNumTuples( ) );
  EXPECT_EQ( numComponents, field->getNumComponents() );
  EXPECT_EQ( f, Field::getDataPtr< double >( field ) );
}

//------------------------------------------------------------------------------
void check_constructor( UniformMesh* m,
                        int expected_dimension,
                        const double* expected_origin,
                        const double* expected_spacing,
                        const int64* expected_extent,
                        const IndexType* expected_dimensions )
{
  EXPECT_TRUE( m != AXOM_NULLPTR );
  EXPECT_TRUE( expected_dimension >= 1 && expected_dimension <= 3 );

  // check Mesh methods
  EXPECT_EQ( m->getDimension(), expected_dimension );
  EXPECT_EQ( m->getMeshType(), STRUCTURED_UNIFORM_MESH );
  EXPECT_FALSE( m->hasExplicitConnectivity() );
  EXPECT_FALSE( m->hasExplicitCoordinates() );
  EXPECT_FALSE( m->hasMixedCellTypes() );

  const int mesh_dimension = m->getDimension();
  CellType expected_cell_type = ( mesh_dimension==3 ) ? HEX :
                                ( ( mesh_dimension==2 ) ? QUAD : SEGMENT );
  EXPECT_EQ( m->getCellType(), expected_cell_type );

  const Extent* extent = m->getExtent();
  EXPECT_TRUE( extent != AXOM_NULLPTR );
  EXPECT_EQ( extent->getNumNodes(), m->getNumberOfNodes() );
  EXPECT_EQ( extent->getNumCells(), m->getNumberOfCells() );

  const double* origin = m->getOrigin();
  EXPECT_TRUE( origin != AXOM_NULLPTR );

  const double* spacing = m->getSpacing();
  EXPECT_TRUE( spacing != AXOM_NULLPTR );

  for ( int i=0 ; i < mesh_dimension ; ++i )
  {
    EXPECT_DOUBLE_EQ( expected_origin[ i ], origin[ i ] );
    EXPECT_DOUBLE_EQ( expected_spacing[ i ], spacing[ i ] );

    EXPECT_EQ( expected_extent[ i*2 ], extent->min( i ) );
    EXPECT_EQ( expected_extent[ i*2+1 ], extent->max( i ) );
    EXPECT_EQ( expected_dimensions[ i ], extent->size( i ) );

  }
}

}

//------------------------------------------------------------------------------
//  UNIT TESTS
//------------------------------------------------------------------------------
TEST( mint_mesh_uniform_mesh_DeathTest, invalid_construction )
{
  const double origin[]  = { 0.0, 0.0, 0.0 };
  const double h[]       = { 0.5, 0.5, 0.5 };
  const int64 ext[]      = { 0,4, 0,4, 0,4 };
  const double lo[]      = { 0.0, 0.0, 0.0 };
  const double hi[]      = { 2.0, 2.0, 2.0 };
  const IndexType Ni     = 5;
  const IndexType Nj     = 5;
  const IndexType Nk     = 5;

  // check 1st native constructor
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(4,origin,h,ext),
                             IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(0,origin, h, ext),
                             IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(3,AXOM_NULLPTR, h, ext),
                             IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(3,origin,AXOM_NULLPTR,ext),
                             IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(3,origin,h,AXOM_NULLPTR),
                             IGNORE_OUTPUT );

  // check 2nd native constructor
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(4,ext,lo,hi),
                             IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(0,ext,lo,hi),
                             IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(3,AXOM_NULLPTR,lo,hi),
                             IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(3,ext,AXOM_NULLPTR,hi),
                             IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(3,ext,lo,AXOM_NULLPTR),
                             IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(3,ext,hi,lo),
                             IGNORE_OUTPUT );

  // check 3rd native constructor
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(4,lo,hi,Ni,Nj,Nk),
                             IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(0,lo,hi,Ni,Nj,Nk),
                             IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(3,AXOM_NULLPTR,hi,Ni,Nj,Nk),
                             IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(3,lo,AXOM_NULLPTR,Ni,Nj,Nk),
                             IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(3,hi,lo,Ni,Nj,Nk),
                             IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(1,hi,lo,-1),
                             IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(2,hi,lo, Ni,-1),
                             IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(3,hi,lo,Ni,Nj,-1),
                             IGNORE_OUTPUT );

#ifdef MINT_USE_SIDRE

  sidre::DataStore ds;
  sidre::Group* root          = ds.getRoot();
  sidre::Group* valid_group   = root->createGroup( "mesh" );
  sidre::Group* particle_mesh = root->createGroup( "particle_mesh" );
  ParticleMesh( 3, 10, particle_mesh );

  // check pull constructor
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh( root, "" ), IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh( particle_mesh, "" ), IGNORE_OUTPUT );


  // check 1st sidre constructor
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(3,lo,hi,ext,particle_mesh),
                             IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(3,lo,hi,ext,AXOM_NULLPTR),
                             IGNORE_OUTPUT );

  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(3,lo,hi,AXOM_NULLPTR,valid_group),
                             IGNORE_OUTPUT );
  EXPECT_EQ( valid_group->getNumGroups(), 0 );
  EXPECT_EQ( valid_group->getNumViews(), 0 );

  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(3,lo,AXOM_NULLPTR,ext,valid_group),
                             IGNORE_OUTPUT );
  EXPECT_EQ( valid_group->getNumGroups(), 0 );
  EXPECT_EQ( valid_group->getNumViews(), 0 );

  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(3,AXOM_NULLPTR,hi,ext,valid_group),
                             IGNORE_OUTPUT );
  EXPECT_EQ( valid_group->getNumGroups(), 0 );
  EXPECT_EQ( valid_group->getNumViews(), 0 );

  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(42,lo,hi,ext,valid_group),
                             IGNORE_OUTPUT );
  EXPECT_EQ( valid_group->getNumGroups(), 0 );
  EXPECT_EQ( valid_group->getNumViews(), 0 );

  // check 2nd sidre constructor
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(3,lo,hi,particle_mesh, Ni, Nj, Nk),
                             IGNORE_OUTPUT );

  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(3,lo,hi,AXOM_NULLPTR, Ni, Nj, Nk),
                             IGNORE_OUTPUT );
  EXPECT_EQ( valid_group->getNumGroups(), 0 );
  EXPECT_EQ( valid_group->getNumViews(), 0 );

  EXPECT_DEATH_IF_SUPPORTED(
    UniformMesh(3,lo,AXOM_NULLPTR,valid_group,Ni,Nj,Nk),
    IGNORE_OUTPUT );
  EXPECT_EQ( valid_group->getNumGroups(), 0 );
  EXPECT_EQ( valid_group->getNumViews(), 0 );

  EXPECT_DEATH_IF_SUPPORTED(
    UniformMesh(3,AXOM_NULLPTR,hi,valid_group,Ni,Nj,Nk),
    IGNORE_OUTPUT );
  EXPECT_EQ( valid_group->getNumGroups(), 0 );
  EXPECT_EQ( valid_group->getNumViews(), 0 );

  EXPECT_DEATH_IF_SUPPORTED(
    UniformMesh(42,lo,hi,valid_group,Ni,Nj,Nk),
    IGNORE_OUTPUT );
  EXPECT_EQ( valid_group->getNumGroups(), 0 );
  EXPECT_EQ( valid_group->getNumViews(), 0 );

  EXPECT_DEATH_IF_SUPPORTED(
    UniformMesh(3,lo,hi,valid_group,Ni,Nj,-1),
    IGNORE_OUTPUT );
  EXPECT_EQ( valid_group->getNumGroups(), 0 );
  EXPECT_EQ( valid_group->getNumViews(), 0 );

  EXPECT_DEATH_IF_SUPPORTED(
    UniformMesh(3,lo,hi,valid_group,Ni,-1,Nk),
    IGNORE_OUTPUT );
  EXPECT_EQ( valid_group->getNumGroups(), 0 );
  EXPECT_EQ( valid_group->getNumViews(), 0 );

  EXPECT_DEATH_IF_SUPPORTED(
    UniformMesh(3,lo,hi,valid_group,-1,Nj,Nk),
    IGNORE_OUTPUT );
  EXPECT_EQ( valid_group->getNumGroups(), 0 );
  EXPECT_EQ( valid_group->getNumViews(), 0 );

#endif
}

//------------------------------------------------------------------------------
TEST( mint_mesh_uniform_mesh, invalid_operations )
{
  const double origin[]  = {0.0, 0.0, 0.0};
  const double h[]       = { 0.5, 0.5, 0.5 };
  const int64 ext[]      = { 0,4, 0,4, 0,4 };

  UniformMesh m( 3, origin, h, ext );
  EXPECT_DEATH_IF_SUPPORTED( m.getCoordinateArray( 0 ), IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( m.getCoordinateArray( 1 ), IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( m.getCoordinateArray( 2 ), IGNORE_OUTPUT );
}

//------------------------------------------------------------------------------
TEST( mint_mesh_uniform_mesh, native_constructor )
{
  constexpr int NDIMS = 3;

  const double lo[]     = { 0.0, 0.0, 0.0 };
  const double hi[]     = { 2.0, 2.0, 2.0 };
  const double h[]      = { 0.5, 0.5, 0.5 };
  const int64 ext[]     = { 0,4, 0,4, 0,4 };
  const IndexType N[]   = {   5,   5,   5 };

  for ( int idim=1 ; idim <= NDIMS ; ++idim )
  {
    UniformMesh* m1 = new UniformMesh( idim, lo, h, ext );
    check_constructor( m1, idim, lo, h, ext, N );
    EXPECT_FALSE( m1->hasSidreGroup() );

    check_create_field( m1, NODE_CENTERED, "n1", 3 );
    check_create_field( m1, CELL_CENTERED, "c1", 2 );

    UniformMesh* m2 = new UniformMesh( idim, ext, lo, hi );
    check_constructor( m2, idim, lo, h, ext, N );
    EXPECT_FALSE( m2->hasSidreGroup() );

    check_create_field( m2, NODE_CENTERED, "n1", 3 );
    check_create_field( m2, CELL_CENTERED, "c1", 2 );

    UniformMesh* m3 = AXOM_NULLPTR;
    switch ( idim )
    {
    case 1:
      m3 = new UniformMesh( idim, lo, hi, N[0] );
      break;
    case 2:
      m3 = new UniformMesh( idim, lo, hi, N[0], N[1] );
      break;
    default:
      EXPECT_TRUE( idim==3 );
      m3 = new UniformMesh( idim, lo, hi, N[0], N[1], N[2] );
    }
    check_constructor( m3, idim, lo, h, ext, N );
    EXPECT_FALSE( m3->hasSidreGroup() );

    check_create_field( m3, NODE_CENTERED, "n1", 3 );
    check_create_field( m3, CELL_CENTERED, "c1", 2 );

    delete m1;
    delete m2;
    delete m3;
    m1 = m2 = m3 = AXOM_NULLPTR;
  } // END for all dimensions

}

#ifdef MINT_USE_SIDRE

//------------------------------------------------------------------------------
TEST( mint_mesh_uniform_mesh, sidre_constructor )
{
  constexpr int NDIMS = 3;
  const double lo[]     = { 0.0, 0.0, 0.0 };
  const double hi[]     = { 2.0, 2.0, 2.0 };
  const double h[]      = { 0.5, 0.5, 0.5 };
  const int64 ext[]     = { 0,4, 0,4, 0,4 };
  const IndexType N[]   = {   5,   5,   5 };

  for ( int idim=1 ; idim <= NDIMS ; ++idim )
  {

    // STEP 0: create a datastore to store two meshes
    sidre::DataStore ds;
    sidre::Group* root     = ds.getRoot();
    sidre::Group* mesh1grp = root->createGroup( "mesh_1" );
    sidre::Group* mesh2grp = root->createGroup( "mesh_2" );

    // STEP 1: create 2 identical uniform mesh in Sidre using the two different
    // flavors of the Sidre constructor
    // BEGIN SCOPE
    {
      UniformMesh m1( idim, lo, hi, ext, mesh1grp );
      check_constructor( &m1, idim, lo, h, ext, N );
      check_create_field( &m1, NODE_CENTERED, "n1", 3 );
      check_create_field( &m1, CELL_CENTERED, "c1", 1 );
      EXPECT_TRUE( m1.hasSidreGroup( ) );

      UniformMesh* m2 = AXOM_NULLPTR;
      switch ( idim )
      {
      case 1:
        m2 = new UniformMesh( idim, lo, hi, mesh2grp, N[0] );
        break;
      case 2:
        m2 = new UniformMesh( idim, lo, hi, mesh2grp, N[0], N[1] );
        break;
      default:
        EXPECT_TRUE( idim==3 );
        m2 = new UniformMesh( idim, lo, hi, mesh2grp, N[0], N[1], N[2] );
      } // END switch

      EXPECT_TRUE( m2->hasSidreGroup( ) );
      check_constructor( m2, idim, lo, h, ext, N );
      check_create_field( m2, NODE_CENTERED, "n1", 3 );
      check_create_field( m2, CELL_CENTERED, "c1", 1 );

      delete m2;
      m2 = AXOM_NULLPTR;
    }
    // END SCOPE

    // STEP 2: pull the mesh from the sidre groups, and check with expected
    {
      UniformMesh* m1 = new UniformMesh( mesh1grp );
      check_constructor( m1, idim, lo, h, ext, N );
      EXPECT_TRUE( m1->hasSidreGroup() );
      EXPECT_TRUE( m1->hasField( "n1", NODE_CENTERED) );
      EXPECT_TRUE( m1->hasField( "c1", CELL_CENTERED) );

      const FieldData* node_data1 = m1->getFieldData( NODE_CENTERED );
      EXPECT_TRUE( node_data1 != AXOM_NULLPTR );
      const FieldData* cell_data1 = m1->getFieldData( CELL_CENTERED );
      EXPECT_TRUE( cell_data1 != AXOM_NULLPTR );

      EXPECT_EQ( 3, node_data1->getField("n1")->getNumComponents() );
      EXPECT_EQ( 1, cell_data1->getField("c1")->getNumComponents() );

      UniformMesh* m2 = new UniformMesh( mesh2grp );
      check_constructor( m2, idim, lo, h, ext, N );
      EXPECT_TRUE( m2->hasSidreGroup() );
      EXPECT_TRUE( m2->hasField( "n1", NODE_CENTERED ) );
      EXPECT_TRUE( m2->hasField( "c1", CELL_CENTERED ) );

      const FieldData* node_data2 = m2->getFieldData( NODE_CENTERED );
      EXPECT_TRUE( node_data2 != AXOM_NULLPTR );
      const FieldData* cell_data2 = m2->getFieldData( CELL_CENTERED );
      EXPECT_TRUE( cell_data2 != AXOM_NULLPTR );

      EXPECT_EQ( 3, node_data1->getField("n1")->getNumComponents() );
      EXPECT_EQ( 1, cell_data1->getField("c1")->getNumComponents() );

      delete m1;
      delete m2;
      m1 = m2 = AXOM_NULLPTR;
    }

    // STEP 3: ensure the data is persistent in sidre
    check_sidre_group( mesh1grp, idim, lo, h, ext );
    check_sidre_group( mesh2grp, idim, lo, h, ext );

  }

}

#endif

//------------------------------------------------------------------------------
TEST( mint_mesh_uniform_mesh, check_evaluate_coordinate )
{
  constexpr int NDIMS   = 3;
  const double lo[]     = { 0.0, 0.0, 0.0 };
  const double hi[]     = { 2.0, 2.0, 2.0 };
  const double h[]      = { 0.5, 0.5, 0.5 };
  const int64 ext[]     = { -2,2, -2,2, -2,2 };
  const IndexType N[]   = {   5,   5,   5 };

  for ( int idim=1 ; idim <= NDIMS ; ++idim )
  {
    UniformMesh m( idim, ext, lo, hi );
    check_constructor( &m, idim, lo, h, ext, N );
    EXPECT_FALSE( m.hasSidreGroup() );

    switch ( idim )
    {
    case 1:
    {
      const IndexType Ni = m.getNumberOfNodesAlongDim( I_DIRECTION );
      EXPECT_EQ( Ni, N[ I_DIRECTION ] );

      for ( IndexType i=0 ; i < Ni ; ++i )
      {
        const double expected_x = lo[ X_COORDINATE ] + i * h[ X_COORDINATE ];
        const double x = m.evaluateCoordinate( i, I_DIRECTION );
        EXPECT_DOUBLE_EQ( expected_x, x );
      } // END for all i

    }   // END 1D
    break;
    case 2:
    {

      const IndexType Ni = m.getNumberOfNodesAlongDim( I_DIRECTION );
      const IndexType Nj = m.getNumberOfNodesAlongDim( J_DIRECTION );
      EXPECT_EQ( Ni, N[ I_DIRECTION ] );
      EXPECT_EQ( Nj, N[ J_DIRECTION ] );

      for ( IndexType j=0 ; j < Nj ; ++j )
      {
        for ( IndexType i=0 ; i < Ni ; ++i )
        {
          const double expected_x = lo[ X_COORDINATE ] + i * h[ X_COORDINATE ];
          const double expected_y = lo[ Y_COORDINATE ] + j * h[ Y_COORDINATE ];

          const double x = m.evaluateCoordinate( i, I_DIRECTION );
          const double y = m.evaluateCoordinate( j, J_DIRECTION );

          EXPECT_DOUBLE_EQ( expected_x, x );
          EXPECT_DOUBLE_EQ( expected_y, y );
        } // END for all i
      } // END for all j

    }   // END 2D
    break;
    default:
      EXPECT_EQ( idim, 3 );
      {

        const IndexType Ni = m.getNumberOfNodesAlongDim( I_DIRECTION );
        const IndexType Nj = m.getNumberOfNodesAlongDim( J_DIRECTION );
        const IndexType Nk = m.getNumberOfNodesAlongDim( K_DIRECTION );
        EXPECT_EQ( Ni, N[ I_DIRECTION ] );
        EXPECT_EQ( Nj, N[ J_DIRECTION ] );
        EXPECT_EQ( Nk, N[ K_DIRECTION ] );

        for ( IndexType k=0 ; k < Nk ; ++k )
        {
          for ( IndexType j=0 ; j < Nj ; ++j )
          {
            for ( IndexType i=0 ; i < Ni ; ++i )
            {
              const double expected_x = lo[ X_COORDINATE ] + i*
                                        h[ X_COORDINATE ];
              const double expected_y = lo[ Y_COORDINATE ] + j*
                                        h[ Y_COORDINATE ];
              const double expected_z = lo[ Z_COORDINATE ] + k*
                                        h[ Z_COORDINATE ];

              const double x = m.evaluateCoordinate( i, I_DIRECTION );
              const double y = m.evaluateCoordinate( j, J_DIRECTION );
              const double z = m.evaluateCoordinate( k, K_DIRECTION );

              EXPECT_DOUBLE_EQ( expected_x, x );
              EXPECT_DOUBLE_EQ( expected_y, y );
              EXPECT_DOUBLE_EQ( expected_z, z );
            } // END for all i
          } // END for all j
        } // END for all k

      } // END 3D

    } // END switch

  }

}

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

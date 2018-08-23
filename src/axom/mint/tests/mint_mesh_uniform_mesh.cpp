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
#include "axom/mint/config.hpp"               // for compile-time type
                                              // definitions

#include "axom/mint/mesh/blueprint.hpp"       // for blueprint functions
#include "axom/mint/mesh/CellTypes.hpp"       // for CellTypes enum definition
#include "axom/mint/mesh/ParticleMesh.hpp"    // for ParticleMesh
#include "axom/mint/mesh/UniformMesh.hpp"     // for UniformMesh
#include "StructuredMesh_helpers.hpp"         // for StructuredMesh test helpers


#include "axom/slic/interface/slic.hpp"       // for slic macros

// Sidre includes
#ifdef AXOM_MINT_USE_SIDRE
#include "axom/sidre/core/sidre.hpp"          // for sidre classes
namespace sidre = axom::sidre;
#endif

#include "gtest/gtest.h"                      // for gtest macros

using namespace axom::mint;

// globals
const char* IGNORE_OUTPUT = ".*";

//------------------------------------------------------------------------------
//  HELPER METHODS
//------------------------------------------------------------------------------
namespace
{

//------------------------------------------------------------------------------
#ifdef AXOM_MINT_USE_SIDRE

void check_sidre_group( sidre::Group* root_group,
                        int expected_dimension,
                        const double* expected_origin,
                        const double* expected_spacing )
{
  EXPECT_TRUE( blueprint::isValidRootGroup( root_group ) );

  const sidre::Group* topology = blueprint::getTopologyGroup( root_group );
  EXPECT_TRUE( blueprint::isValidTopologyGroup( topology ) );
  EXPECT_EQ( topology->getParent()->getNumGroups(), 1 );
  EXPECT_EQ( topology->getParent()->getNumViews(), 0 );

  const sidre::Group* coordset =
    blueprint::getCoordsetGroup( root_group, topology );
  EXPECT_TRUE( blueprint::isValidCoordsetGroup( coordset ) );
  EXPECT_EQ( coordset->getParent()->getNumGroups(), 1 );
  EXPECT_EQ( coordset->getParent()->getNumViews(), 0 );

  int mesh_type = UNDEFINED_MESH;
  int dimension = -1;
  blueprint::getMeshTypeAndDimension( mesh_type, dimension, root_group );
  EXPECT_EQ( mesh_type, STRUCTURED_UNIFORM_MESH );
  EXPECT_EQ( dimension, expected_dimension );

  double mesh_origin[ 3 ];
  double mesh_spacing[ 3 ];
  blueprint::getUniformMesh( dimension, mesh_origin, mesh_spacing, coordset );

  for ( int i=0 ; i < dimension ; ++i )
  {
    EXPECT_DOUBLE_EQ( expected_origin[ i ], mesh_origin[ i ] );
    EXPECT_DOUBLE_EQ( expected_spacing[ i ], mesh_spacing[ i ] );
  }

}

#endif

}

//------------------------------------------------------------------------------
//  UNIT TESTS
//------------------------------------------------------------------------------
TEST( mint_mesh_uniform_mesh_DeathTest, invalid_construction )
{
  const double origin[]  = { 0.0, 0.0, 0.0 };
  const double h[]       = { 0.5, 0.5, 0.5 };
  const double lo[]      = { 0.0, 0.0, 0.0 };
  const double hi[]      = { 2.0, 2.0, 2.0 };
  const IndexType Ni     = 5;
  const IndexType Nj     = 5;
  const IndexType Nk     = 5;
  const IndexType N[]    = { Ni, Nj, Nk };

  // check 1st native constructor
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(4,origin,h,N),
                             IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(0,origin, h, N),
                             IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(3,nullptr, h, N),
                             IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(3,origin,nullptr,N),
                             IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(3,origin,h,nullptr),
                             IGNORE_OUTPUT );

  // check 2nd native constructor
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(4,N,lo,hi),
                             IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(0,N,lo,hi),
                             IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(3,nullptr,lo,hi),
                             IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(3,N,nullptr,hi),
                             IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(3,N,lo,nullptr),
                             IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(3,N,hi,lo),
                             IGNORE_OUTPUT );

  // check 3rd native constructor
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(nullptr,hi,Ni,Nj,Nk),
                             IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(lo,nullptr,Ni,Nj,Nk),
                             IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(lo, hi,-1),
                             IGNORE_OUTPUT );


#ifdef AXOM_MINT_USE_SIDRE

  sidre::DataStore ds;
  sidre::Group* root          = ds.getRoot();
  sidre::Group* valid_group   = root->createGroup( "mesh" );
  sidre::Group* particle_mesh = root->createGroup( "particle_mesh" );
  ParticleMesh( 3, 10, particle_mesh );

  // check pull constructor
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh( root, "" ), IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh( particle_mesh, "" ), IGNORE_OUTPUT );


  // check 1st sidre constructor
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(3,lo,hi,N,particle_mesh),
                             IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(3,lo,hi,N,nullptr),
                             IGNORE_OUTPUT );

  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(3,lo,hi,nullptr,valid_group),
                             IGNORE_OUTPUT );
  EXPECT_EQ( valid_group->getNumGroups(), 0 );
  EXPECT_EQ( valid_group->getNumViews(), 0 );

  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(3,lo,nullptr,N,valid_group),
                             IGNORE_OUTPUT );
  EXPECT_EQ( valid_group->getNumGroups(), 0 );
  EXPECT_EQ( valid_group->getNumViews(), 0 );

  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(3,nullptr,hi,N,valid_group),
                             IGNORE_OUTPUT );
  EXPECT_EQ( valid_group->getNumGroups(), 0 );
  EXPECT_EQ( valid_group->getNumViews(), 0 );

  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(42,lo,hi,N,valid_group),
                             IGNORE_OUTPUT );
  EXPECT_EQ( valid_group->getNumGroups(), 0 );
  EXPECT_EQ( valid_group->getNumViews(), 0 );

  // check 2nd sidre constructor
  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(lo,hi,particle_mesh, Ni, Nj, Nk),
                             IGNORE_OUTPUT );

  EXPECT_DEATH_IF_SUPPORTED( UniformMesh(lo,hi,nullptr, Ni, Nj, Nk),
                             IGNORE_OUTPUT );
  EXPECT_EQ( valid_group->getNumGroups(), 0 );
  EXPECT_EQ( valid_group->getNumViews(), 0 );

  EXPECT_DEATH_IF_SUPPORTED(
    UniformMesh(lo,nullptr,valid_group,Ni,Nj,Nk),
    IGNORE_OUTPUT );
  EXPECT_EQ( valid_group->getNumGroups(), 0 );
  EXPECT_EQ( valid_group->getNumViews(), 0 );

  EXPECT_DEATH_IF_SUPPORTED(
    UniformMesh(nullptr,hi,valid_group,Ni,Nj,Nk),
    IGNORE_OUTPUT );
  EXPECT_EQ( valid_group->getNumGroups(), 0 );
  EXPECT_EQ( valid_group->getNumViews(), 0 );

  EXPECT_DEATH_IF_SUPPORTED(
    UniformMesh(lo,hi,valid_group,-1),
    IGNORE_OUTPUT );
  EXPECT_EQ( valid_group->getNumGroups(), 0 );
  EXPECT_EQ( valid_group->getNumViews(), 0 );

#endif
}

//------------------------------------------------------------------------------
TEST( mint_mesh_uniform_mesh_DeathTest, invalid_operations )
{
  const double origin[]  = {0.0, 0.0, 0.0};
  const double h[]       = { 0.5, 0.5, 0.5 };
  const IndexType N[]    = { 4, 4, 4 };

  UniformMesh m( 3, origin, h, N );
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
  const IndexType N[]   = {   5,   5,   5 };
  const int64 node_ext[] = { -55, -51, 10, 14, 0, 4 };

  for ( int idim=1 ; idim <= NDIMS ; ++idim )
  {
    UniformMesh* m = new UniformMesh( idim, lo, h, N );
    internal::check_constructor( m, idim, lo, h, N );
    EXPECT_FALSE( m->hasSidreGroup() );

    internal::check_create_fields( m );
    internal::check_fields( m, false );
    m->setNodeExtent( idim, node_ext );
    internal::check_node_extent( m, node_ext );
    delete m;

    m = new UniformMesh( idim, N, lo, hi );
    internal::check_constructor( m, idim, lo, h, N );
    EXPECT_FALSE( m->hasSidreGroup() );

    internal::check_create_fields( m );
    internal::check_fields( m, false );
    m->setNodeExtent( idim, node_ext );
    internal::check_node_extent( m, node_ext );
    delete m;

    switch ( idim )
    {
    case 1:
      m = new UniformMesh( lo, hi, N[0] );
      break;
    case 2:
      m = new UniformMesh( lo, hi, N[0], N[1] );
      break;
    default:
      EXPECT_TRUE( idim==3 );
      m = new UniformMesh( lo, hi, N[0], N[1], N[2] );
    }

    internal::check_constructor( m, idim, lo, h, N );
    EXPECT_FALSE( m->hasSidreGroup() );

    internal::check_create_fields( m );
    internal::check_fields( m, false );
    m->setNodeExtent( idim, node_ext );
    internal::check_node_extent( m, node_ext );
    delete m;
  } // END for all dimensions

}

#ifdef AXOM_MINT_USE_SIDRE

//------------------------------------------------------------------------------
TEST( mint_mesh_uniform_mesh, sidre_constructor )
{
  constexpr int NDIMS = 3;
  const double lo[]     = { 0.0, 0.0, 0.0 };
  const double hi[]     = { 2.0, 2.0, 2.0 };
  const double h[]      = { 0.5, 0.5, 0.5 };
  const IndexType N[]   = {   5,   5,   5 };
  const int64 node_ext[] = { -55, -51, 10, 14, 0, 4 };

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
      UniformMesh* m = new UniformMesh( idim, lo, hi, N, mesh1grp );
      EXPECT_TRUE( m->hasSidreGroup() );
      internal::check_constructor( m, idim, lo, h, N );
      internal::check_create_fields( m );
      m->setNodeExtent( idim, node_ext );
      internal::check_node_extent( m, node_ext );

      delete m;
      m = nullptr;

      switch ( idim )
      {
      case 1:
        m = new UniformMesh( lo, hi, mesh2grp, N[0] );
        break;
      case 2:
        m = new UniformMesh( lo, hi, mesh2grp, N[0], N[1] );
        break;
      default:
        EXPECT_TRUE( idim==3 );
        m = new UniformMesh( lo, hi, mesh2grp, N[0], N[1], N[2] );
      } // END switch

      EXPECT_TRUE( m->hasSidreGroup() );
      internal::check_constructor( m, idim, lo, h, N );
      internal::check_create_fields( m );
      m->setNodeExtent( idim, node_ext );
      internal::check_node_extent( m, node_ext );

      delete m;
    }
    // END SCOPE

    // STEP 2: pull the mesh from the sidre groups, and check with expected
    {
      UniformMesh* m = new UniformMesh( mesh1grp );
      internal::check_constructor( m, idim, lo, h, N );
      EXPECT_TRUE( m->hasSidreGroup() );
      internal::check_fields( m, true );
      internal::check_node_extent( m, node_ext );
      delete m;

      m = new UniformMesh( mesh2grp );
      internal::check_constructor( m, idim, lo, h, N );
      internal::check_fields( m, true );
      internal::check_node_extent( m, node_ext );
      delete m;
    }

    // STEP 3: ensure the data is persistent in sidre
    check_sidre_group( mesh1grp, idim, lo, h );
    check_sidre_group( mesh2grp, idim, lo, h );
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
  const IndexType N[]   = {   5,   5,   5 };

  for ( int idim=1 ; idim <= NDIMS ; ++idim )
  {
    UniformMesh m( idim, N, lo, hi );
    internal::check_constructor( &m, idim, lo, h, N );
    EXPECT_FALSE( m.hasSidreGroup() );

    switch ( idim )
    {
    case 1:
    {
      const IndexType Ni = m.getNodeDimension( I_DIRECTION );
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

      const IndexType Ni = m.getNodeDimension( I_DIRECTION );
      const IndexType Nj = m.getNodeDimension( J_DIRECTION );
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

        const IndexType Ni = m.getNodeDimension( I_DIRECTION );
        const IndexType Nj = m.getNodeDimension( J_DIRECTION );
        const IndexType Nk = m.getNodeDimension( K_DIRECTION );
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

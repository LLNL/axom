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
#include "mint/config.hpp"

#include "mint/Mesh.hpp"            // for Mesh base class
#include "mint/ParticleMesh.hpp"    // for ParticleMesh definition

// Sidre includes
#ifdef MINT_USE_SIDRE
#include "sidre/sidre.hpp"
namespace sidre = axom::sidre;
#endif

#include "gtest/gtest.h"        // for gtest macros

//namespace aliases
namespace mint = axom::mint;

//------------------------------------------------------------------------------
//  UNIT TESTS
//------------------------------------------------------------------------------
#ifdef MINT_USE_SIDRE

TEST( mint_mesh, get_particle_mesh_from_sidre )
{
  constexpr int DIMENSION = 3;
  constexpr mint::IndexType NUM_PARTICLES = 10;
  constexpr double MAGIC_NUMBER = 42.0;

  // STEP 0: get empty Sidre group where to store the mesh
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  // STEP 1: populate the group with a particle mesh (positions + fields)
  mint::ParticleMesh particles( DIMENSION, 9, 9, NUM_PARTICLES, root );
  double* x = particles.getParticlePositions( mint::X_COORDINATE );
  double* y = particles.getParticlePositions( mint::Y_COORDINATE );
  double* z = particles.getParticlePositions( mint::Z_COORDINATE );

  double* phi = particles.createField< double >( "phi", mint::NODE_CENTERED );
  mint::IndexType* id =
      particles.createField< mint::IndexType >( "id", mint::NODE_CENTERED );

  for ( mint::IndexType ipart=0; ipart < NUM_PARTICLES; ++ipart )
  {
    const double val = static_cast< double >( ipart+ 1 );
    x[ ipart ]       = y[ ipart ] = z[ ipart ] = val*val;
    phi[ ipart ]     = MAGIC_NUMBER;
    id[ ipart ]      = ipart;
  }

  // STEP 2: get a mesh object from the group
  mint::Mesh* m = mint::Mesh::getMesh( root );

  // STEP 3: test the object
  EXPECT_EQ( m->getMeshType(), mint::PARTICLE_MESH );
  EXPECT_EQ( m->getBlockId(), 9 );
  EXPECT_EQ( m->getPartitionId(), 9 );
  EXPECT_EQ( m->getDimension(), DIMENSION );
  EXPECT_EQ( m->getNumberOfNodes(), NUM_PARTICLES );
  EXPECT_TRUE( m->hasField( "phi", mint::NODE_CENTERED ) );
  EXPECT_TRUE( m->hasField( "id", mint::NODE_CENTERED ) );
  EXPECT_EQ( m->getCoordinateArray( mint::X_COORDINATE), x );
  EXPECT_EQ( m->getCoordinateArray( mint::Y_COORDINATE), y );
  EXPECT_EQ( m->getCoordinateArray( mint::Z_COORDINATE), z );

  double* phi_test = m->getFieldPtr< double >( "phi", mint::NODE_CENTERED );
  EXPECT_EQ( phi, phi_test );

  mint::IndexType* id_test  =
      m->getFieldPtr< mint::IndexType >( "id", mint::NODE_CENTERED );
  EXPECT_EQ( id, id_test );

  // STEP 4: de-allocate
  delete m;
  m = AXOM_NULLPTR;
}

#endif

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



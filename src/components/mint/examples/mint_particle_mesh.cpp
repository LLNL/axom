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

/*!
 * \file
 *
 * \brief Illustrates how to construct and use a ParticleMesh to perform
 *  operations on a set of particles.
 */

// Mint includes
#include "mint/config.hpp"
#include "mint/ParticleMesh.hpp"
#include "mint/vtk_utils.hpp"

// C/C++ includes
#include <cstdlib>               // for rand()
#include <ctime>                 // for time()

// namespace aliases
namespace mint = axom::mint;

inline double random_double( double min, double max )
{
  const double t  = ( double( rand() ) / RAND_MAX );
  const double dx = max-min;
  return (min + t*dx );
}

//------------------------------------------------------------------------------
int main( int AXOM_NOT_USED(argc), char** AXOM_NOT_USED(argv) )
{
  using int64 = mint::IndexType;
  const mint::IndexType NUM_PARTICLES = 100;
  const int DIMENSION = 3;

  const double HI  = -10.0;
  const double LO  = 10.0;
  const double VLO = 0.0;
  const double VHI = 1.0;

  srand( time(0) );

  // STEP 0: create the ParticleMesh
  mint::ParticleMesh particles( DIMENSION, NUM_PARTICLES );

  // STEP 1: Add fields to the Particles
  double* vx = particles.createField< double >( "vx", mint::NODE_CENTERED );
  double* vy = particles.createField< double >( "vy", mint::NODE_CENTERED );
  double* vz = particles.createField< double >( "vz", mint::NODE_CENTERED );
  int64* id  = particles.createField< int64 >( "id", mint::NODE_CENTERED );

  // STEP 2: grab handle to the particle position arrays
  double* px = particles.getCoordinateArray( mint::X_COORDINATE );
  double* py = particles.getCoordinateArray( mint::Y_COORDINATE );
  double* pz = particles.getCoordinateArray( mint::Z_COORDINATE );

  // STEP 3: loop over the particle data
  const int64 numParticles = particles.getNumberOfNodes();
  for ( int64 i=0; i < numParticles; ++i )
  {

    px[ i ] = random_double( LO, HI );
    py[ i ] = random_double( LO, HI );
    pz[ i ] = random_double( LO, HI );

    vx[ i ] = random_double( VLO, VHI );
    vy[ i ] = random_double( VLO, VHI );
    vz[ i ] = random_double( VLO, VHI );

    id[ i ] = i;

  } // END

  // STEP 4: write the particle mesh in VTK format for visualization
  mint::write_vtk( &particles, "particles.vtk" );

  return 0;
}

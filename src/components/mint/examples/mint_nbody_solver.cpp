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
 * \brief A simple, naive, N-Body solver to illustrate the use of the
 *  ParticleMesh class
 *
 */

// Axom includes
#include "mint/ParticleMesh.hpp"   // for mint::ParticleMesh
#include "mint/vtk_utils.hpp"      // for mint::write_vtk()

#include "slic/slic.hpp"           // for SLIC macros
#include "slic/UnitTestLogger.hpp" // for slic logger

// Axom namespace aliases
namespace mint = axom::mint;
namespace slic = axom::slic;

// C/C++ includes
#include <cstdlib>  // for srand()/rand()
#include <cstring>  // for strcmp()
#include <ctime>    // for time()
#include <fstream>  // for std::ofstream
#include <sstream>  // for std::ostringstream

/*!
 * \brief Holds program arguments
 */
static struct
{
  double domain_min;
  double domain_max;
  double dt;
  int ncycles;
  int dumpFrequency;
  mint::IndexType numParticles;
} Arguments;

//------------------------------------------------------------------------------
// HELPER FUNCTION DEFINITIONS
//------------------------------------------------------------------------------
void showHelp();
void parse_arguments( int argc, char** argv );
void initialize( mint::ParticleMesh& particles );
void apply_forces( mint::ParticleMesh& particles, double dt );
inline double random_double( double min, double max );

//------------------------------------------------------------------------------
int main ( int argc, char** argv )
{
  // STEP 0: initialize logger & parse arguments
  srand( time(0) );
  slic::UnitTestLogger logger;
  parse_arguments( argc, argv );

  // STEP 1: construct particle mesh
  mint::ParticleMesh particles( 3, Arguments.numParticles );
  double* vx = particles.createField< double >( "vx", mint::NODE_CENTERED );
  double* vy = particles.createField< double >( "vy", mint::NODE_CENTERED );
  double* vz = particles.createField< double >( "vz", mint::NODE_CENTERED );

  initialize( particles );

  double* px = particles.getParticlePositions( mint::X_COORDINATE );
  double* py = particles.getParticlePositions( mint::Y_COORDINATE );
  double* pz = particles.getParticlePositions( mint::Z_COORDINATE );

  const mint::IndexType numParticles = particles.getNumParticles();

  // STEP 2: time march
  const double& dt = Arguments.dt;
  double t=0.0;
  for ( int cycle=0; cycle < Arguments.ncycles; ++cycle, t+=dt )
  {
    SLIC_INFO( "cycle=" << cycle << " t=" << t );

    apply_forces( particles, dt );

    // integrate positions
    for ( mint::IndexType i=0; i < numParticles; ++i )
    {
      px[ i ] += vx[ i ]*dt;
      py[ i ] += vy[ i ]*dt;
      pz[ i ] += vz[ i ]*dt;
    } // END for all particles

    if ( (cycle % Arguments.dumpFrequency) == 0 )
    {
      std::ostringstream oss;
      oss << "particles_" << cycle << ".vtk";
      std::string fileName = oss.str();
      mint::write_vtk( &particles, fileName );
      SLIC_INFO( "wrote " << fileName );
    }

  }

}

//------------------------------------------------------------------------------
// HELPER FUNCTION IMPLEMENTATION
//------------------------------------------------------------------------------
void apply_forces( mint::ParticleMesh& particles, double dt )
{
  const double ptiny = 1.e-9;

  const mint::IndexType numParticles = particles.getNumParticles();

  double* vx = particles.getFieldPtr< double >( "vx", mint::NODE_CENTERED );
  double* vy = particles.getFieldPtr< double >( "vy", mint::NODE_CENTERED );
  double* vz = particles.getFieldPtr< double >( "vz", mint::NODE_CENTERED );
  double* px = particles.getParticlePositions( mint::X_COORDINATE );
  double* py = particles.getParticlePositions( mint::Y_COORDINATE );
  double* pz = particles.getParticlePositions( mint::Z_COORDINATE );

  for ( mint::IndexType i=0; i < numParticles; ++i )
  {
    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;

    for ( mint::IndexType j=0; j < numParticles; ++j )
    {
      const double dx       = px[ j ] - px[ i ];
      const double dy       = py[ j ] - py[ i ];
      const double dz       = pz[ j ] - pz[ i ];
      const double d2       = dx*dx + dy*dy + dz*dz + ptiny;
      const double invdist  = 1.0f / sqrt( d2 );
      const double invdist3 = invdist * invdist * invdist;

      // add force contributions
      fx += dx*invdist3;
      fy += dy*invdist3;
      fz += dz*invdist3;
    }

    // udpate velocities
    vx[ i ] += dt*fx;
    vy[ i ] += dt*fy;
    vz[ i ] += dt*fz;
  }

}

//------------------------------------------------------------------------------
void initialize( mint::ParticleMesh& particles )
{
  const double LO = Arguments.domain_min;
  const double HI = Arguments.domain_max;

  mint::IndexType numParticles = particles.getNumParticles();

  double* vx = particles.getFieldPtr< double >( "vx", mint::NODE_CENTERED );
  double* vy = particles.getFieldPtr< double >( "vy", mint::NODE_CENTERED );
  double* vz = particles.getFieldPtr< double >( "vz", mint::NODE_CENTERED );
  double* px = particles.getParticlePositions( mint::X_COORDINATE );
  double* py = particles.getParticlePositions( mint::Y_COORDINATE );
  double* pz = particles.getParticlePositions( mint::Z_COORDINATE );

  for ( mint::IndexType i=0; i < numParticles; ++i )
  {
     px[ i ] = random_double( LO, HI );
     py[ i ] = random_double( LO, HI );
     pz[ i ] = random_double( LO, HI );

     vx[ i ] = random_double( LO, HI );
     vy[ i ] = random_double( LO, HI );
     vz[ i ] = random_double( LO, HI );
  }

}

//------------------------------------------------------------------------------
inline double random_double( double min, double max )
{
  const double t  = ( double( rand() ) / RAND_MAX );
  const double dx = max-min;
  return (min + t*dx );
}

//------------------------------------------------------------------------------
void showHelp()
{
  SLIC_INFO( "Usage: ./mint_nbody_solver_ex [options]" );
  SLIC_INFO( "--dt <dt> specifies the time increment at each cycle" );
  SLIC_INFO( "--ncycles <N> specifies the number of cycles to run" );
  SLIC_INFO( "--nparts <N> specifies number of particles in the system" );
  SLIC_INFO( "--dump <N> specifies the dump frequency" );
  SLIC_INFO( "--help prints this help information" );
}
//------------------------------------------------------------------------------
void parse_arguments( int argc, char** argv )
{
  Arguments.dt            = 0.01;
  Arguments.ncycles       = 25;
  Arguments.numParticles  = 1000;
  Arguments.domain_max    =  10;
  Arguments.domain_min    = -10;
  Arguments.dumpFrequency = 2;

  for ( int i=1; i < argc; ++i )
  {

    if ( strcmp( argv[i], "--dt")==0 )
    {
      Arguments.dt = atof( argv[ ++i ] );
    }
    else if ( strcmp( argv[i], "--ncycles" )==0 )
    {
      Arguments.ncycles = atoi( argv[ ++i ] );
    }
    else if ( strcmp( argv[i], "--nparts" )==0 )
    {
      Arguments.numParticles = atoi( argv[ ++i ] );
    }
    else if ( strcmp( argv[i], "--dump" )==0 )
    {
      Arguments.dumpFrequency = atoi( argv[ ++i ] );
    }
    else if ( strcmp( argv[i], "--help" )==0 )
    {
      showHelp( );
      exit( 0 );
    }

  }

}


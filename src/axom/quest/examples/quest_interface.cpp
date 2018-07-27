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
 * \file quest_interface.cpp
 *
 * \brief Simple example that exercises the quest interface for point
 *        containment and signed distance queries.
 */

#include "axom_utils/Utilities.hpp"

// Quest includes
#include "quest/quest.hpp"
#include "quest/signed_distance.hpp"

// SLIC includes
#include "slic/slic.hpp"

#ifdef AXOM_USE_MPI
  #include <mpi.h>

  #ifdef AXOM_USE_LUMBERJACK
    #include "slic/LumberjackStream.hpp"
  #else
    #include "slic/SynchronizedStream.hpp"
  #endif
#else
  #include "slic/GenericOutputStream.hpp"
#endif  // AXOM_USE_MPI


// C/C++ includes
#include <cmath>     // for signbit()
#include <cstdlib>   // for std::rand(), RAND_MAX
#include <ctime>     // for time
#include <iostream>  // for std::cout


using namespace axom;

typedef std::vector<double> CoordsVec;

//------------------------------------------------------------------------------
template < typename T >
T getRandomDouble( T low, T high )
{
  const T delta = high-low;
  const T c = static_cast< T >( std::rand() ) / static_cast< T >( RAND_MAX );
  return ( delta*c+low );
}

template < typename T >
T scaleAndOffset( T rangeMin, T rangeMax, T val)
{
  return rangeMin + (rangeMax - rangeMin) * val;
}


void outputMeshStats()
{
  // Obtain and log the mesh bounding box
  double bbMin[3], bbMax[3];
  quest::mesh_min_bounds(bbMin);
  quest::mesh_max_bounds(bbMax);
  SLIC_INFO(
    "Mesh bounding box: "
    << "{ lower: (" << bbMin[0] <<"," << bbMin[1] <<","<< bbMin[2] << ")"
    << "; upper: (" << bbMax[0] <<"," << bbMax[1] <<","<< bbMax[2] << ")}" );

  // Obtain and log the mesh center of mass
  double cm[3];
  quest::mesh_center_of_mass(cm);
  SLIC_INFO(
    "Mesh center of mass: "
    << "(" << cm[0] <<"," << cm[1] <<","<< cm[2] << ")"    );
}

void runQuestDistance(const std::string& fileName, const CoordsVec& points)
{
  constexpr int DIMENSION     = 3;
  constexpr int MAX_LEVELS    = 20;
  constexpr int MAX_OCCUPANCY = 25;

  // set signed distance parameters
  quest::signed_distance_set_dimension( DIMENSION );
  quest::signed_distance_set_max_levels( MAX_LEVELS );
  quest::signed_distance_set_max_occupancy( MAX_OCCUPANCY );

  // initialize the signed distance query
#ifdef AXOM_USE_MPI
  quest::signed_distance_init( fileName, MPI_COMM_WORLD );
#else
  quest::signed_distance_init( fileName );
#endif

  CoordsVec coords[3];
  {
    int nOrigPts = points.size()/3;

    double bbMin[3], bbMax[3];
    quest::signed_distance_get_mesh_bounds( bbMin, bbMax );

    // Reserve space and add the mesh BB center
    for(int i=0 ; i< 3 ; ++i)
    {
      coords[i].reserve(nOrigPts+2);
      coords[i].push_back( scaleAndOffset(bbMin[i], bbMax[i], 0.5));
      coords[i].push_back( scaleAndOffset(bbMin[i], bbMax[i], 0.5));
    }

    // Scale and add the random points from input parameter
    for(int j=0 ; j< nOrigPts ; ++j)
    {
      for(int i=0 ; i< 3 ; ++i)
      {
        coords[i].push_back(
          scaleAndOffset(bbMin[i], bbMax[i], points[j*3 + i]));
      }
    }
  }


  int nPoints = coords[0].size();
  for(int i=0 ; i< nPoints ; ++i)
  {
    // scale points within bounding box of mesh
    const double x = coords[0][i];
    const double y = coords[1][i];
    const double z = coords[2][i];

    const double phi = quest::signed_distance_evaluate( x,y,z );
    const int ins    = ( std::signbit( phi ) != 0) ? -1 : 1;

    SLIC_INFO(
      "Point (" << x << ", " << y << ", " << z << ") "
                << "is " << (ins ? "inside" : "outside") << " surface."
                <<" Distance is " << phi << ".");
  }

  quest::signed_distance_finalize();
}

void runQuestContainment(const std::string& fileName, const CoordsVec& points)
{
  const bool useDistance = false;
  const int unusedVar = -1;

  #ifdef AXOM_USE_MPI
  quest::initialize( MPI_COMM_WORLD, fileName, useDistance, 3, unusedVar,
                     unusedVar);
  #else
  quest::initialize( fileName, useDistance, 3, unusedVar, unusedVar );
  #endif

  outputMeshStats();

  CoordsVec coords[3];
  {
    int nOrigPts = points.size()/3;

    double bbMin[3], bbMax[3], cMass[3];
    quest::mesh_min_bounds(bbMin);
    quest::mesh_max_bounds(bbMax);
    quest::mesh_center_of_mass(cMass);

    // Reserve space and add the mesh BB center and center of mass
    for(int i=0 ; i< 3 ; ++i)
    {
      coords[i].reserve(nOrigPts+2);
      coords[i].push_back( scaleAndOffset(bbMin[i], bbMax[i], 0.5));
      coords[i].push_back( scaleAndOffset(cMass[i], cMass[i], 0.5));
    }

    // Scale and add the random points from input parameter
    for(int j=0 ; j< nOrigPts ; ++j)
    {
      for(int i=0 ; i< 3 ; ++i)
      {
        coords[i].push_back(
          scaleAndOffset(bbMin[i], bbMax[i], points[j*3 + i]));
      }
    }
  }


  int nPoints = coords[0].size();
  for(int i=0 ; i< nPoints ; ++i)
  {
    // scale points within bounding box of mesh
    const double x = coords[0][i];
    const double y = coords[1][i];
    const double z = coords[2][i];
    const int ins = quest::inside( x,y,z);

    SLIC_INFO(
      "Point (" << x << ", " << y << ", " << z << ") "
                << "is " << (ins ? "inside" : "outside") << " surface." );
  }

  quest::finalize();
}


/*!
 * \brief A simple example illustrating the use of the Quest C-Style interface.
 * \param [in] argc argument counter.
 * \param [in] argv argument vector.
 * \return rc return code.
 *
 * \note To run the example
 * \verbatim
 *
 *   [mpirun -np N] ./quest_interface_ex <stl_file>
 *
 * \endverbatim
 */
int main( int argc, char** argv )
{
#ifdef AXOM_USE_MPI
  // Initialize MPI
  MPI_Init( &argc, &argv );
#endif

  // Initialize Logger
  axom::slic::initialize();
  axom::slic::setLoggingMsgLevel( axom::slic::message::Info );

  axom::slic::LogStream* logStream;

#ifdef AXOM_USE_MPI
  std::string fmt = "[<RANK>][<LEVEL>]: <MESSAGE>\n";
  #ifdef AXOM_USE_LUMBERJACK
  const int RLIMIT = 8;
  logStream =
    new axom::slic::LumberjackStream(&std::cout,MPI_COMM_WORLD, RLIMIT, fmt);
  #else
  logStream =
    new axom::slic::SynchronizedStream(&std::cout,MPI_COMM_WORLD, fmt);
  #endif
#else
  std::string fmt = "[<LEVEL>]: <MESSAGE>\n";
  logStream = new axom::slic::GenericOutputStream(&std::cout, fmt);
#endif // AXOM_USE_MPI

  axom::slic::addStreamToAllMsgLevels( logStream );

  if(argc != 2)
  {
  #ifdef AXOM_USE_MPI
    SLIC_WARNING("Usage: [mpirun -np N] ./quest_interface_ex <stl_file>");
  #else
    SLIC_WARNING("Usage: ./quest_interface_ex <stl_file>");
  #endif
    axom::slic::finalize();
    axom::utilities::processAbort();
  }

  // Generate the query points
  std::string fileName = std::string(argv[1]);

  std::srand( time(NULL) );
  const int npoints = 10;
  const double lb = 0.;
  const double ub = 1.;
  CoordsVec ptVec;
  ptVec.reserve(3*npoints);
  for ( int ipnt=0 ; ipnt < npoints ; ++ipnt )
  {
    ptVec.push_back( getRandomDouble(lb,ub) );
    ptVec.push_back( getRandomDouble(lb,ub) );
    ptVec.push_back( getRandomDouble(lb,ub) );
  }

  SLIC_INFO("**Running distance queries using quest interface...");
  runQuestDistance(fileName, ptVec);
  SLIC_INFO("--");

  SLIC_INFO("**Running containment queries using quest interface...");
  runQuestContainment(fileName, ptVec);
  SLIC_INFO("--");


  axom::slic::finalize();

#ifdef AXOM_USE_MPI
  MPI_Finalize();
#endif

  return 0;
}

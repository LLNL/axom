/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

/*!
 *******************************************************************************
 * \file quest_interface.cpp
 *
 * \date Mar 16, 2016
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */

// Quest includes
#include "quest/quest.hpp"

// SLIC includes
#include "slic/slic.hpp"
#include "slic/SynchronizedStream.hpp"

// C/C++ includes
#include <cstdlib>   // for std::rand(), RAND_MAX
#include <ctime>     // for time
#include <iostream>  // for std::cout

// MPI includes
#include "mpi.h"

using namespace asctoolkit;

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
    quest::mesh_bounds_min(bbMin);
    quest::mesh_bounds_max(bbMax);
    SLIC_INFO("Mesh bounding box: "
            << "{ lower: (" << bbMin[0] <<"," << bbMin[1] <<","<< bbMin[2] << ")"
            << "; upper: (" << bbMax[0] <<"," << bbMax[1] <<","<< bbMax[2] << ")}"
            );

    // Obtain and log the mesh center of mass
    double cm[3];
    quest::mesh_center_of_mass(cm);
    SLIC_INFO("Mesh center of mass: "
            << "(" << cm[0] <<"," << cm[1] <<","<< cm[2] << ")"
            );
}

void runQuestDistance(const std::string& fileName, const CoordsVec& points)
{
    const bool useDistance = true;
    quest::initialize( MPI_COMM_WORLD, fileName, useDistance, 3, 25, 20 );

    outputMeshStats();

    CoordsVec coords[3];
    {
        int nOrigPts = points.size()/3;

        double bbMin[3], bbMax[3], cMass[3];
        quest::mesh_bounds_min(bbMin);
        quest::mesh_bounds_max(bbMax);
        quest::mesh_center_of_mass(cMass);

        // Reserve space and add the mesh BB center and center of mass
        for(int i=0; i< 3; ++i)
        {
            coords[i].reserve(nOrigPts+2);
            coords[i].push_back( scaleAndOffset(bbMin[i], bbMax[i], 0.5));
            coords[i].push_back( scaleAndOffset(cMass[i], cMass[i], 0.5));
        }

        // Scale and add the random points from input parameter
        for(int j=0; j< nOrigPts; ++j)
        {
            for(int i=0; i< 3; ++i)
            {
                coords[i].push_back( scaleAndOffset(bbMin[i], bbMax[i], points[j*3 + i]));
            }
        }
    }


    int nPoints = coords[0].size();
    for(int i=0; i< nPoints; ++i)
    {
        // scale points within bounding box of mesh
        const double x = coords[0][i];
        const double y = coords[1][i];
        const double z = coords[2][i];

        const double phi = quest::distance(x,y,z);
        const int ins = quest::inside( x,y,z );

        SLIC_INFO("Point (" << x << ", " << y << ", " << z << ") "
                  << "is " << (ins? "inside" : "outside") << " surface."
                  <<" Distance is " << phi << ".");
    }

    quest::finalize();
}

void runQuestContainment(const std::string& fileName, const CoordsVec& points)
{
    const bool useDistance = false;
    const int unusedVar = -1;
    quest::initialize( MPI_COMM_WORLD, fileName, useDistance, 3, unusedVar, unusedVar);

    outputMeshStats();

    CoordsVec coords[3];
    {
        int nOrigPts = points.size()/3;

        double bbMin[3], bbMax[3], cMass[3];
        quest::mesh_bounds_min(bbMin);
        quest::mesh_bounds_max(bbMax);
        quest::mesh_center_of_mass(cMass);

        // Reserve space and add the mesh BB center and center of mass
        for(int i=0; i< 3; ++i)
        {
            coords[i].reserve(nOrigPts+2);
            coords[i].push_back( scaleAndOffset(bbMin[i], bbMax[i], 0.5));
            coords[i].push_back( scaleAndOffset(cMass[i], cMass[i], 0.5));
        }

        // Scale and add the random points from input parameter
        for(int j=0; j< nOrigPts; ++j)
        {
            for(int i=0; i< 3; ++i)
            {
                coords[i].push_back( scaleAndOffset(bbMin[i], bbMax[i], points[j*3 + i]));
            }
        }
    }


    int nPoints = coords[0].size();
    for(int i=0; i< nPoints; ++i)
    {
        // scale points within bounding box of mesh
        const double x = coords[0][i];
        const double y = coords[1][i];
        const double z = coords[2][i];
        const int ins = quest::inside( x,y,z);

        SLIC_INFO("Point (" << x << ", " << y << ", " << z << ") "
                  << "is " << (ins? "inside" : "outside") << " surface."
                  );
    }

    quest::finalize();
}


/*!
 *******************************************************************************
 * \brief A simple example illustrating the use of the Quest C-Style interface.
 * \param [in] argc argument counter.
 * \param [in] argv argument vector.
 * \return rc return code.
 *
 * \note To run the example
 * \verbatim
 *
 *   [mpirun -np N] ./quest_interface <stl_file>
 *
 * \endverbatim
 *******************************************************************************
 */
int main( int argc, char**argv )
{
  // Initialize MPI
  MPI_Init( &argc, &argv );

  // Initialize Logger
  slic::initialize();
  slic::setLoggingMsgLevel( slic::message::Info );

  std::string fmt = "[<RANK>][<LEVEL>]: <MESSAGE>";
  slic::SynchronizedStream* sstream =
          new slic::SynchronizedStream(&std::cout,MPI_COMM_WORLD);
  sstream->setFormatString( fmt );
  slic::addStreamToAllMsgLevels( sstream );

  // Generate the query points
  std::string fileName = std::string(argv[1]);

  std::srand( time(NULL) );
  const int npoints = 10;
  const double lb = 0.;
  const double ub = 1.;
  CoordsVec ptVec;
  ptVec.reserve(3*npoints);
  for ( int ipnt=0; ipnt < npoints; ++ipnt )
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


  slic::finalize();

  MPI_Finalize();

  return 0;
}


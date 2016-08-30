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
    double bb0[3] = {quest::mesh_min_x(), quest::mesh_min_y(), quest::mesh_min_z()};
    double bb1[3] = {quest::mesh_max_x(), quest::mesh_max_y(), quest::mesh_max_z()};
    SLIC_INFO("Mesh bounding box: "
            << "{ lower: (" << bb0[0] <<"," << bb0[1] <<","<< bb0[2] << ")"
            << "; upper: (" << bb1[0] <<"," << bb1[1] <<","<< bb1[2] << ")}"
            );

    // Obtain and log the mesh center of mass
    double cm[3] =  { quest::mesh_center_of_mass_x()
                    , quest::mesh_center_of_mass_y()
                    , quest::mesh_center_of_mass_z()};
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
        coords[0].reserve(nOrigPts+2);
        coords[1].reserve(nOrigPts+2);
        coords[2].reserve(nOrigPts+2);

        CoordsVec::value_type xMin = quest::mesh_min_x();
        CoordsVec::value_type xMax = quest::mesh_max_x();
        CoordsVec::value_type yMin = quest::mesh_min_y();
        CoordsVec::value_type yMax = quest::mesh_max_y();
        CoordsVec::value_type zMin = quest::mesh_min_z();
        CoordsVec::value_type zMax = quest::mesh_max_z();

        // Add the BB center
        coords[0].push_back( scaleAndOffset(xMin, xMax, 0.5));
        coords[1].push_back( scaleAndOffset(yMin, yMax, 0.5));
        coords[2].push_back( scaleAndOffset(zMin, zMax, 0.5));

        // Add the mesh center of mass
        coords[0].push_back( quest::mesh_center_of_mass_x() );
        coords[1].push_back( quest::mesh_center_of_mass_y() );
        coords[2].push_back( quest::mesh_center_of_mass_z() );

        // Scale and add the random points
        for(int i=0; i< nOrigPts; ++i)
        {
            coords[0].push_back( scaleAndOffset(xMin, xMax, points[i*3]));
            coords[1].push_back( scaleAndOffset(yMin, yMax, points[i*3 +1]));
            coords[2].push_back( scaleAndOffset(zMin, zMax, points[i*3 +2]));
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
        coords[0].reserve(nOrigPts+2);
        coords[1].reserve(nOrigPts+2);
        coords[2].reserve(nOrigPts+2);

        CoordsVec::value_type xMin = quest::mesh_min_x();
        CoordsVec::value_type xMax = quest::mesh_max_x();
        CoordsVec::value_type yMin = quest::mesh_min_y();
        CoordsVec::value_type yMax = quest::mesh_max_y();
        CoordsVec::value_type zMin = quest::mesh_min_z();
        CoordsVec::value_type zMax = quest::mesh_max_z();

        // Add the BB center
        coords[0].push_back( scaleAndOffset(xMin, xMax, 0.5));
        coords[1].push_back( scaleAndOffset(yMin, yMax, 0.5));
        coords[2].push_back( scaleAndOffset(zMin, zMax, 0.5));

        // Add the mesh center of mass
        coords[0].push_back( quest::mesh_center_of_mass_x() );
        coords[1].push_back( quest::mesh_center_of_mass_y() );
        coords[2].push_back( quest::mesh_center_of_mass_z() );

        // Scale and add the random points
        for(int i=0; i< nOrigPts; ++i)
        {
            coords[0].push_back( scaleAndOffset(xMin, xMax, points[i*3]));
            coords[1].push_back( scaleAndOffset(yMin, yMax, points[i*3 +1]));
            coords[2].push_back( scaleAndOffset(zMin, zMax, points[i*3 +2]));
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


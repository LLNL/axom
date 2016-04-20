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


void runQuestDistance(const std::string& fileName, const CoordsVec& points)
{
    const bool useDistance = true;
    quest::initialize( MPI_COMM_WORLD, fileName, useDistance, 3, 25, 20 );

    int nPoints = points.size()/3;
    for(int i=0; i< nPoints; ++i)
    {
        const double x = points[i*3];
        const double y = points[i*3+1];
        const double z = points[i*3+2];

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

    int nPoints = points.size()/3;
    for(int i=0; i< nPoints; ++i)
    {
        const double x = points[i*3];
        const double y = points[i*3+1];
        const double z = points[i*3+2];
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

  // Initialize quest

  std::string fileName = std::string(argv[1]);

  std::srand( time(NULL) );
  const int npoints = 10;
  const double lb = -4.5;
  const double ub =  4.5;
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


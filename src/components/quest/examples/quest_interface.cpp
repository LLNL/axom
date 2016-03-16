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
//------------------------------------------------------------------------------
template < typename T >
T getRandomDouble( T low, T high )
{
  const T delta = high-low;
  const T c = static_cast< T >( std::rand() ) / static_cast< T >( RAND_MAX );
  return ( delta*c+low );
}

//------------------------------------------------------------------------------
int main( int argc, char**argv )
{
  // Initialize MPI
  MPI_Init( &argc, &argv );

  // Initialize Logger
  slic::initialize();
  slic::setLoggingMsgLevel( slic::message::Debug );

  std::string fmt = "[<RANK>][<LEVEL>]: <MESSAGE>";
  slic::SynchronizedStream* sstream =
          new slic::SynchronizedStream(&std::cout,MPI_COMM_WORLD);
  sstream->setFormatString( fmt );
  slic::addStreamToAllMsgLevels( sstream );

  // Initialize quest
  quest::initialize( MPI_COMM_WORLD, std::string(argv[1]), 3, 25, 10 );

  std::srand( time(NULL) );

  const int npoints = 10;
  for ( int ipnt=0; ipnt < npoints; ++ipnt ) {

      const double x = getRandomDouble( 0.0, 10.0 );
      const double y = getRandomDouble( 0.0, 10.0 );
      const double z = getRandomDouble( 0.0, 10.0 );

      const double phi = quest::distance( x, y, z );
      SLIC_INFO( "distance(" << x << ", " << y << ", " << z << ") = " << phi );

  } // END for all points

  // Finalize
  quest::finalize();

  slic::finalize();

  MPI_Finalize();
}


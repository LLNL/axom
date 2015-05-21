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
 * \file parallel_logging_example.cc
 *
 * \date May 7, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *
 *******************************************************************************
 */

// C/C++ includes
#include <cstdlib> // for rand()
#include <sstream> // for ostringstream

// Logging includes
#include "slic/Logger.hpp"
#include "slic/SynchronizedStream.hpp"

// MPI
#include <mpi.h>

using namespace asctoolkit;

#define N 20

slic::message::Level getRandomEvent( const int start, const int end )
{
  return( static_cast<slic::message::Level>(std::rand() % (end-start) + start));
}

//------------------------------------------------------------------------------
int main( int argc, char** argv )
{
  // STEP 0: initialize MPI & logging environment
  MPI_Init( &argc, &argv );

  int rank=-1;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  std::string format = std::string( "<MESSAGE>\n") +
                       std::string( "\t<TIMESTAMP>" ) +
                       std::string( "\tLEVEL=<LEVEL>\n") +
                       std::string( "\tFILE=<FILE>\n") +
                       std::string( "\tLINE=<LINE>\n");

  slic::Logger::initialize();


  slic::Logger::setLogLevel( slic::message::Debug );
  slic::Logger::addStream(
      new slic::SynchronizedStream( &std::cout, MPI_COMM_WORLD, format ) );

  // STEP 3: loop N times and generate a random logging event
  for ( int i=0; i < N; ++i ) {

    std::ostringstream oss;
    oss << "[ " << rank << "]: message " << i << "/" << N-1;

    slic::Logger::log( getRandomEvent(0,slic::message::Num_Levels),
                         oss.str(),
                         __FILE__,
                         __LINE__
                         );

    // Flush every 5 cycles
    if ( (i % 5)==0 ) {

      slic::Logger::flushStreams();

    } // END if

  }

  // STEP 4: shutdown logging environment
  slic::Logger::finalize();

  // STEP 5: Finalize MPI
  MPI_Finalize();

  return 0;
}

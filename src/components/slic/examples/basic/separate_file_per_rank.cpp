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
#include <fstream> // for ofstream
#include <sstream> // for ostringstream

// Logging includes
#include "slic/slic.hpp"
#include "slic/GenericOutputStream.hpp"

// MPI
#include <mpi.h>

using namespace axom;

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
  slic::initialize();

  int rank = -1;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank);

  std::ostringstream oss;
  oss << "logfile_" << rank << ".dat";

  std::ofstream ofs;
  ofs.open( oss.str().c_str() );

  slic::setLoggingMsgLevel( slic::message::Debug );
  slic::disableAbortOnError();
  slic::addStreamToAllMsgLevels( new slic::GenericOutputStream(&ofs) );


  // STEP 3: loop N times and generate a random logging event
  for ( int i=0; i < N; ++i ) {

    slic::logMessage(
        getRandomEvent(0,slic::message::Num_Levels),
        "a random message", __FILE__,  __LINE__  );

    // Flush every 5 cycles
    if ( (i % 5)==0 ) {

      slic::flushStreams();

    } // END if

  }

  ofs.close();

  // STEP 4: shutdown logging environment
  slic::finalize();


  // STEP 5: Finalize MPI
  MPI_Finalize();

  return 0;
}

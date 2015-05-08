/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
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

// Logging includes
#include "logapi/Logger.hpp"
#include "logapi/SeperateFilePerRankStream.hpp"

// MPI
#include <mpi.h>

using namespace asctoolkit;

#define N 20

logapi::MessageType getRandomEvent( const int start, const int end )
{
  return( static_cast<logapi::MessageType>(std::rand() % (end-start) + start));
}

//------------------------------------------------------------------------------
int main( int argc, char** argv )
{
  // STEP 0: initialize MPI & logging environment
  MPI_Init( &argc, &argv );
  logapi::Logger::initialize();

  // STEP 1: Setup the log stream
  logapi::SeperateFilePerRankStream* scls =
      new logapi::SeperateFilePerRankStream("logfile",MPI_COMM_WORLD);

  // STEP 2: enable logging of all messages
  for ( int mtype=0; mtype < logapi::Num_Msg_Types; ++mtype ) {

    logapi::Logger::getInstance()->enable(
        static_cast<logapi::MessageType>(mtype) );
    logapi::Logger::getInstance()->setLogStream(
        static_cast<logapi::MessageType>(mtype), scls);

  } // END for all message types

  // STEP 3: loop N times and generate a random logging event
  for ( int i=0; i < N; ++i ) {

    logapi::Logger::getInstance()->logMessage(
            getRandomEvent(0,logapi::Num_Msg_Types),
            "a random message", __FILE__,  __LINE__  );

    // Flush every 5 cycles
    if ( (i % 5)==0 ) {

      logapi::Logger::getInstance()->flushAllStreams();

    } // END if

  }

  // STEP 4: shutdown logging environment
  logapi::Logger::finalize();

  delete scls;

  // STEP 5: Finalize MPI
  MPI_Finalize();

  return 0;
}

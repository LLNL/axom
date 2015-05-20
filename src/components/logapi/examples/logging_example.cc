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
 * \file logging_example.cc
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
#include "logapi/Console.hpp"

using namespace asctoolkit;

#define N 10

logapi::message::Level getRandomEvent( const int start, const int end )
{
  return( static_cast<logapi::message::Level>(std::rand() % (end-start) + start));
}

//------------------------------------------------------------------------------
int main( int argc, char** argv )
{
  static_cast<void>(argc);
  static_cast<void>(argv);

  // STEP 0: initialize logging environment
  logapi::Logger::initialize();

  logapi::Logger::setLogLevel( logapi::message::Debug );
  logapi::Logger::addStream( new logapi::Console() );


  // STEP 3: loop N times and generate a random logging event
  for ( int i=0; i < N; ++i ) {

    logapi::Logger::log( getRandomEvent(0,logapi::message::Num_Levels),
            "a random message", __FILE__,  __LINE__  );

  }

  // STEP 4: shutdown logging environment
  logapi::Logger::finalize();

  return 0;
}



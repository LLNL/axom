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
#include "slic/slic.hpp"
#include "slic/GenericOutputStream.hpp"

using namespace axom;

#define N 10

slic::message::Level getRandomEvent( const int start, const int end )
{
  return( static_cast<slic::message::Level>(std::rand() % (end-start) + start));
}

//------------------------------------------------------------------------------
int main( int argc, char** argv )
{
  static_cast<void>(argc);
  static_cast<void>(argv);

  // STEP 0: initialize logging environment
  slic::initialize();
  slic::disableAbortOnError();

  std::string format =  std::string( "***********************************\n" )+
                        std::string( "* <TIMESTAMP>\n\n" ) +
                        std::string( "* LEVEL=<LEVEL>\n" ) +
                        std::string( "* MESSAGE=<MESSAGE>\n" ) +
                        std::string( "* FILE=<FILE>\n" ) +
                        std::string( "* LINE=<LINE>\n" ) +
                        std::string( "***********************************\n" );

  slic::setLoggingMsgLevel( slic::message::Debug );
  slic::addStreamToAllMsgLevels(
      new slic::GenericOutputStream( &std::cout, format ) );


  // STEP 1: loop N times and generate a random logging event
  for ( int i=0; i < N; ++i ) {

    slic::logMessage( getRandomEvent(0,slic::message::Num_Levels),
            "a random message", __FILE__,  __LINE__  );

  }

  // STEP 2: shutdown logging environment
  slic::finalize();

  return 0;
}



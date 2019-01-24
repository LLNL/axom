/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
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

// C/C++ includes
#include <cstdlib> // for rand()

// Logging includes
#include "axom/slic/interface/slic.hpp"
#include "axom/slic/streams/GenericOutputStream.hpp"

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
  for ( int i=0 ; i < N ; ++i )
  {

    slic::logMessage( getRandomEvent(0,slic::message::Num_Levels),
                      "a random message", __FILE__,  __LINE__  );

  }

  // STEP 2: shutdown logging environment
  slic::finalize();

  return 0;
}

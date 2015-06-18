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
 * \file cpplogger_example.cc
 *
 * \date Jun 18, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *
 *******************************************************************************
 */

// C/C++ includes
#include <fstream>  // for file stream
#include <cstdlib>  // for rand()

// Logging includes
#include "slic/Logger.hpp"
#include "slic/LogStream.hpp"
#include "slic/GenericOutputStream.hpp"

using namespace asctoolkit;
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

  //----------------------------------------------------------------------------
  // STEP 0: Initialize logger
  //----------------------------------------------------------------------------
  slic::Logger::initialize();
  slic::Logger::getInstance()->setLoggingLevel(slic::message::Debug);

  //----------------------------------------------------------------------------
  // STEP 1: Create log streams
  //----------------------------------------------------------------------------

  // setup log stream for FATAL, ERROR and WARNING messages
  std::ofstream hspStream;
  hspStream.open( "HSP.dat" );

  std::string hsp_format =
      std::string( "***********************************\n" )+
      std::string( "* <TIMESTAMP>\n\n" ) +
      std::string( "* LEVEL=<LEVEL>\n" ) +
      std::string( "* MESSAGE=<MESSAGE>\n" ) +
      std::string( "* FILE=<FILE>\n" ) +
      std::string( "* LINE=<LINE>\n" ) +
      std::string( "***********************************\n" );

  slic::LogStream* hspLogStream =
      new slic::GenericOutputStream(&hspStream, hsp_format);

  // setup log stream for ALL messages, including FATAL, ERROR and WARNING
  std::string console_format = std::string("[<LEVEL>]: <MESSAGE>\n");
  slic::LogStream* console =
      new slic::GenericOutputStream( &std::cerr, console_format );

  //----------------------------------------------------------------------------
  // STEP 2: add streams to logger
  //----------------------------------------------------------------------------
  slic::Logger::getInstance()->addStreamToLevel(
      hspLogStream,slic::message::Fatal);
  slic::Logger::getInstance()->addStreamToLevel(
      hspLogStream,slic::message::Error);
  slic::Logger::getInstance()->addStreamToLevel(
      hspLogStream,slic::message::Warning);

  slic::Logger::getInstance()->addStreamToAllLevels( console );

  //----------------------------------------------------------------------------
  // STEP 3: Loop N times and generate random logging events
  //----------------------------------------------------------------------------
  for ( int i=0; i < N; ++i ) {

    slic::Logger::getInstance()->logMessage(
        getRandomEvent(0,slic::message::Num_Levels),
        "a random message", __FILE__, __LINE__  );

  }

  //----------------------------------------------------------------------------
  // STEP 4: Loop N times and generate random logging events
  //----------------------------------------------------------------------------

  // finalize logging & close HSP file
  slic::Logger::finalize();
  hspStream.close();

  return 0;
}


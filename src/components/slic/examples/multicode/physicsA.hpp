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
 * \file physicsA.hpp
 *
 * \date Jul 28, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *
 *******************************************************************************
 */

#ifndef PHYSICSA_HPP_
#define PHYSICSA_HPP_

// SLIC includes
#include "slic/slic.hpp"
#include "slic/GenericOutputStream.hpp"

using namespace asctoolkit;

namespace physicsA {

std::ofstream physicsA_log;

int getRandInt( const int start, const int end )
{
  return( std::rand() % (end-start) + start );
}

slic::message::Level getRandomLevel()
{
  return( static_cast< slic::message::Level >(
                getRandInt(0,slic::message::Num_Levels)) );
}
//------------------------------------------------------------------------------
void init()
{
  std::string physicsA_format =
      std::string( "====\n" ) +
      std::string( "<TIMESTAMP>" ) +
      std::string( "====\n" ) +
      std::string( "[<LEVEL>]: <MESSAGE>\n" ) +
      std::string( "\t FILE:<FILE>\n" ) +
      std::string( "\t LINE:<LINE>\n" );

  physicsA_log.open( "physicsA.log" );
  slic::LogStream* ls = new slic::GenericOutputStream(&physicsA_log, physicsA_format);

  slic::createLogger( "physicsA", slic::inherit::errors_and_warnings );
  slic::activateLogger( "physicsA" );
  slic::setLoggingLevel( slic::message::Debug );
  slic::addStreamToAllLevels( ls );

  slic::activateLogger("root");
}

//------------------------------------------------------------------------------
void timestep(int step, int n)
{
  slic::activateLogger( "physicsA" );

  std::ostringstream oss;
  oss << "n=" << n << " physicsA cycles";
  slic::logMessage( slic::message::Info, oss.str(), __FILE__, __LINE__ );

  for ( int i=0; i < n; ++i ) {

    slic::message::Level random = getRandomLevel();

    oss.str("");
    oss << "cycle=" << step << " subcycle=" << i << " a random message!";
    slic::logMessage(random,oss.str(),__FILE__,__LINE__);
  }

  slic::activateLogger( "root" );
}

//------------------------------------------------------------------------------
void finalize()
{
  physicsA_log.close();
}

} /* namespace physicsA */



#endif /* PHYSICSA_HPP_ */

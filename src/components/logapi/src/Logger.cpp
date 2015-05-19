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
 * \file Logger.cpp
 *
 * \date May 7, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *
 *******************************************************************************
 */

#include "Logger.hpp"

#include "LogStream.hpp"

// C/C++ includes
#include <cassert> // for assert()
#include <cstddef> // for NULL

namespace asctoolkit {

namespace logapi {

Logger* Logger::s_Logger = NULL;


//------------------------------------------------------------------------------
Logger::Logger()
{
  // by default, all message streams are disabled
  for ( int i=0; i < message::Num_Levels; ++i ) {

    m_StreamState[ i ] = false;
    m_Streams[ i ]     = NULL;
  }

}

//------------------------------------------------------------------------------
Logger::~Logger()
{

}

//------------------------------------------------------------------------------
void Logger::enableStreamsBelow( message::Level level )
{
  assert("pre: invalid message type" &&
          (level >= 0) && (level < message::Num_Levels));

  for ( int i=0; i <= message::Num_Levels; ++i ) {
    m_StreamState[ i ] = true;
  }

}

//------------------------------------------------------------------------------
void Logger::setStreamsBelow( message::Level level, LogStream* ls)
{
  assert("pre: invalid message type" &&
         (level >= 0) && (level < message::Num_Levels));
  assert("pre: supplied log stream is NULL!" && (ls != NULL) );

  for ( int i=0; i <= message::Num_Levels; ++i ) {
    m_Streams[ i ] = ls;
  }

}

//------------------------------------------------------------------------------
void Logger::enable( message::Level level)
{
  assert("pre: invalid message type" &&
          (level >= 0) && (level < message::Num_Levels));
  m_StreamState[ level ] = true;
}

//------------------------------------------------------------------------------
void Logger::disable( message::Level level)
{
  assert("pre: invalid message type" &&
          (level >= 0) && (level < message::Num_Levels));
  m_StreamState[ level ] = false;
}

//------------------------------------------------------------------------------
void Logger::setLogStream( message::Level level, LogStream* ls )
{
  assert("pre: supplied log stream is NULL!" && ls != NULL );
  assert("pre: invalid message type" &&
          (level >= 0) && (level < message::Num_Levels));

  m_Streams[ level ] = ls;
}

//------------------------------------------------------------------------------
void Logger::logMessage( message::Level level,
                         const std::string& message,
                         const std::string& fileName,
                         int line )
{
  if ( m_StreamState[ level ]==false  ) {

    /* short-circuit */
    return;

  } // END if

  assert( "pre: no stream set for type!" && m_Streams[ level ] != NULL );

  m_Streams[ level ]->append( level,message,fileName,line );
}

//------------------------------------------------------------------------------
void Logger::flushAllStreams()
{
  for ( int i=0; i < message::Num_Levels; ++i ) {

    if ( m_StreamState[ i ] ) {

      m_Streams[ i ]->flush( );

    } // END if

  }

}

//------------------------------------------------------------------------------
void Logger::initialize()
{
  if( s_Logger == NULL ) {
    s_Logger = new Logger();
  }
}

//------------------------------------------------------------------------------
void Logger::finalize()
{
  delete s_Logger;
  s_Logger = NULL;
}

//------------------------------------------------------------------------------
Logger* Logger::getInstance()
{
  assert( "pre: uninitialized logger!" && (s_Logger != NULL) );
  return s_Logger;
}

} /* namespace logapi */

} /* namespace asctoolkit */

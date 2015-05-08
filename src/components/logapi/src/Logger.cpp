/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
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


const std::string MsgTypeName[ Num_Msg_Types ] = {
    "FATAL",
    "ERROR",
    "WARNING",
    "INFO",
    "DEBUG"
};


//------------------------------------------------------------------------------
Logger::Logger()
{
  // by default, all message streams are disabled
  for ( int i=0; i < Num_Msg_Types; ++i ) {

    m_StreamState[ i ] = false;
    m_Streams[ i ]     = NULL;
  }

}

//------------------------------------------------------------------------------
Logger::~Logger()
{

}

//------------------------------------------------------------------------------
void Logger::enable( int type)
{
  assert("pre: invalid message type" && (type >= 0) && (type < Num_Msg_Types));
  m_StreamState[ type ] = true;
}

//------------------------------------------------------------------------------
void Logger::disable( int type)
{
  assert("pre: invalid message type" && (type >= 0) && (type < Num_Msg_Types));
  m_StreamState[ type ] = false;
}

//------------------------------------------------------------------------------
void Logger::setLogStream( int type, LogStream* ls )
{
  assert("pre: supplied log stream is NULL!" && ls != NULL );
  assert("pre: invalid message type" && (type >= 0) && (type < Num_Msg_Types));

  m_Streams[ type ] = ls;
}

//------------------------------------------------------------------------------
void Logger::logMessage( int type,
                         const std::string& message,
                         const std::string& fileName,
                         int line )
{
  if ( m_StreamState[ type ]==false  ) {

    /* short-circuit */
    return;

  } // END if

  assert( "pre: no stream set for type!" && m_Streams[ type ] != NULL );

  m_Streams[ type ]->append( type,MsgTypeName[ type ],message,fileName,line );
}

//------------------------------------------------------------------------------
void Logger::flushAllStreams()
{
  for ( int i=0; i < Num_Msg_Types; ++i ) {

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

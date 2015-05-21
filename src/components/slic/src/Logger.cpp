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

#include "common/Types.hpp"

// C/C++ includes
#include <iostream> // for std::cout, std::cerr

namespace asctoolkit {

namespace slic {

Logger* Logger::s_Logger = static_cast< Logger* >( ATK_NULLPTR );


//------------------------------------------------------------------------------
Logger::Logger()
{
  // by default, all message streams are disabled
  for ( int i=0; i < message::Num_Levels; ++i ) {

    m_isEnabled[ i ] = false;

  }

}

//------------------------------------------------------------------------------
Logger::~Logger()
{
  unsigned nstreams = m_logStreams.size( );
  for ( unsigned istream=0; istream < nstreams; ++istream ) {

    delete m_logStreams[ istream ];
    m_logStreams[ istream ] = static_cast< LogStream* >( ATK_NULLPTR );

  } // END for all streams

  m_logStreams.clear() ;

}

//------------------------------------------------------------------------------
void Logger::setLoggingLevel( message::Level level )
{
  for ( int i=message::Fatal; i < level; ++i ) {
    m_isEnabled[ i ] = true;
  }

}

//------------------------------------------------------------------------------
void Logger::addLogStream( LogStream* ls )
{
  if ( ls == ATK_NULLPTR ) {

    std::cerr << "WARNING: supplied log stream is NULL!\n";
    return;

  }

  m_logStreams.push_back( ls );

}

//------------------------------------------------------------------------------
void Logger::logMessage( message::Level level,
                         const std::string& message )
{
  this->logMessage(
      level,message,MSG_IGNORE_TAG,MSG_IGNORE_FILE,MSG_IGNORE_LINE );
}

//------------------------------------------------------------------------------
void Logger::logMessage( message::Level level,
                         const std::string& message,
                         const std::string& tagName )
{
  this->logMessage( level,message,tagName,MSG_IGNORE_FILE,MSG_IGNORE_LINE );
}

//------------------------------------------------------------------------------
void Logger::logMessage( message::Level level,
                         const std::string& message,
                         const std::string& fileName,
                         int line )
{
  this->logMessage( level, message, MSG_IGNORE_TAG, fileName, line );
}

//------------------------------------------------------------------------------
void Logger::logMessage( message::Level level,
                         const std::string& message,
                         const std::string& tagName,
                         const std::string& fileName,
                         int line )
{
  if ( m_isEnabled[ level ]==false  ) {

    /* short-circuit */
    return;

  } // END if

  unsigned nstreams = m_logStreams.size();
  for ( unsigned istream=0; istream < nstreams; ++istream ) {

    m_logStreams[ istream ]->append( level, message, tagName, fileName, line );

  } // END for all streams

}

//------------------------------------------------------------------------------
void Logger::flushAllStreams()
{
  unsigned nstreams = m_logStreams.size();
  for ( unsigned istream=0; istream < nstreams; ++istream ) {

    m_logStreams[ istream ]->flush( );

  } // END for all streams

}

//------------------------------------------------------------------------------
void Logger::initialize()
{
  if( s_Logger == ATK_NULLPTR ) {
    s_Logger = new Logger();
  }
}

//------------------------------------------------------------------------------
void Logger::setLogLevel( message::Level level )
{
  if ( s_Logger == ATK_NULLPTR ) {

    std::cerr << "ERROR: `" << __FUNCTION__
              << "` is called with a NULL logger!";

  }
  s_Logger->setLoggingLevel( level );
}

//------------------------------------------------------------------------------
void Logger::addStream( LogStream* ls )
{
  if ( s_Logger == ATK_NULLPTR ) {

    std::cerr << "ERROR: `" << __FUNCTION__
              << "` is called with a NULL logger!";

  }
  s_Logger->addLogStream( ls );
}

//------------------------------------------------------------------------------
void Logger::log( message::Level level,
                  const std::string& message )
{
  if ( s_Logger == ATK_NULLPTR ) {

    std::cerr << "ERROR: `" << __FUNCTION__
              << "` is called with a NULL logger!";

  }
  s_Logger->logMessage( level,message );
}

//------------------------------------------------------------------------------
void Logger::log( message::Level level,
                  const std::string& message,
                  const std::string& tag )
{
  if ( s_Logger == ATK_NULLPTR ) {

    std::cerr << "ERROR: `" << __FUNCTION__
              << "` is called with a NULL logger!";

  }
  s_Logger->logMessage( level,message,tag );
}

//------------------------------------------------------------------------------
void Logger::log( message::Level level,
                  const std::string& message,
                  const std::string& fileName,
                  int line )
{
  if ( s_Logger == ATK_NULLPTR ) {

    std::cerr << "ERROR: `"
              << __FUNCTION__ << "` is called with a NULL logger!";

  }
  s_Logger->logMessage( level,message,fileName,line );
}

//------------------------------------------------------------------------------
void Logger::log( message::Level level,
                  const std::string& message,
                  const std::string& tag,
                  const std::string& fileName,
                  int line )
{
  if ( s_Logger == ATK_NULLPTR ) {

    std::cerr << "ERROR: `" << __FUNCTION__
              << "` is called with a NULL logger!";

  }
  s_Logger->logMessage( level,message,tag,fileName,line );
}

//------------------------------------------------------------------------------
void Logger::flushStreams()
{
  if ( s_Logger == ATK_NULLPTR ) {

    std::cerr << "ERROR: `" << __FUNCTION__
              << "` is called with a NULL logger!";

  }
  s_Logger->flushAllStreams();
}

//------------------------------------------------------------------------------
void Logger::finalize()
{
  s_Logger->flushAllStreams();

  delete s_Logger;
  s_Logger = static_cast< Logger* >( ATK_NULLPTR );
}

//------------------------------------------------------------------------------
Logger* Logger::getInstance()
{
  if ( s_Logger == ATK_NULLPTR ) {

    std::cerr << "ERROR: `" << __FUNCTION__
              << "` is called with a NULL logger!";

  }
  return s_Logger;
}

} /* namespace slic */

} /* namespace asctoolkit */

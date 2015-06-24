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

#include "common/CommonTypes.hpp"

// C/C++ includes
#include <iostream> // for std::cout, std::cerr

namespace asctoolkit {

namespace slic {

Logger* Logger::s_Logger = ATK_NULLPTR;


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
  std::map< LogStream*, LogStream* >::iterator it =
      m_streamObjectsManager.begin();
  for ( ; it != m_streamObjectsManager.end(); ++it ) {
    delete it->second;
  } // END for all logStreams

  for ( int level=message::Fatal; level < message::Num_Levels; ++level ) {

    m_logStreams[ level ].clear();

  } // END for all levels

}

//------------------------------------------------------------------------------
void Logger::setLoggingLevel( message::Level level )
{
  for ( int i=message::Fatal; i < level; ++i ) {
    m_isEnabled[ i ] = true;
  }

}

//------------------------------------------------------------------------------
void Logger::addStreamToLevel( LogStream* ls, message::Level level )
{
  if ( ls == ATK_NULLPTR ) {

    std::cerr << "WARNING: supplied log stream is NULL!\n";
    return;

  }

  m_logStreams[ level ].push_back( ls );
  m_streamObjectsManager[ ls ] = ls;

}

//------------------------------------------------------------------------------
void Logger::addStreamToAllLevels( LogStream* ls )
{
  if ( ls == ATK_NULLPTR ) {

    std::cerr << "WARNING: supplied log stream is NULL!\n";
    return;

  }

  m_streamObjectsManager[ ls ] = ls;
  for ( int level=message::Fatal; level < message::Num_Levels; ++level ) {

    m_logStreams[ level ].push_back( ls );

  } // END for all levels

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

  unsigned nstreams = m_logStreams[ level ].size();
  for ( unsigned istream=0; istream < nstreams; ++istream ) {

    m_logStreams[ level ][ istream ]->append(
                              level,message,tagName,fileName,line);

  } // END for all streams

}

//------------------------------------------------------------------------------
void Logger::flushStreams()
{
  for ( int level=message::Fatal; level < message::Num_Levels; ++level ) {

    unsigned nstreams = m_logStreams[ level ].size();
    for ( unsigned istream=0; istream < nstreams; ++istream ) {

      m_logStreams[ level ][ istream ]->flush( );

    } // END for all streams

  } // END for all levels

}

//------------------------------------------------------------------------------
void Logger::initialize()
{
  if( s_Logger == ATK_NULLPTR ) {
    s_Logger = new Logger();
  }
}

//------------------------------------------------------------------------------
void Logger::finalize()
{
  s_Logger->flushStreams();

  delete s_Logger;
  s_Logger = ATK_NULLPTR;
}

//------------------------------------------------------------------------------
Logger* Logger::getInstance()
{
  return s_Logger;
}

} /* namespace slic */

} /* namespace asctoolkit */

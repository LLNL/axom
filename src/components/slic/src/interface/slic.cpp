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
 * \file slic.cpp
 *
 * \date Jun 18, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *
 *******************************************************************************
 */

#include "slic/slic.hpp"

#include "common/Utilities.hpp"   // for utilities::processAbort()

#include <cstdlib>    // for free
#include <sstream>    // for std::ostringstream
#include <execinfo.h> // for backtrace()

namespace asctoolkit {

namespace slic {


//------------------------------------------------------------------------------
// Initialize static variables for controlling runtime behavior of asserts and
// error macros.
//------------------------------------------------------------------------------
bool RuntimeAbortBehavior::willAbortOnAssert = true;
bool RuntimeAbortBehavior::willAbortOnError = true;

void initialize()
{
  Logger::initialize();
}

//------------------------------------------------------------------------------
bool isInitialized()
{
  return ( Logger::getActiveLogger() != ATK_NULLPTR );
}

//------------------------------------------------------------------------------
void createLogger( const std::string& name, char imask )
{
  Logger::createLogger( name, imask );
}

//------------------------------------------------------------------------------
void activateLogger( const std::string& name )
{
  Logger::activateLogger( name );
}

//------------------------------------------------------------------------------
std::string getActiveLoggerName()
{
  return ( Logger::getActiveLoggerName() );
}

//------------------------------------------------------------------------------
void setLoggingMsgLevel( message::Level level )
{
  Logger::getActiveLogger()->setLoggingMsgLevel( level );
}

//------------------------------------------------------------------------------
void setAbortOnAssert( bool willAbort )
{
  RuntimeAbortBehavior::willAbortOnAssert = willAbort;
}

//------------------------------------------------------------------------------
bool getAbortOnAssert()
{
  return RuntimeAbortBehavior::willAbortOnAssert;
}


//------------------------------------------------------------------------------
void setAbortOnError( bool willAbort )
{
  RuntimeAbortBehavior::willAbortOnError = willAbort;
}

//------------------------------------------------------------------------------
bool getAbortOnError()
{
  return RuntimeAbortBehavior::willAbortOnError;
}

//------------------------------------------------------------------------------
void addStreamToMsgLevel( LogStream* ls, message::Level level )
{
  Logger::getActiveLogger()->addStreamToMsgLevel( ls, level );
}

//------------------------------------------------------------------------------
void addStreamToAllMsgLevels( LogStream* ls )
{
  Logger::getActiveLogger()->addStreamToAllMsgLevels( ls );
}

//------------------------------------------------------------------------------
void logMessage( message::Level level, const std::string& message,
                 bool filter_duplicates )
{
  if ( !isInitialized() ) {
    return;
  }
  Logger::getActiveLogger()->logMessage( level, message, filter_duplicates );
}

//------------------------------------------------------------------------------
void logMessage( message::Level level,
                 const std::string& message,
                 const std::string& tag,
                 bool filter_duplicates )
{
  if ( !isInitialized() ) {
     return;
  }
  Logger::getActiveLogger()->logMessage( level, message, tag,
                                         filter_duplicates );
}

//------------------------------------------------------------------------------
void logMessage( message::Level level,
                const std::string& message,
                const std::string& fileName,
                int line,
                bool filter_duplicates )
{
  if ( !isInitialized() ) {
     return;
  }
  Logger::getActiveLogger()->logMessage( level, message, fileName, line,
                                         filter_duplicates );
}

//------------------------------------------------------------------------------
void logMessage( message::Level level,
                 const std::string& message,
                 const std::string& tag,
                 const std::string& fileName,
                 int line,
                 bool filter_duplicates )
{
  if ( !isInitialized() ) {
     return;
  }
  Logger::getActiveLogger()->logMessage( level, message, tag, fileName, line,
                                         filter_duplicates );
}

//------------------------------------------------------------------------------
void logErrorMessage( const std::string& message,
                      const std::string& fileName,
                      int line,
                      bool shouldAbort )
{
  std::ostringstream oss;
  oss << message << slic::stacktrace();

  slic::logMessage( message::Fatal, oss.str(), fileName, line );
  if ( shouldAbort ) {
    asctoolkit::utilities::processAbort();
  }
}

//------------------------------------------------------------------------------
void logWarningMessage( const std::string& message,
                        const std::string& fileName,
                        int line )
{
  slic::logMessage( message::Warning, message, fileName, line );
}

//------------------------------------------------------------------------------
void flushStreams()
{
  Logger::getActiveLogger()->flushStreams();
}

//------------------------------------------------------------------------------
void finalize()
{
  Logger::finalize();
}

//------------------------------------------------------------------------------
std::string stacktrace( )
{
  void* array[10];
  const int size = backtrace( array, 10 );
  char** strings = backtrace_symbols( array, size );

  std::ostringstream oss;
  oss << "\n** StackTrace of " << size << " frames **\n";
  for ( int i=0; i < size; ++i ) {
     oss << strings[ i ] << std::endl;
  }
  oss << "=====\n\n";


  free( strings );

  return ( oss.str() );
}

} /* namespace slic */

} /* namespace asctoolkit */


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
void setLoggingLevel( message::Level level )
{
  Logger::getActiveLogger()->setLoggingLevel( level );
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
void addStreamToLevel( LogStream* ls, message::Level level )
{
  Logger::getActiveLogger()->addStreamToLevel( ls, level );
}

//------------------------------------------------------------------------------
void addStreamToAllLevels( LogStream* ls )
{
  Logger::getActiveLogger()->addStreamToAllLevels( ls );
}

//------------------------------------------------------------------------------
void logMessage( message::Level level, const std::string& message )
{
  if ( !isInitialized() ) {
    return;
  }
  Logger::getActiveLogger()->logMessage( level, message );
}

//------------------------------------------------------------------------------
void logMessage( message::Level level,
                 const std::string& message,
                 const std::string& tag )
{
  if ( !isInitialized() ) {
     return;
  }
  Logger::getActiveLogger()->logMessage( level, message, tag );
}

//------------------------------------------------------------------------------
void logMessage( message::Level level,
                const std::string& message,
                const std::string& fileName,
                int line )
{
  if ( !isInitialized() ) {
     return;
  }
  Logger::getActiveLogger()->logMessage( level, message, fileName, line );
}

//------------------------------------------------------------------------------
void logMessage( message::Level level,
                 const std::string& message,
                 const std::string& tag,
                 const std::string& fileName,
                 int line )
{
  if ( !isInitialized() ) {
     return;
  }
  Logger::getActiveLogger()->logMessage( level, message, tag, fileName, line );
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

} /* namespace slic */

} /* namespace asctoolkit */


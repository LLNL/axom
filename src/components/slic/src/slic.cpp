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

void initialize()
{
  Logger::initialize();
}

//------------------------------------------------------------------------------
bool isInitialized()
{
  return ( Logger::getInstance() != ATK_NULLPTR );
}

//------------------------------------------------------------------------------
void setLoggingLevel( message::Level level )
{
  Logger::getInstance()->setLoggingLevel( level );
}

//------------------------------------------------------------------------------
void addStreamToLevel( LogStream* ls, message::Level level )
{
  Logger::getInstance()->addStreamToLevel( ls, level );
}

//------------------------------------------------------------------------------
void addStreamToAllLevels( LogStream* ls )
{
  Logger::getInstance()->addStreamToAllLevels( ls );
}

//------------------------------------------------------------------------------
void logMessage( message::Level level, const std::string& message )
{
  if ( !isInitialized() ) {
    return;
  }
  Logger::getInstance()->logMessage( level, message );
}

//------------------------------------------------------------------------------
void logMessage( message::Level level,
                 const std::string& message,
                 const std::string& tag )
{
  if ( !isInitialized() ) {
     return;
  }
  Logger::getInstance()->logMessage( level, message, tag );
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
  Logger::getInstance()->logMessage( level, message, fileName, line );
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
  Logger::getInstance()->logMessage( level, message, tag, fileName, line );
}

//------------------------------------------------------------------------------
void flushStreams()
{
  Logger::getInstance()->flushStreams();
}

//------------------------------------------------------------------------------
void finalize()
{
  Logger::finalize();
}

} /* namespace slic */

} /* namespace asctoolkit */


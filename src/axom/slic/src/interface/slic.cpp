/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
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

#include "slic/slic.hpp"

#include <cstdlib>    // for free
#include <sstream>    // for std::ostringstream

#ifdef WIN32
  #define NOMINMAX
  #include <Windows.h>
  #include <WinBase.h>
  #include <DbgHelp.h>
#else
  #include <execinfo.h> // for backtrace()
#endif


namespace axom
{
namespace slic
{

//------------------------------------------------------------------------------
// Initialize static variables for controlling runtime behavior of asserts and
// error macros.
//------------------------------------------------------------------------------
bool debug::checksAreErrors = false;

void initialize()
{
  Logger::initialize();
}

//------------------------------------------------------------------------------
bool isInitialized()
{
  return ( Logger::getActiveLogger() != AXOM_NULLPTR );
}

//------------------------------------------------------------------------------
void createLogger( const std::string& name, char imask )
{
  if ( !isInitialized() )
  {
    std::cerr << "[ERROR]: slic::initialize() must be called first "
              << "before making any other calls to SLIC.";
    return;
  }
  Logger::createLogger( name, imask );
}

//------------------------------------------------------------------------------
bool activateLogger( const std::string& name )
{
  if ( !isInitialized() )
  {
    std::cerr << "[ERROR]: slic::initialize() must be called first "
              << "before making any other calls to SLIC.";
    return false;
  }
  return Logger::activateLogger( name );
}

//------------------------------------------------------------------------------
std::string getActiveLoggerName()
{
  if ( !isInitialized() )
  {
    std::cerr << "[ERROR]: slic::initialize() must be called first "
              << "before making any other calls to SLIC.";
    return "";
  }
  return ( Logger::getActiveLoggerName() );
}

//------------------------------------------------------------------------------
void setLoggingMsgLevel( message::Level level )
{
  if ( !isInitialized() )
  {
    std::cerr << "[ERROR]: slic::initialize() must be called first "
              << "before making any other calls to SLIC.";
    return;
  }
  Logger::getActiveLogger()->setLoggingMsgLevel( level );
}

//------------------------------------------------------------------------------
void setAbortOnError( bool status )
{
  if ( !isInitialized() )
  {
    std::cerr << "[ERROR]: slic::initialize() must be called first "
              << "before making any other calls to SLIC.";
    return;
  }

  Logger::getActiveLogger()->setAbortOnError( status );
}

//------------------------------------------------------------------------------
void enableAbortOnError() {
  setAbortOnError( true );
}

//------------------------------------------------------------------------------
void disableAbortOnError() {
  setAbortOnError( false );
}

//------------------------------------------------------------------------------
bool isAbortOnErrorsEnabled()
{
  if ( !isInitialized() )
  {
    std::cerr << "[ERROR]: slic::initialize() must be called first "
              << "before making any other calls to SLIC.";
    return false;
  }

  return ( Logger::getActiveLogger()->isAbortOnErrorsEnabled() );
}

//------------------------------------------------------------------------------
void setAbortOnWarning( bool status )
{
  if ( !isInitialized() )
  {
    std::cerr << "[ERROR]: slic::initialize() must be called first "
              << "before making any other calls to SLIC.";
    return;
  }

  Logger::getActiveLogger()->setAbortOnWarning( status );
}

//------------------------------------------------------------------------------
void enableAbortOnWarning() {
  setAbortOnWarning(true);
}

//------------------------------------------------------------------------------
void disableAbortOnWarning() {
  setAbortOnWarning(false);
}

//------------------------------------------------------------------------------
bool isAbortOnWarningsEnabled()
{
  if ( !isInitialized() )
  {
    std::cerr << "[ERROR]: slic::initialize() must be called first "
              << "before making any other calls to SLIC.";
    return false;
  }

  return ( Logger::getActiveLogger()->isAbortOnWarningsEnabled() );
}

//------------------------------------------------------------------------------
void addStreamToMsgLevel( LogStream* ls, message::Level level )
{
  if ( !isInitialized() )
  {
    std::cerr << "[ERROR]: slic::initialize() must be called first "
              << "before making any other calls to SLIC.";
    return;
  }
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
  if ( !isInitialized() )
  {
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
  if ( !isInitialized() )
  {
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
  if ( !isInitialized() )
  {
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
  if ( !isInitialized() )
  {
    return;
  }
  Logger::getActiveLogger()->logMessage( level, message, tag, fileName, line,
                                         filter_duplicates );
}

//------------------------------------------------------------------------------
void logErrorMessage( const std::string& message,
                      const std::string& fileName,
                      int line )
{
  std::ostringstream oss;
  oss << message << slic::stacktrace();

  slic::logMessage( message::Error, oss.str(), fileName, line );
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
  if ( !isInitialized() )
  {
    std::cerr << "[ERROR]: slic::initialize() must be called first "
              << "before making any other calls to SLIC.";
    return;
  }
  Logger::getActiveLogger()->flushStreams();
}

//------------------------------------------------------------------------------
void pushStreams()
{
  if ( !isInitialized() )
  {
    std::cerr << "[ERROR]: slic::initialize() must be called first "
              << "before making any other calls to SLIC.";
    return;
  }
  Logger::getActiveLogger()->pushStreams();
}

//------------------------------------------------------------------------------
void finalize()
{
  Logger::finalize();
}

//------------------------------------------------------------------------------
#ifdef WIN32
std::string stacktrace( )
{
  void* stack[10];
  std::ostringstream oss;

  unsigned short frames;
  SYMBOL_INFO* symbol;
  HANDLE process;

  process = GetCurrentProcess();

  SymInitialize( process, NULL, TRUE );

  frames               = CaptureStackBackTrace( 0, 10, stack, NULL );
  symbol               = ( SYMBOL_INFO* )calloc(
    sizeof( SYMBOL_INFO ) + 256 * sizeof( char ), 1 );
  symbol->MaxNameLen   = 255;
  symbol->SizeOfStruct = sizeof( SYMBOL_INFO );

  oss << "\n** StackTrace of " << frames << " frames **\n";
  for(int i = 0 ; i < frames ; i++ )
  {
    char outString[512];
    SymFromAddr( process, ( DWORD64 )( stack[ i ] ), 0, symbol );

    sprintf_s(outString, "%i: %s - 0x%0X", frames - i - 1, symbol->Name,
              symbol->Address );
    oss << outString << std::endl;
  }

  free( symbol );
  oss << "=====\n\n";

  return ( oss.str() );
}

#else
std::string stacktrace( )
{
  void* array[10];
  const int size = backtrace( array, 10 );
  char** strings = backtrace_symbols( array, size );

  std::ostringstream oss;
  oss << "\n** StackTrace of " << size << " frames **\n";
  for ( int i=0 ; i < size ; ++i )
  {
    oss << strings[ i ] << std::endl;
  }
  oss << "=====\n\n";

  free( strings );

  return ( oss.str() );
}
#endif

} /* namespace slic */

} /* namespace axom */

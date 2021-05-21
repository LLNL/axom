// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/slic/interface/slic.hpp"
#include "axom/slic/internal/stacktrace.hpp"

#include <sstream>  // for std::ostringstream

namespace axom
{
namespace slic
{
//------------------------------------------------------------------------------
// Initialize static variables for controlling runtime behavior of asserts and
// error macros.
//------------------------------------------------------------------------------
bool debug::checksAreErrors = false;

void initialize() { Logger::initialize(); }

//------------------------------------------------------------------------------
bool isInitialized() { return (Logger::getActiveLogger() != nullptr); }

//------------------------------------------------------------------------------
void createLogger(const std::string& name, char imask)
{
  if(!isInitialized())
  {
    std::cerr << "[ERROR]: slic::initialize() must be called "
              << "before making any other calls to SLIC.";
    return;
  }
  Logger::createLogger(name, imask);
}

//------------------------------------------------------------------------------
bool activateLogger(const std::string& name)
{
  if(!isInitialized())
  {
    std::cerr << "[ERROR]: slic::initialize() must be called "
              << "before making any other calls to SLIC.";
    return false;
  }
  return Logger::activateLogger(name);
}

//------------------------------------------------------------------------------
std::string getActiveLoggerName()
{
  if(!isInitialized())
  {
    std::cerr << "[ERROR]: slic::initialize() must be called "
              << "before making any other calls to SLIC.";
    return "";
  }
  return (Logger::getActiveLoggerName());
}

//------------------------------------------------------------------------------
message::Level getLoggingMsgLevel()
{
  if(!isInitialized())
  {
    std::cerr << "[ERROR]: slic::initialize() must be called "
              << "before making any other calls to SLIC.";
    return message::Num_Levels;
  }
  return Logger::getActiveLogger()->getLoggingMsgLevel();
}

//------------------------------------------------------------------------------
void setLoggingMsgLevel(message::Level level)
{
  if(!isInitialized())
  {
    std::cerr << "[ERROR]: slic::initialize() must be called "
              << "before making any other calls to SLIC.";
    return;
  }
  Logger::getActiveLogger()->setLoggingMsgLevel(level);
}

//------------------------------------------------------------------------------
void setAbortOnError(bool status)
{
  if(!isInitialized())
  {
    std::cerr << "[ERROR]: slic::initialize() must be called "
              << "before making any other calls to SLIC.";
    return;
  }

  Logger::getActiveLogger()->setAbortOnError(status);
}

//------------------------------------------------------------------------------
void enableAbortOnError() { setAbortOnError(true); }

//------------------------------------------------------------------------------
void disableAbortOnError() { setAbortOnError(false); }

//------------------------------------------------------------------------------
bool isAbortOnErrorsEnabled()
{
  if(!isInitialized())
  {
    std::cerr << "[ERROR]: slic::initialize() must be called "
              << "before making any other calls to SLIC.";
    return false;
  }

  return (Logger::getActiveLogger()->isAbortOnErrorsEnabled());
}

//------------------------------------------------------------------------------
void setAbortOnWarning(bool status)
{
  if(!isInitialized())
  {
    std::cerr << "[ERROR]: slic::initialize() must be called "
              << "before making any other calls to SLIC.";
    return;
  }

  Logger::getActiveLogger()->setAbortOnWarning(status);
}

//------------------------------------------------------------------------------
void enableAbortOnWarning() { setAbortOnWarning(true); }

//------------------------------------------------------------------------------
void disableAbortOnWarning() { setAbortOnWarning(false); }

//------------------------------------------------------------------------------
bool isAbortOnWarningsEnabled()
{
  if(!isInitialized())
  {
    std::cerr << "[ERROR]: slic::initialize() must be called "
              << "before making any other calls to SLIC.";
    return false;
  }

  return (Logger::getActiveLogger()->isAbortOnWarningsEnabled());
}

//------------------------------------------------------------------------------
void setAbortFunction(AbortFunctionPtr abort_func)
{
  if(!isInitialized())
  {
    std::cerr << "[ERROR]: slic::initialize() must be called "
              << "before making any other calls to SLIC.";
    return;
  }

  Logger::getActiveLogger()->setAbortFunction(abort_func);
}

//------------------------------------------------------------------------------
void addStreamToMsgLevel(LogStream* ls, message::Level level)
{
  if(!isInitialized())
  {
    std::cerr << "[ERROR]: slic::initialize() must be called "
              << "before making any other calls to SLIC.";
    return;
  }
  Logger::getActiveLogger()->addStreamToMsgLevel(ls, level);
}

//------------------------------------------------------------------------------
void addStreamToAllMsgLevels(LogStream* ls)
{
  Logger::getActiveLogger()->addStreamToAllMsgLevels(ls);
}

//------------------------------------------------------------------------------
void logMessage(message::Level level,
                const std::string& message,
                bool filter_duplicates)
{
  if(!isInitialized())
  {
    return;
  }
  Logger::getActiveLogger()->logMessage(level, message, filter_duplicates);
}

//------------------------------------------------------------------------------
void logMessage(message::Level level,
                const std::string& message,
                const std::string& tag,
                bool filter_duplicates)
{
  if(!isInitialized())
  {
    return;
  }
  Logger::getActiveLogger()->logMessage(level, message, tag, filter_duplicates);
}

//------------------------------------------------------------------------------
void logMessage(message::Level level,
                const std::string& message,
                const std::string& fileName,
                int line,
                bool filter_duplicates)
{
  if(!isInitialized())
  {
    return;
  }
  Logger::getActiveLogger()->logMessage(level,
                                        message,
                                        fileName,
                                        line,
                                        filter_duplicates);
}

//------------------------------------------------------------------------------
void logMessage(message::Level level,
                const std::string& message,
                const std::string& tag,
                const std::string& fileName,
                int line,
                bool filter_duplicates)
{
  if(!isInitialized())
  {
    return;
  }
  Logger::getActiveLogger()
    ->logMessage(level, message, tag, fileName, line, filter_duplicates);
}

//------------------------------------------------------------------------------
void logErrorMessage(const std::string& message,
                     const std::string& fileName,
                     int line)
{
  std::ostringstream oss;
  oss << message << slic::internal::stacktrace();

  slic::logMessage(message::Error, oss.str(), fileName, line);
}

//------------------------------------------------------------------------------
void logWarningMessage(const std::string& message,
                       const std::string& fileName,
                       int line)
{
  slic::logMessage(message::Warning, message, fileName, line);
}

//------------------------------------------------------------------------------
void flushStreams()
{
  if(!isInitialized())
  {
    std::cerr << "[ERROR]: slic::initialize() must be called "
              << "before making any other calls to SLIC.";
    return;
  }
  Logger::getActiveLogger()->flushStreams();
}

//------------------------------------------------------------------------------
void pushStreams()
{
  if(!isInitialized())
  {
    std::cerr << "[ERROR]: slic::initialize() must be called "
              << "before making any other calls to SLIC.";
    return;
  }
  Logger::getActiveLogger()->pushStreams();
}

//------------------------------------------------------------------------------
void finalize() { Logger::finalize(); }

} /* namespace slic */

} /* namespace axom */

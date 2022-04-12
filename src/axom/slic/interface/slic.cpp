// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
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
void ensureInitialized()
{
  if(!isInitialized())
  {
    initialize();
    disableAbortOnError();
    disableAbortOnWarning();
    setLoggingMsgLevel(slic::message::Debug);
    std::string format = std::string("<TIMESTAMP>\n") +
      std::string("[<LEVEL>]: <MESSAGE> \n") + std::string("FILE=<FILE>\n") +
      std::string("LINE=<LINE>\n\n");
    addStreamToAllMsgLevels(new slic::GenericOutputStream(&std::cout, format));

    logMessage(
      message::Warning,
      "slic::initialize() must be called before any other calls to SLIC\n."
      "The SLIC library called slic::initialize() for you and set up a minimal"
      "configuration\nto allow log messages to print.\n"
      "Please call slic::initialize() near the beginning of the code\n"
      "to fix this error and get rid of this message.\n"
      "Please call slic::finalize() after all other calls to SLIC.\n");
  }
}

//------------------------------------------------------------------------------
void createLogger(const std::string& name, char imask)
{
  ensureInitialized();
  Logger::createLogger(name, imask);
}

//------------------------------------------------------------------------------
bool activateLogger(const std::string& name)
{
  ensureInitialized();
  return Logger::activateLogger(name);
}

//------------------------------------------------------------------------------
std::string getActiveLoggerName()
{
  ensureInitialized();
  return (Logger::getActiveLoggerName());
}

//------------------------------------------------------------------------------
message::Level getLoggingMsgLevel()
{
  ensureInitialized();
  return Logger::getActiveLogger()->getLoggingMsgLevel();
}

//------------------------------------------------------------------------------
void setLoggingMsgLevel(message::Level level)
{
  ensureInitialized();
  Logger::getActiveLogger()->setLoggingMsgLevel(level);
}

//------------------------------------------------------------------------------
void setAbortOnError(bool status)
{
  ensureInitialized();

  Logger::getActiveLogger()->setAbortOnError(status);
}

//------------------------------------------------------------------------------
void enableAbortOnError() { setAbortOnError(true); }

//------------------------------------------------------------------------------
void disableAbortOnError() { setAbortOnError(false); }

//------------------------------------------------------------------------------
bool isAbortOnErrorsEnabled()
{
  ensureInitialized();
  return (Logger::getActiveLogger()->isAbortOnErrorsEnabled());
}

//------------------------------------------------------------------------------
void setAbortOnWarning(bool status)
{
  ensureInitialized();
  Logger::getActiveLogger()->setAbortOnWarning(status);
}

//------------------------------------------------------------------------------
void enableAbortOnWarning() { setAbortOnWarning(true); }

//------------------------------------------------------------------------------
void disableAbortOnWarning() { setAbortOnWarning(false); }

//------------------------------------------------------------------------------
bool isAbortOnWarningsEnabled()
{
  ensureInitialized();
  return (Logger::getActiveLogger()->isAbortOnWarningsEnabled());
}

//------------------------------------------------------------------------------
void setAbortFunction(AbortFunctionPtr abort_func)
{
  ensureInitialized();
  Logger::getActiveLogger()->setAbortFunction(abort_func);
}

//------------------------------------------------------------------------------
void addStreamToMsgLevel(LogStream* ls, message::Level level)
{
  ensureInitialized();
  Logger::getActiveLogger()->addStreamToMsgLevel(ls, level);
}

//------------------------------------------------------------------------------
void addStreamToMsgLevel(GenericOutputStream* ls, message::Level level)
{
  ensureInitialized();
  Logger::getActiveLogger()->addStreamToMsgLevel(ls, level);
}

//------------------------------------------------------------------------------
void addStreamToAllMsgLevels(LogStream* ls)
{
  ensureInitialized();
  Logger::getActiveLogger()->addStreamToAllMsgLevels(ls);
}

//------------------------------------------------------------------------------
void addStreamToAllMsgLevels(GenericOutputStream* ls)
{
  ensureInitialized();
  Logger::getActiveLogger()->addStreamToAllMsgLevels(ls);
}

//------------------------------------------------------------------------------
void logMessage(message::Level level,
                const std::string& message,
                bool filter_duplicates)
{
  ensureInitialized();
  Logger::getActiveLogger()->logMessage(level, message, filter_duplicates);
}

//------------------------------------------------------------------------------
void logMessage(message::Level level,
                const std::string& message,
                const std::string& tag,
                bool filter_duplicates)
{
  ensureInitialized();
  Logger::getActiveLogger()->logMessage(level, message, tag, filter_duplicates);
}

//------------------------------------------------------------------------------
void logMessage(message::Level level,
                const std::string& message,
                const std::string& fileName,
                int line,
                bool filter_duplicates)
{
  ensureInitialized();
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
  ensureInitialized();
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
  ensureInitialized();
  Logger::getActiveLogger()->flushStreams();
}

//------------------------------------------------------------------------------
void pushStreams()
{
  ensureInitialized();
  Logger::getActiveLogger()->pushStreams();
}

//------------------------------------------------------------------------------
void finalize() { Logger::finalize(); }

} /* namespace slic */

} /* namespace axom */

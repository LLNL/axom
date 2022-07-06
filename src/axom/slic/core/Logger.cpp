// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/slic/core/Logger.hpp"

#include "axom/slic/core/LogStream.hpp"

#include "axom/core/utilities/Utilities.hpp"  // for utilities::processAbort()

// C/C++ includes
#include <iostream>  // for std::cout, std::cerr

namespace axom
{
namespace slic
{
using Loggermap = std::map<std::string, Logger*>;

//------------------------------------------------------------------------------
// This is a singleton, scope-limited to this file.
Loggermap& getLoggers()
{
  static Loggermap s_loggers;
  return s_loggers;
}

//------------------------------------------------------------------------------
// This is a singleton, scope-limited to this file.
Logger*& getLogger()
{
  static Logger* s_Logger = nullptr;
  return s_Logger;
}

//------------------------------------------------------------------------------
Logger::Logger()
  : m_abortOnError(true)
  , m_abortOnWarning(false)
  , m_abortFunction(axom::utilities::processAbort)
{
  // by default, all message streams are disabled
  for(int i = 0; i < message::Num_Levels; ++i)
  {
    m_isEnabled[i] = false;
  }
}

//------------------------------------------------------------------------------
Logger::Logger(const std::string& name)
  : m_name(name)
  , m_abortOnError(true)
  , m_abortOnWarning(false)
  , m_abortFunction(axom::utilities::processAbort)
{
  // by default, all message streams are disabled
  for(int i = 0; i < message::Num_Levels; ++i)
  {
    m_isEnabled[i] = false;
  }
}

//------------------------------------------------------------------------------
Logger::~Logger()
{
  std::map<LogStream*, LogStream*>::iterator it = m_streamObjectsManager.begin();
  for(; it != m_streamObjectsManager.end(); ++it)
  {
    delete it->second;
  }  // END for all logStreams

  for(int level = message::Error; level < message::Num_Levels; ++level)
  {
    m_logStreams[level].clear();

  }  // END for all levels
}

//------------------------------------------------------------------------------
void Logger::setAbortFlag(bool val, message::Level level)
{
  if(m_isEnabled[level] == false)
  {
    return;
  }

  // Set flag for the given stream and level
  unsigned nstreams = static_cast<unsigned>(m_logStreams[level].size());
  for(unsigned istream = 0; istream < nstreams; ++istream)
  {
    m_logStreams[level][istream]->setAbortFlag(val);
  }
}

//------------------------------------------------------------------------------
void Logger::abortIfEnabled(message::Level level)
{
  if(this->confirmAbortStreams(level) &&
     ((m_abortOnError && (level == message::Error)) ||
      (m_abortOnWarning && (level == message::Warning))))
  {
    this->flushStreams();
    m_abortFunction();
  }
}

//------------------------------------------------------------------------------
void Logger::setAbortFunction(AbortFunctionPtr abort_func)
{
  if(abort_func == nullptr)
  {
    std::cerr << "WARNING: slic::Logger::setAbortFunction() -- supplied abort "
                 "function is invalid!\n";
    return;
  }

  m_abortFunction = abort_func;
}

//------------------------------------------------------------------------------
message::Level Logger::getLoggingMsgLevel()
{
  int lev = 0;
  const int maxLevels = static_cast<int>(message::Num_Levels);

  // find first disabled level
  while(lev < maxLevels && m_isEnabled[lev])
  {
    ++lev;
  }

  // Subtract one level (with clamp) since we overshot by one
  lev = axom::utilities::clampVal(lev - 1, 0, maxLevels - 1);

  return static_cast<message::Level>(lev);
}

//------------------------------------------------------------------------------
void Logger::setLoggingMsgLevel(message::Level level)
{
  for(int i = 0; i < message::Num_Levels; ++i)
  {
    m_isEnabled[i] = (i <= level) ? true : false;
  }
}

//------------------------------------------------------------------------------
void Logger::addStreamToMsgLevel(LogStream* ls,
                                 message::Level level,
                                 bool pass_ownership)
{
  if(ls == nullptr)
  {
    std::cerr << "WARNING: supplied log stream is NULL!\n";
    return;
  }

  m_logStreams[level].push_back(ls);

  if(pass_ownership)
  {
    m_streamObjectsManager[ls] = ls;
  }
}

//------------------------------------------------------------------------------
void Logger::addStreamToAllMsgLevels(LogStream* ls)
{
  if(ls == nullptr)
  {
    std::cerr << "WARNING: supplied log stream is NULL!\n";
    return;
  }

  for(int level = message::Error; level < message::Num_Levels; ++level)
  {
    this->addStreamToMsgLevel(ls, static_cast<message::Level>(level));

  }  // END for all levels
}

//------------------------------------------------------------------------------
int Logger::getNumStreamsAtMsgLevel(message::Level level)
{
  return static_cast<int>(m_logStreams[level].size());
}

//------------------------------------------------------------------------------
LogStream* Logger::getStream(message::Level level, int i)
{
  if(i < 0 || i >= static_cast<int>(m_logStreams[level].size()))
  {
    std::cerr << "ERROR: stream index is out-of-bounds!\n";
    return nullptr;
  }

  return m_logStreams[level][i];
}

//------------------------------------------------------------------------------
void Logger::logMessage(message::Level level,
                        const std::string& message,
                        bool filter_duplicates)
{
  this->logMessage(level,
                   message,
                   MSG_IGNORE_TAG,
                   MSG_IGNORE_FILE,
                   MSG_IGNORE_LINE,
                   filter_duplicates);
}

//------------------------------------------------------------------------------
void Logger::logMessage(message::Level level,
                        const std::string& message,
                        const std::string& tagName,
                        bool filter_duplicates)
{
  this->logMessage(level,
                   message,
                   tagName,
                   MSG_IGNORE_FILE,
                   MSG_IGNORE_LINE,
                   filter_duplicates);
}

//------------------------------------------------------------------------------
void Logger::logMessage(message::Level level,
                        const std::string& message,
                        const std::string& fileName,
                        int line,
                        bool filter_duplicates)
{
  this->logMessage(level, message, MSG_IGNORE_TAG, fileName, line, filter_duplicates);
}

//------------------------------------------------------------------------------
void Logger::logMessage(message::Level level,
                        const std::string& message,
                        const std::string& tagName,
                        const std::string& fileName,
                        int line,
                        bool filter_duplicates)
{
  if(m_isEnabled[level] == false)
  {
    return;
  }

  unsigned nstreams = static_cast<unsigned>(m_logStreams[level].size());
  for(unsigned istream = 0; istream < nstreams; ++istream)
  {
    m_logStreams[level][istream]
      ->append(level, message, tagName, fileName, line, filter_duplicates);
  }

  if((m_abortOnError && (level == message::Error)) ||
     (m_abortOnWarning && (level == message::Warning)))
  {
    setAbortFlag(true, level);
  }

  if(level == message::Error || level == message::Warning)
  {
    abortIfEnabled(level);
  }
}

//------------------------------------------------------------------------------
void Logger::flushStreams()
{
  for(int level = message::Error; level < message::Num_Levels; ++level)
  {
    unsigned nstreams = static_cast<unsigned>(m_logStreams[level].size());
    for(unsigned istream = 0; istream < nstreams; ++istream)
    {
      m_logStreams[level][istream]->flush();
    }  // END for all streams

  }  // END for all levels
}

bool Logger::confirmAbortStreams(message::Level level)
{
  bool ret = false;

  // This also needs to also consider message::Error and message::Warning, definitely.
  unsigned nstreams = static_cast<unsigned>(m_logStreams[level].size());
  for(unsigned istream = 0; istream < nstreams; ++istream)
  {
    ret |= m_logStreams[level][istream]->confirmAbort();
  }  // END for all streams

  return ret;
}

//------------------------------------------------------------------------------
void Logger::pushStreams()
{
  for(int level = message::Error; level < message::Num_Levels; ++level)
  {
    unsigned nstreams = static_cast<unsigned>(m_logStreams[level].size());
    for(unsigned istream = 0; istream < nstreams; ++istream)
    {
      m_logStreams[level][istream]->push();

    }  // END for all streams

  }  // END for all levels
}

//------------------------------------------------------------------------------
//                Static Method Implementations
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void Logger::initialize()
{
  Logger::createLogger("root");
  Logger::activateLogger("root");
}

//------------------------------------------------------------------------------
std::map<std::string, Logger*>& Logger::getSLoggers()
{
  static std::map<std::string, Logger*>* s_loggers =
    new std::map<std::string, Logger*>();
  return *s_loggers;
}

//------------------------------------------------------------------------------
bool Logger::createLogger(const std::string& name, char imask)
{
  Loggermap& loggers = getLoggers();
  if(loggers.find(name) != loggers.end())
  {
    std::cerr << "ERROR: " << name << " logger is duplicated!\n";
    return false;
  }

  loggers[name] = new Logger(name);

  if(imask == inherit::nothing)
  {
    /* short-circuit */
    return true;
  }  // END if inherit nothing

  Logger* rootLogger = Logger::getRootLogger();
  if(rootLogger == nullptr)
  {
    std::cerr << "ERROR: no root logger found!\n";
    return false;
  }

  for(int level = message::Error; level < message::Num_Levels; ++level)
  {
    message::Level current_level = static_cast<message::Level>(level);

    int nstreams = rootLogger->getNumStreamsAtMsgLevel(current_level);
    if(nstreams == 0)
    {
      continue;
    }

    if(imask & inherit::masks[level])
    {
      for(int istream = 0; istream < nstreams; ++istream)
      {
        loggers[name]->addStreamToMsgLevel(
          rootLogger->getStream(current_level, istream),
          current_level,
          /* pass_ownership */ false);

      }  // END for all streams at this level

    }  // END if inherit streams at the given level

  }  // END for all message levels

  return true;
}

//------------------------------------------------------------------------------
bool Logger::activateLogger(const std::string& name)
{
  Loggermap& loggers = getLoggers();
  if(loggers.find(name) != loggers.end())
  {
    getLogger() = loggers[name];
    return true;
  }

  return false;
}

//------------------------------------------------------------------------------
void Logger::finalize()
{
  Loggermap& loggers = getLoggers();
  for(Loggermap::iterator it = loggers.begin(); it != loggers.end(); ++it)
  {
    it->second->flushStreams();
  }

  for(Loggermap::iterator it = loggers.begin(); it != loggers.end(); ++it)
  {
    delete it->second;
  }

  loggers.clear();

  getLogger() = nullptr;
}

//------------------------------------------------------------------------------
std::string Logger::getActiveLoggerName() { return getLogger()->getName(); }

//------------------------------------------------------------------------------
Logger* Logger::getActiveLogger() { return getLogger(); }

//------------------------------------------------------------------------------
Logger* Logger::getRootLogger()
{
  Loggermap& loggers = getLoggers();
  if(loggers.find("root") == loggers.end())
  {
    // no root logger
    return nullptr;
  }

  return (loggers["root"]);
}

} /* namespace slic */

} /* namespace axom */

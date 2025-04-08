// wrapSLIC.cpp
// This file is generated by Shroud 0.13.0. Do not edit.
//
// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/slic/interface/slic.hpp"
#include <string>
#include "axom/slic/streams/GenericOutputStream.hpp"
#include <cstring>
#include "wrapSLIC.h"

// splicer begin CXX_definitions
// splicer end CXX_definitions

extern "C" {

// helper ShroudLenTrim
// Returns the length of character string src with length nsrc,
// ignoring any trailing blanks.
static int ShroudLenTrim(const char *src, int nsrc)
{
  int i;

  for(i = nsrc - 1; i >= 0; i--)
  {
    if(src[i] != ' ')
    {
      break;
    }
  }

  return i + 1;
}

// helper ShroudStrCopy
// Copy src into dest, blank fill to ndest characters
// Truncate if dest is too short.
// dest will not be NULL terminated.
static void ShroudStrCopy(char *dest, int ndest, const char *src, int nsrc)
{
  if(src == NULL)
  {
    std::memset(dest, ' ', ndest);  // convert NULL pointer to blank filled string
  }
  else
  {
    if(nsrc < 0) nsrc = std::strlen(src);
    int nm = nsrc < ndest ? nsrc : ndest;
    std::memcpy(dest, src, nm);
    if(ndest > nm) std::memset(dest + nm, ' ', ndest - nm);  // blank fill
  }
}
// splicer begin C_definitions
// splicer end C_definitions

void SLIC_initialize(void)
{
  // splicer begin function.initialize
  axom::slic::initialize();
  // splicer end function.initialize
}

bool SLIC_isInitialized(void)
{
  // splicer begin function.isInitialized
  bool SHC_rv = axom::slic::isInitialized();
  return SHC_rv;
  // splicer end function.isInitialized
}

void SLIC_createLogger(const char *name, char imask)
{
  // splicer begin function.createLogger
  const std::string SHCXX_name(name);
  axom::slic::createLogger(SHCXX_name, imask);
  // splicer end function.createLogger
}

void SLIC_createLogger_bufferify(char *name, int SHT_name_len, char imask)
{
  // splicer begin function.createLogger_bufferify
  const std::string SHCXX_name(name, ShroudLenTrim(name, SHT_name_len));
  axom::slic::createLogger(SHCXX_name, imask);
  // splicer end function.createLogger_bufferify
}

bool SLIC_activateLogger(const char *name)
{
  // splicer begin function.activateLogger
  const std::string SHCXX_name(name);
  bool SHC_rv = axom::slic::activateLogger(SHCXX_name);
  return SHC_rv;
  // splicer end function.activateLogger
}

bool SLIC_activateLogger_bufferify(char *name, int SHT_name_len)
{
  // splicer begin function.activateLogger_bufferify
  const std::string SHCXX_name(name, ShroudLenTrim(name, SHT_name_len));
  bool SHC_rv = axom::slic::activateLogger(SHCXX_name);
  return SHC_rv;
  // splicer end function.activateLogger_bufferify
}

void SLIC_getActiveLoggerName_bufferify(char *name, int SHT_name_len)
{
  // splicer begin function.getActiveLoggerName_bufferify
  std::string SHCXX_rv = axom::slic::getActiveLoggerName();
  if(SHCXX_rv.empty())
  {
    ShroudStrCopy(name, SHT_name_len, nullptr, 0);
  }
  else
  {
    ShroudStrCopy(name, SHT_name_len, SHCXX_rv.data(), SHCXX_rv.size());
  }
  // splicer end function.getActiveLoggerName_bufferify
}

int SLIC_getLoggingMsgLevel(void)
{
  // splicer begin function.getLoggingMsgLevel
  axom::slic::message::Level SHCXX_rv = axom::slic::getLoggingMsgLevel();
  int SHC_rv = static_cast<int>(SHCXX_rv);
  return SHC_rv;
  // splicer end function.getLoggingMsgLevel
}

void SLIC_setLoggingMsgLevel(int level)
{
  // splicer begin function.setLoggingMsgLevel
  axom::slic::message::Level SHCXX_level = static_cast<axom::slic::message::Level>(level);
  axom::slic::setLoggingMsgLevel(SHCXX_level);
  // splicer end function.setLoggingMsgLevel
}

void SLIC_addStreamToMsgLevel(SLIC_GenericOutputStream *ls, int level)
{
  // splicer begin function.addStreamToMsgLevel
  axom::slic::GenericOutputStream *SHCXX_ls =
    static_cast<axom::slic::GenericOutputStream *>(ls->addr);
  axom::slic::message::Level SHCXX_level = static_cast<axom::slic::message::Level>(level);
  axom::slic::addStreamToMsgLevel(SHCXX_ls, SHCXX_level);
  // splicer end function.addStreamToMsgLevel
}

void SLIC_addStreamToAllMsgLevels(SLIC_GenericOutputStream *ls)
{
  // splicer begin function.addStreamToAllMsgLevels
  axom::slic::GenericOutputStream *SHCXX_ls =
    static_cast<axom::slic::GenericOutputStream *>(ls->addr);
  axom::slic::addStreamToAllMsgLevels(SHCXX_ls);
  // splicer end function.addStreamToAllMsgLevels
}

void SLIC_setAbortOnError(bool status)
{
  // splicer begin function.setAbortOnError
  axom::slic::setAbortOnError(status);
  // splicer end function.setAbortOnError
}

void SLIC_enableAbortOnError(void)
{
  // splicer begin function.enableAbortOnError
  axom::slic::enableAbortOnError();
  // splicer end function.enableAbortOnError
}

void SLIC_disableAbortOnError(void)
{
  // splicer begin function.disableAbortOnError
  axom::slic::disableAbortOnError();
  // splicer end function.disableAbortOnError
}

bool SLIC_isAbortOnErrorsEnabled(void)
{
  // splicer begin function.isAbortOnErrorsEnabled
  bool SHC_rv = axom::slic::isAbortOnErrorsEnabled();
  return SHC_rv;
  // splicer end function.isAbortOnErrorsEnabled
}

void SLIC_setAbortOnWarning(bool status)
{
  // splicer begin function.setAbortOnWarning
  axom::slic::setAbortOnWarning(status);
  // splicer end function.setAbortOnWarning
}

void SLIC_enableAbortOnWarning(void)
{
  // splicer begin function.enableAbortOnWarning
  axom::slic::enableAbortOnWarning();
  // splicer end function.enableAbortOnWarning
}

void SLIC_disableAbortOnWarning(void)
{
  // splicer begin function.disableAbortOnWarning
  axom::slic::disableAbortOnWarning();
  // splicer end function.disableAbortOnWarning
}

bool SLIC_isAbortOnWarningsEnabled(void)
{
  // splicer begin function.isAbortOnWarningsEnabled
  bool SHC_rv = axom::slic::isAbortOnWarningsEnabled();
  return SHC_rv;
  // splicer end function.isAbortOnWarningsEnabled
}

void SLIC_logMessage_file_line(int level, const char *message, const char *fileName, int line)
{
  // splicer begin function.logMessage_file_line
  axom::slic::message::Level SHCXX_level = static_cast<axom::slic::message::Level>(level);
  const std::string SHCXX_message(message);
  const std::string SHCXX_fileName(fileName);
  axom::slic::logMessage(SHCXX_level, SHCXX_message, SHCXX_fileName, line);
  // splicer end function.logMessage_file_line
}

void SLIC_logMessage_file_line_bufferify(int level,
                                         char *message,
                                         int SHT_message_len,
                                         char *fileName,
                                         int SHT_fileName_len,
                                         int line)
{
  // splicer begin function.logMessage_file_line_bufferify
  axom::slic::message::Level SHCXX_level = static_cast<axom::slic::message::Level>(level);
  const std::string SHCXX_message(message, ShroudLenTrim(message, SHT_message_len));
  const std::string SHCXX_fileName(fileName, ShroudLenTrim(fileName, SHT_fileName_len));
  axom::slic::logMessage(SHCXX_level, SHCXX_message, SHCXX_fileName, line);
  // splicer end function.logMessage_file_line_bufferify
}

void SLIC_logMessage_file_line_filter(int level,
                                      const char *message,
                                      const char *fileName,
                                      int line,
                                      bool filter_duplicates)
{
  // splicer begin function.logMessage_file_line_filter
  axom::slic::message::Level SHCXX_level = static_cast<axom::slic::message::Level>(level);
  const std::string SHCXX_message(message);
  const std::string SHCXX_fileName(fileName);
  axom::slic::logMessage(SHCXX_level, SHCXX_message, SHCXX_fileName, line, filter_duplicates);
  // splicer end function.logMessage_file_line_filter
}

void SLIC_logMessage_file_line_filter_bufferify(int level,
                                                char *message,
                                                int SHT_message_len,
                                                char *fileName,
                                                int SHT_fileName_len,
                                                int line,
                                                bool filter_duplicates)
{
  // splicer begin function.logMessage_file_line_filter_bufferify
  axom::slic::message::Level SHCXX_level = static_cast<axom::slic::message::Level>(level);
  const std::string SHCXX_message(message, ShroudLenTrim(message, SHT_message_len));
  const std::string SHCXX_fileName(fileName, ShroudLenTrim(fileName, SHT_fileName_len));
  axom::slic::logMessage(SHCXX_level, SHCXX_message, SHCXX_fileName, line, filter_duplicates);
  // splicer end function.logMessage_file_line_filter_bufferify
}

void SLIC_logMessage(int level, const char *message)
{
  // splicer begin function.logMessage
  axom::slic::message::Level SHCXX_level = static_cast<axom::slic::message::Level>(level);
  const std::string SHCXX_message(message);
  axom::slic::logMessage(SHCXX_level, SHCXX_message);
  // splicer end function.logMessage
}

void SLIC_logMessage_bufferify(int level, char *message, int SHT_message_len)
{
  // splicer begin function.logMessage_bufferify
  axom::slic::message::Level SHCXX_level = static_cast<axom::slic::message::Level>(level);
  const std::string SHCXX_message(message, ShroudLenTrim(message, SHT_message_len));
  axom::slic::logMessage(SHCXX_level, SHCXX_message);
  // splicer end function.logMessage_bufferify
}

void SLIC_logMessage_filter(int level, const char *message, bool filter_duplicates)
{
  // splicer begin function.logMessage_filter
  axom::slic::message::Level SHCXX_level = static_cast<axom::slic::message::Level>(level);
  const std::string SHCXX_message(message);
  axom::slic::logMessage(SHCXX_level, SHCXX_message, filter_duplicates);
  // splicer end function.logMessage_filter
}

void SLIC_logMessage_filter_bufferify(int level,
                                      char *message,
                                      int SHT_message_len,
                                      bool filter_duplicates)
{
  // splicer begin function.logMessage_filter_bufferify
  axom::slic::message::Level SHCXX_level = static_cast<axom::slic::message::Level>(level);
  const std::string SHCXX_message(message, ShroudLenTrim(message, SHT_message_len));
  axom::slic::logMessage(SHCXX_level, SHCXX_message, filter_duplicates);
  // splicer end function.logMessage_filter_bufferify
}

void SLIC_flushStreams(void)
{
  // splicer begin function.flushStreams
  axom::slic::flushStreams();
  // splicer end function.flushStreams
}

void SLIC_finalize(void)
{
  // splicer begin function.finalize
  axom::slic::finalize();
  // splicer end function.finalize
}

}  // extern "C"

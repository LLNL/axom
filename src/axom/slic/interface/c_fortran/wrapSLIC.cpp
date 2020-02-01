// wrapSLIC.cpp
// This is generated code, do not edit
//
// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "wrapSLIC.h"
#include <cstring>
#include <stdlib.h>
#include <string>
#include "axom/slic/interface/slic.hpp"
#include "typesSLIC.h"

// splicer begin CXX_definitions
// splicer end CXX_definitions

extern "C" {


// helper function
// Copy src into dest, blank fill to ndest characters
// Truncate if dest is too short.
// dest will not be NULL terminated.
static void ShroudStrCopy(char* dest, int ndest, const char* src, int nsrc)
{
  int nm = nsrc < ndest ? nsrc : ndest;
  std::memcpy(dest,src,nm);
  if(ndest > nm)
    std::memset(dest+nm,' ',ndest-nm);
}
// splicer begin C_definitions
// splicer end C_definitions

void SLIC_initialize()
{
// splicer begin function.initialize
  axom::slic::initialize();
  return;
// splicer end function.initialize
}

bool SLIC_is_initialized()
{
// splicer begin function.is_initialized
  bool SHC_rv = axom::slic::isInitialized();
  return SHC_rv;
// splicer end function.is_initialized
}

void SLIC_create_logger(const char* name, char imask)
{
// splicer begin function.create_logger
  const std::string SH_name(name);
  axom::slic::createLogger(SH_name, imask);
  return;
// splicer end function.create_logger
}

void SLIC_create_logger_bufferify(const char* name, int Lname, char imask)
{
// splicer begin function.create_logger_bufferify
  const std::string SH_name(name, Lname);
  axom::slic::createLogger(SH_name, imask);
  return;
// splicer end function.create_logger_bufferify
}

bool SLIC_activate_logger(const char* name)
{
// splicer begin function.activate_logger
  const std::string SH_name(name);
  bool SHC_rv = axom::slic::activateLogger(SH_name);
  return SHC_rv;
// splicer end function.activate_logger
}

bool SLIC_activate_logger_bufferify(const char* name, int Lname)
{
// splicer begin function.activate_logger_bufferify
  const std::string SH_name(name, Lname);
  bool SHC_rv = axom::slic::activateLogger(SH_name);
  return SHC_rv;
// splicer end function.activate_logger_bufferify
}

void SLIC_get_active_logger_name_bufferify(char* name, int Nname)
{
// splicer begin function.get_active_logger_name_bufferify
  std::string SHCXX_rv = axom::slic::getActiveLoggerName();
  if (SHCXX_rv.empty())
  {
    std::memset(name, ' ', Nname);
  }
  else
  {
    ShroudStrCopy(name, Nname, SHCXX_rv.data(), SHCXX_rv.size());
  }
  return;
// splicer end function.get_active_logger_name_bufferify
}

int SLIC_get_logging_msg_level()
{
// splicer begin function.get_logging_msg_level
  axom::slic::message::Level SHCXX_rv = axom::slic::getLoggingMsgLevel();
  int SHC_rv = static_cast<int>(SHCXX_rv);
  return SHC_rv;
// splicer end function.get_logging_msg_level
}

void SLIC_set_logging_msg_level(int level)
{
// splicer begin function.set_logging_msg_level
  axom::slic::message::Level SHCXX_level =
    static_cast<axom::slic::message::Level>(level);
  axom::slic::setLoggingMsgLevel(SHCXX_level);
  return;
// splicer end function.set_logging_msg_level
}

void SLIC_set_abort_on_error(bool status)
{
// splicer begin function.set_abort_on_error
  axom::slic::setAbortOnError(status);
  return;
// splicer end function.set_abort_on_error
}

void SLIC_enable_abort_on_error()
{
// splicer begin function.enable_abort_on_error
  axom::slic::enableAbortOnError();
  return;
// splicer end function.enable_abort_on_error
}

void SLIC_disable_abort_on_error()
{
// splicer begin function.disable_abort_on_error
  axom::slic::disableAbortOnError();
  return;
// splicer end function.disable_abort_on_error
}

bool SLIC_is_abort_on_errors_enabled()
{
// splicer begin function.is_abort_on_errors_enabled
  bool SHC_rv = axom::slic::isAbortOnErrorsEnabled();
  return SHC_rv;
// splicer end function.is_abort_on_errors_enabled
}

void SLIC_set_abort_on_warning(bool status)
{
// splicer begin function.set_abort_on_warning
  axom::slic::setAbortOnWarning(status);
  return;
// splicer end function.set_abort_on_warning
}

void SLIC_enable_abort_on_warning()
{
// splicer begin function.enable_abort_on_warning
  axom::slic::enableAbortOnWarning();
  return;
// splicer end function.enable_abort_on_warning
}

void SLIC_disable_abort_on_warning()
{
// splicer begin function.disable_abort_on_warning
  axom::slic::disableAbortOnWarning();
  return;
// splicer end function.disable_abort_on_warning
}

bool SLIC_is_abort_on_warnings_enabled()
{
// splicer begin function.is_abort_on_warnings_enabled
  bool SHC_rv = axom::slic::isAbortOnWarningsEnabled();
  return SHC_rv;
// splicer end function.is_abort_on_warnings_enabled
}

void SLIC_log_message(int level, const char* message, const char* fileName,
                      int line, bool filter)
{
// splicer begin function.log_message
  axom::slic::message::Level SHCXX_level =
    static_cast<axom::slic::message::Level>(level);
  const std::string SH_message(message);
  const std::string SH_fileName(fileName);
  axom::slic::logMessage(SHCXX_level, SH_message, SH_fileName, line, filter);
  return;
// splicer end function.log_message
}

void SLIC_log_message_bufferify(int level, const char* message, int Lmessage,
                                const char* fileName, int LfileName, int line,
                                bool filter)
{
// splicer begin function.log_message_bufferify
  axom::slic::message::Level SHCXX_level =
    static_cast<axom::slic::message::Level>(level);
  const std::string SH_message(message, Lmessage);
  const std::string SH_fileName(fileName, LfileName);
  axom::slic::logMessage(SHCXX_level, SH_message, SH_fileName, line, filter);
  return;
// splicer end function.log_message_bufferify
}

void SLIC_finalize()
{
// splicer begin function.finalize
  axom::slic::finalize();
  return;
// splicer end function.finalize
}

// Release C++ allocated memory.
void SLIC_SHROUD_memory_destructor(SLI_SHROUD_capsule_data* cap)
{
  cap->addr = NULL;
  cap->idtor = 0;    // avoid deleting again
}

}  // extern "C"

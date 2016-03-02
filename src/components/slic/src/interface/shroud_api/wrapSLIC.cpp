// wrapSLIC.cpp
// This is generated code, do not edit
//
// Copyright (c) 2015, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory.
//
// All rights reserved.
//
// This source code cannot be distributed without permission and
// further review from Lawrence Livermore National Laboratory.
//
// wrapSLIC.cpp
#include "wrapSLIC.h"
#include <string>
#include "shroudrt.hpp"
#include "slic/slic.hpp"

extern "C" {
namespace asctoolkit {
namespace slic {

void SLIC_initialize()
{
// splicer begin function.initialize
initialize();
return;
// splicer end function.initialize
}

bool SLIC_is_initialized()
{
// splicer begin function.is_initialized
bool rv = isInitialized();
return rv;
// splicer end function.is_initialized
}

void SLIC_finalize()
{
// splicer begin function.finalize
finalize();
return;
// splicer end function.finalize
}

void SLIC_create_logger(const char * name, char imask)
{
// splicer begin function.create_logger
std::string SH_name(name);
createLogger(SH_name, imask);
return;
// splicer end function.create_logger
}

void SLIC_create_logger_bufferify(const char * name, int Lname, char imask)
{
// splicer begin function.create_logger_bufferify
std::string SH_name(name, Lname);
createLogger(SH_name, imask);
return;
// splicer end function.create_logger_bufferify
}

void SLIC_activate_logger(const char * name)
{
// splicer begin function.activate_logger
std::string SH_name(name);
activateLogger(SH_name);
return;
// splicer end function.activate_logger
}

void SLIC_activate_logger_bufferify(const char * name, int Lname)
{
// splicer begin function.activate_logger_bufferify
std::string SH_name(name, Lname);
activateLogger(SH_name);
return;
// splicer end function.activate_logger_bufferify
}

void SLIC_get_active_logger_name_bufferify(char * name, int Lname)
{
// splicer begin function.get_active_logger_name_bufferify
std::string rv = getActiveLoggerName();
asctoolkit::shroud::FccCopy(name, Lname, rv.c_str());
return;
// splicer end function.get_active_logger_name_bufferify
}

void SLIC_set_logging_msg_level(int level)
{
// splicer begin function.set_logging_msg_level
setLoggingMsgLevel(static_cast<message::Level>(level));
return;
// splicer end function.set_logging_msg_level
}

void SLIC_log_message(int level, const char * message, const char * fileName, int line, bool filter)
{
// splicer begin function.log_message
std::string SH_message(message);
std::string SH_fileName(fileName);
logMessage(static_cast<message::Level>(level), SH_message, SH_fileName, line, filter);
return;
// splicer end function.log_message
}

void SLIC_log_message_bufferify(int level, const char * message, int Lmessage, const char * fileName, int LfileName, int line, bool filter)
{
// splicer begin function.log_message_bufferify
std::string SH_message(message, Lmessage);
std::string SH_fileName(fileName, LfileName);
logMessage(static_cast<message::Level>(level), SH_message, SH_fileName, line, filter);
return;
// splicer end function.log_message_bufferify
}

// splicer begin additional_functions
// splicer end additional_functions

}  // namespace asctoolkit
}  // namespace slic
}  // extern "C"

// wrapSLIC.h
// This file is generated by Shroud 0.12.2. Do not edit.
//
// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
/**
 * \file wrapSLIC.h
 * \brief Shroud generated wrapper for slic namespace
 */
// For C users and C++ implementation

#ifndef WRAPSLIC_H
#define WRAPSLIC_H

#include "typesSLIC.h"
#ifndef __cplusplus
#include <stdbool.h>
#endif

// splicer begin CXX_declarations
// splicer end CXX_declarations

#ifdef __cplusplus
extern "C" {
#endif

//  axom::slic::message::Level
enum SLIC_message_Level {
    SLIC_message_Error,
    SLIC_message_Warning,
    SLIC_message_Info,
    SLIC_message_Debug,
    SLIC_message_Num_Levels
};

// splicer begin C_declarations
// splicer end C_declarations

void SLIC_initialize(void);

bool SLIC_is_initialized(void);

void SLIC_create_logger(const char * name, char imask);

void SLIC_create_logger_bufferify(const char * name, int Lname, char imask);

bool SLIC_activate_logger(const char * name);

bool SLIC_activate_logger_bufferify(const char * name, int Lname);

void SLIC_get_active_logger_name_bufferify(char * name, int Nname);

int SLIC_get_logging_msg_level(void);

void SLIC_set_logging_msg_level(int level);

void SLIC_add_stream_to_all_msg_levels(SLIC_LogStream * ls);

void SLIC_set_abort_on_error(bool status);

void SLIC_enable_abort_on_error(void);

void SLIC_disable_abort_on_error(void);

bool SLIC_is_abort_on_errors_enabled(void);

void SLIC_set_abort_on_warning(bool status);

void SLIC_enable_abort_on_warning(void);

void SLIC_disable_abort_on_warning(void);

bool SLIC_is_abort_on_warnings_enabled(void);

void SLIC_log_message(int level, const char * message, const char * fileName, int line, bool filter);

void SLIC_log_message_bufferify(int level, const char * message, int Lmessage, const char * fileName, int LfileName, int line, bool filter);

void SLIC_finalize(void);

#ifdef __cplusplus
}
#endif

#endif  // WRAPSLIC_H

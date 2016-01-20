// wrapstrings.h
// This is generated code, do not edit
/**
 * \file wrapstrings.h
 * \brief Shroud generated wrapper for strings library
 */
// For C users and C++ implementation

#ifndef WRAPSTRINGS_H
#define WRAPSTRINGS_H

#ifdef __cplusplus
extern "C" {
#endif

// declaration of wrapped types

// splicer begin C_definition
// splicer end C_definition

const char * STR_get_name();

void STR_get_name_bufferify(char * output, int Loutput);

void STR_accept_string_const_reference(const char * arg1);

void STR_accept_string_const_reference_bufferify(const char * arg1, int Larg1);

#ifdef __cplusplus
}
#endif

#endif  // WRAPSTRINGS_H

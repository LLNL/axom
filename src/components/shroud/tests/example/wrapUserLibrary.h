// wrapUserLibrary.h
// This is generated code, do not edit
// blah blah
// yada yada
//
/**
 * \file wrapUserLibrary.h
 * \brief Shroud generated wrapper for UserLibrary library
 */
// For C users and C++ implementation

#ifndef WRAPUSERLIBRARY_H
#define WRAPUSERLIBRARY_H

#include "mpi.h"
#include "sidre/wrapGroup.h"

// splicer begin CXX_declarations
// splicer end CXX_declarations

#ifdef __cplusplus
extern "C" {
#endif

// declaration of wrapped types
struct s_SIDRE_group;
typedef struct s_SIDRE_group SIDRE_group;

// splicer begin C_declarations
// splicer end C_declarations

void AA_local_function1();

bool AA_is_name_valid(const char * name);

bool AA_is_name_valid_bufferify(const char * name, int Lname);

bool AA_is_initialized();

void AA_check_bool(bool arg1, bool * arg2, bool * arg3);

void AA_test_names(const char * name);

void AA_test_names_bufferify(const char * name, int Lname);

void AA_test_names_flag(const char * name, int flag);

void AA_test_names_flag_bufferify(const char * name, int Lname, int flag);

void AA_testoptional_0();

void AA_testoptional_1(int i);

void AA_testoptional_2(int i, long j);

void AA_testmpi(MPI_Fint comm);

void AA_testgroup1(SIDRE_group * grp);

void AA_testgroup2(const SIDRE_group * grp);

#ifdef __cplusplus
}
#endif

#endif  // WRAPUSERLIBRARY_H

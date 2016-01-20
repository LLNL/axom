// wrapExClass1.h
// This is generated code, do not edit
// blah blah
// yada yada
//
/**
 * \file wrapExClass1.h
 * \brief Shroud generated wrapper for ExClass1 class
 */
// For C users and C++ implementation

#ifndef WRAPEXCLASS1_H
#define WRAPEXCLASS1_H

#ifdef __cplusplus
extern "C" {
#endif

// declaration of wrapped types
struct s_AA_exclass1;
typedef struct s_AA_exclass1 AA_exclass1;
struct s_AA_exclass2;
typedef struct s_AA_exclass2 AA_exclass2;

// splicer begin class.ExClass1.C_definition
// splicer end class.ExClass1.C_definition

AA_exclass1 * AA_exclass1_new(const char * name);

AA_exclass1 * AA_exclass1_new_bufferify(const char * name, int Lname);

void AA_exclass1_delete(AA_exclass1 * self);

int AA_exclass1_increment_count(AA_exclass1 * self, int incr);

const char * AA_exclass1_get_name(const AA_exclass1 * self);

int AA_exclass1_get_name_length(AA_exclass1 * self);

const char * AA_exclass1_get_name_error_check(const AA_exclass1 * self);

const char * AA_exclass1_get_name_arg(const AA_exclass1 * self);

void AA_exclass1_get_name_arg_bufferify(AA_exclass1 * self, char * name, int Lname);

AA_exclass2 * AA_exclass1_get_root(AA_exclass1 * self);

int AA_exclass1_get_value_from_int(AA_exclass1 * self, int value);

long AA_exclass1_get_value_1(AA_exclass1 * self, long value);

void * AA_exclass1_get_addr(AA_exclass1 * self);

bool AA_exclass1_has_addr(AA_exclass1 * self, bool in);

void AA_exclass1_splicer_special(AA_exclass1 * self);

#ifdef __cplusplus
}
#endif

#endif  // WRAPEXCLASS1_H

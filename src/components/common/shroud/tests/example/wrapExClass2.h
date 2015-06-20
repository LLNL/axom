// wrapExClass2.h
// blah blah
// yada yada
//
// wrapExClass2.h
// For C users and C++ implementation

#ifndef WRAPEXCLASS2_H
#define WRAPEXCLASS2_H

#include "sidre/SidreTypes.h"

#ifdef __cplusplus
extern "C" {
#endif

// declaration of wrapped types
#ifdef EXAMPLE_WRAPPER_IMPL
typedef void AA_exclass1;
typedef void AA_exclass2;
#else
struct s_AA_exclass1;
typedef struct s_AA_exclass1 AA_exclass1;
struct s_AA_exclass2;
typedef struct s_AA_exclass2 AA_exclass2;
#endif

AA_exclass2 * AA_exclass2_ex_class2(const char * name);

void AA_exclass2_ex_class1(AA_exclass2 * self);

const char * AA_exclass2_get_name(const AA_exclass2 * self);

const int AA_exclass2_get_name_length(const AA_exclass2 * self);

AA_exclass1 * AA_exclass2_get_class1(AA_exclass2 * self, AA_exclass1 * in);

void AA_exclass2_declare(AA_exclass2 * self, int type, ATK_SidreLength len);

void AA_exclass2_destroyall(AA_exclass2 * self);

int AA_exclass2_get_type_id(AA_exclass2 * self);

#ifdef __cplusplus
}
#endif

#endif  // WRAPEXCLASS2_H

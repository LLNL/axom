// blah blah
// yada yada
//
// wrapUserLibrary.h
// For C users and C++ implementation

#ifndef WRAPUSERLIBRARY_H
#define WRAPUSERLIBRARY_H

#ifdef __cplusplus
extern "C" {
#endif

// declaration of wrapped types
#ifdef EXAMPLE_WRAPPER_IMPL
#else
#endif

void AA_local_function1();

bool AA_is_name_valid(const char * name);

#ifdef __cplusplus
}
#endif

#endif  // WRAPUSERLIBRARY_H

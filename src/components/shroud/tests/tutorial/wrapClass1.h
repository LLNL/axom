// wrapClass1.h
// This is generated code, do not edit
/**
 * \file wrapClass1.h
 * \brief Shroud generated wrapper for Class1 class
 */
// For C users and C++ implementation

#ifndef WRAPCLASS1_H
#define WRAPCLASS1_H

// splicer begin class.Class1.CXX_declarations
// splicer end class.Class1.CXX_declarations

#ifdef __cplusplus
extern "C" {
#endif

// declaration of wrapped types
struct s_TUT_class1;
typedef struct s_TUT_class1 TUT_class1;

// splicer begin class.Class1.C_declarations
// splicer end class.Class1.C_declarations

TUT_class1 * TUT_class1_new();

void TUT_class1_delete(TUT_class1 * self);

void TUT_class1_method1(TUT_class1 * self);

#ifdef __cplusplus
}
#endif

#endif  // WRAPCLASS1_H

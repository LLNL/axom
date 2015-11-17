// wrapTutorial.h
// This is generated code, do not edit
// wrapTutorial.h
// For C users and C++ implementation

#ifndef WRAPTUTORIAL_H
#define WRAPTUTORIAL_H

#ifdef __cplusplus
extern "C" {
#endif

// declaration of wrapped types

// splicer begin C_definition
// splicer end C_definition

void TUT_function1();

double TUT_function2(double arg1, int arg2);

int TUT_sum(int len, int * values);

bool TUT_function3(bool arg);

const char * TUT_function4a(const char * arg1, const char * arg2);

const char * TUT_function4a_bufferify(const char * arg1, int Larg1, const char * arg2, int Larg2);

const char * TUT_function4b(const char * arg1, const char * arg2);

const char * TUT_function4b_bufferify(const char * arg1, int Larg1, const char * arg2, int Larg2);

double TUT_function5(double arg1, int arg2);

void TUT_function6_from_name(const char * name);

void TUT_function6_bufferify(const char * name, int Lname);

void TUT_function6_from_index(int indx);

void TUT_function7_int(int arg);

void TUT_function7_double(double arg);

int TUT_function8_int();

double TUT_function8_double();

void TUT_function9(double arg);

const char * TUT_last_function_called();

#ifdef __cplusplus
}
#endif

#endif  // WRAPTUTORIAL_H

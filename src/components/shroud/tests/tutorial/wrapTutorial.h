// wrapTutorial.h
// This is generated code, do not edit
/**
 * \file wrapTutorial.h
 * \brief Shroud generated wrapper for Tutorial library
 */
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

void TUT_sum(int len, int * values, int * result);

bool TUT_function3(bool arg);

const char * TUT_function4a(const char * arg1, const char * arg2);

const char * TUT_function4a_bufferify(const char * arg1, int Larg1, const char * arg2, int Larg2);

void TUT_function4b_bufferify(const char * arg1, int Larg1, const char * arg2, int Larg2, char * output, int Loutput);

double TUT_function5();

double TUT_function5_arg1(double arg1);

double TUT_function5_arg1_arg2(double arg1, bool arg2);

void TUT_function6_from_name(const char * name);

void TUT_function6_from_name_bufferify(const char * name, int Lname);

void TUT_function6_from_index(int indx);

void TUT_function7_int(int arg);

void TUT_function7_double(double arg);

int TUT_function8_int();

double TUT_function8_double();

void TUT_function9(double arg);

void TUT_function10_0();

void TUT_function10_1(const char * name, double arg2);

void TUT_function10_1_bufferify(const char * name, int Lname, double arg2);

int TUT_overload1_num(int num);

int TUT_overload1_num_offset(int num, int offset);

int TUT_overload1_num_offset_stride(int num, int offset, int stride);

int TUT_overload1_3(double type, int num);

int TUT_overload1_4(double type, int num, int offset);

int TUT_overload1_5(double type, int num, int offset, int stride);

int TUT_typefunc(int arg);

int TUT_enumfunc(int arg);

const char * TUT_last_function_called();

#ifdef __cplusplus
}
#endif

#endif  // WRAPTUTORIAL_H

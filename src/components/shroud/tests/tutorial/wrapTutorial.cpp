// wrapTutorial.cpp
// This is generated code, do not edit
// wrapTutorial.cpp
#include "wrapTutorial.h"
#include "tutorial.hpp"

extern "C" {
namespace tutorial {

// void Function1()
// function_index=3
void TUT_function1()
{
// splicer begin function.function1
Function1();
return;
// splicer end function.function1
}

// double Function2(double arg1+intent(in)+value, int arg2+intent(in)+value)
// function_index=4
double TUT_function2(double arg1, int arg2)
{
// splicer begin function.function2
double rv = Function2(arg1, arg2);
return rv;
// splicer end function.function2
}

// void Sum(int len+intent(in)+value, int * values+intent(in)+dimension((*)), int * result+intent(out))
// function_index=5
void TUT_sum(int len, int * values, int * result)
{
// splicer begin function.sum
Sum(len, values, result);
return;
// splicer end function.sum
}

// bool Function3(bool arg+intent(in)+value)
// function_index=6
bool TUT_function3(bool arg)
{
// splicer begin function.function3
bool rv = Function3(arg);
return rv;
// splicer end function.function3
}

// const std::string & Function4a(const std::string & arg1+intent(in), const std::string & arg2+intent(in))+pure
// function_index=7
const char * TUT_function4a(const char * arg1, const char * arg2)
{
// splicer begin function.function4a
const std::string & rv = Function4a(arg1, arg2);
return rv.c_str();
// splicer end function.function4a
}

// const std::string & Function4a(const std::string & arg1+intent(in)+len_trim, const std::string & arg2+intent(in)+len_trim)+pure
// function_index=24
const char * TUT_function4a_bufferify(const char * arg1, int Larg1, const char * arg2, int Larg2)
{
// splicer begin function.function4a_bufferify
const std::string & rv = Function4a(std::string(arg1, Larg1), std::string(arg2, Larg2));
return rv.c_str();
// splicer end function.function4a_bufferify
}

// const std::string & Function4b(const std::string & arg1+intent(in), const std::string & arg2+intent(in))
// function_index=8
const char * TUT_function4b(const char * arg1, const char * arg2)
{
// splicer begin function.function4b
const std::string & rv = Function4b(arg1, arg2);
return rv.c_str();
// splicer end function.function4b
}

// const std::string & Function4b(const std::string & arg1+intent(in)+len_trim, const std::string & arg2+intent(in)+len_trim)
// function_index=25
const char * TUT_function4b_bufferify(const char * arg1, int Larg1, const char * arg2, int Larg2)
{
// splicer begin function.function4b_bufferify
const std::string & rv = Function4b(std::string(arg1, Larg1), std::string(arg2, Larg2));
return rv.c_str();
// splicer end function.function4b_bufferify
}

// double Function5(double arg1+default(3.13)+intent(in)+value, int arg2+default(5)+intent(in)+value)
// function_index=9
double TUT_function5(double arg1, int arg2)
{
// splicer begin function.function5
double rv = Function5(arg1, arg2);
return rv;
// splicer end function.function5
}

// void Function6(const std::string & name+intent(in))
// function_index=10
void TUT_function6_from_name(const char * name)
{
// splicer begin function.function6_from_name
Function6(name);
return;
// splicer end function.function6_from_name
}

// void Function6(const std::string & name+intent(in)+len_trim)
// function_index=26
void TUT_function6_from_name_bufferify(const char * name, int Lname)
{
// splicer begin function.function6_from_name_bufferify
Function6(std::string(name, Lname));
return;
// splicer end function.function6_from_name_bufferify
}

// void Function6(int indx+intent(in)+value)
// function_index=11
void TUT_function6_from_index(int indx)
{
// splicer begin function.function6_from_index
Function6(indx);
return;
// splicer end function.function6_from_index
}

// void Function7(int arg+intent(in)+value)
// function_index=20
void TUT_function7_int(int arg)
{
// splicer begin function.function7_int
Function7<int>(arg);
return;
// splicer end function.function7_int
}

// void Function7(double arg+intent(in)+value)
// function_index=21
void TUT_function7_double(double arg)
{
// splicer begin function.function7_double
Function7<double>(arg);
return;
// splicer end function.function7_double
}

// int Function8()
// function_index=22
int TUT_function8_int()
{
// splicer begin function.function8_int
int rv = Function8<int>();
return rv;
// splicer end function.function8_int
}

// double Function8()
// function_index=23
double TUT_function8_double()
{
// splicer begin function.function8_double
double rv = Function8<double>();
return rv;
// splicer end function.function8_double
}

// void Function9(double arg+intent(in)+value)
// function_index=14
void TUT_function9(double arg)
{
// splicer begin function.function9
Function9(arg);
return;
// splicer end function.function9
}

// void Function10()
// function_index=15
void TUT_function10_0()
{
// splicer begin function.function10_0
Function10();
return;
// splicer end function.function10_0
}

// void Function10(const std::string & name+intent(in), double arg2+intent(in)+value)
// function_index=16
void TUT_function10_1(const char * name, double arg2)
{
// splicer begin function.function10_1
Function10(name, arg2);
return;
// splicer end function.function10_1
}

// void Function10(const std::string & name+intent(in)+len_trim, double arg2+intent(in)+value)
// function_index=27
void TUT_function10_1_bufferify(const char * name, int Lname, double arg2)
{
// splicer begin function.function10_1_bufferify
Function10(std::string(name, Lname), arg2);
return;
// splicer end function.function10_1_bufferify
}

// int overload1(int num+intent(in)+value, int offset+default(0)+intent(in)+value, int stride+default(1)+intent(in)+value)
// function_index=17
int TUT_overload1_0(int num, int offset, int stride)
{
// splicer begin function.overload1_0
int rv = overload1(num, offset, stride);
return rv;
// splicer end function.overload1_0
}

// int overload1(double type+intent(in)+value, int num+intent(in)+value, int offset+default(0)+intent(in)+value, int stride+default(1)+intent(in)+value)
// function_index=18
int TUT_overload1_1(double type, int num, int offset, int stride)
{
// splicer begin function.overload1_1
int rv = overload1(type, num, offset, stride);
return rv;
// splicer end function.overload1_1
}

// const std::string & LastFunctionCalled()+pure
// function_index=19
const char * TUT_last_function_called()
{
// splicer begin function.last_function_called
const std::string & rv = LastFunctionCalled();
return rv.c_str();
// splicer end function.last_function_called
}

// splicer begin additional_functions
// splicer end additional_functions

}  // namespace tutorial
}  // extern "C"

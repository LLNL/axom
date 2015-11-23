// wrapTutorial.cpp
// This is generated code, do not edit
// wrapTutorial.cpp
#include "wrapTutorial.h"
#include "tutorial.hpp"

extern "C" {
namespace tutorial {

void TUT_function1()
{
// splicer begin function.function1
Function1();
return;
// splicer end function.function1
}

double TUT_function2(double arg1, int arg2)
{
// splicer begin function.function2
double rv = Function2(arg1, arg2);
return rv;
// splicer end function.function2
}

void TUT_sum(int len, int * values, int * result)
{
// splicer begin function.sum
Sum(len, values, result);
return;
// splicer end function.sum
}

bool TUT_function3(bool arg)
{
// splicer begin function.function3
bool rv = Function3(arg);
return rv;
// splicer end function.function3
}

const char * TUT_function4a(const char * arg1, const char * arg2)
{
// splicer begin function.function4a
const std::string & rv = Function4a(arg1, arg2);
return rv.c_str();
// splicer end function.function4a
}

const char * TUT_function4a_bufferify(const char * arg1, int Larg1, const char * arg2, int Larg2)
{
// splicer begin function.function4a_bufferify
const std::string & rv = Function4a(std::string(arg1, Larg1), std::string(arg2, Larg2));
return rv.c_str();
// splicer end function.function4a_bufferify
}

const char * TUT_function4b(const char * arg1, const char * arg2)
{
// splicer begin function.function4b
const std::string & rv = Function4b(arg1, arg2);
return rv.c_str();
// splicer end function.function4b
}

const char * TUT_function4b_bufferify(const char * arg1, int Larg1, const char * arg2, int Larg2)
{
// splicer begin function.function4b_bufferify
const std::string & rv = Function4b(std::string(arg1, Larg1), std::string(arg2, Larg2));
return rv.c_str();
// splicer end function.function4b_bufferify
}

double TUT_function5(double arg1, int arg2)
{
// splicer begin function.function5
double rv = Function5(arg1, arg2);
return rv;
// splicer end function.function5
}

void TUT_function6_from_name(const char * name)
{
// splicer begin function.function6_from_name
Function6(name);
return;
// splicer end function.function6_from_name
}

void TUT_function6_from_name_bufferify(const char * name, int Lname)
{
// splicer begin function.function6_from_name_bufferify
Function6(std::string(name, Lname));
return;
// splicer end function.function6_from_name_bufferify
}

void TUT_function6_from_index(int indx)
{
// splicer begin function.function6_from_index
Function6(indx);
return;
// splicer end function.function6_from_index
}

void TUT_function7_int(int arg)
{
// splicer begin function.function7_int
Function7<int>(arg);
return;
// splicer end function.function7_int
}

void TUT_function7_double(double arg)
{
// splicer begin function.function7_double
Function7<double>(arg);
return;
// splicer end function.function7_double
}

int TUT_function8_int()
{
// splicer begin function.function8_int
int rv = Function8<int>();
return rv;
// splicer end function.function8_int
}

double TUT_function8_double()
{
// splicer begin function.function8_double
double rv = Function8<double>();
return rv;
// splicer end function.function8_double
}

void TUT_function9(double arg)
{
// splicer begin function.function9
Function9(arg);
return;
// splicer end function.function9
}

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

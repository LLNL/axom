// wrapTutorial.cpp
// This is generated code, do not edit
// wrapTutorial.cpp
#define EXAMPLE_WRAPPER_IMPL
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

// splicer begin additional_functions
// splicer end additional_functions

}  // namespace tutorial
}  // extern "C"

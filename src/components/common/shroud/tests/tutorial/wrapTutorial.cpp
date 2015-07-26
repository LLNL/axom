// wrapTutorial.cpp
// This is generated code, do not edit
// wrapTutorial.cpp
#define EXAMPLE_WRAPPER_IMPL
#include "wrapTutorial.h"
#include "tutorial.hpp"

extern "C" {

void TUT_function1()
{
// splicer begin function.function1
Function1();
return;
// splicer end function.function1
}

// splicer begin additional_functions
// splicer end additional_functions

}  // extern "C"

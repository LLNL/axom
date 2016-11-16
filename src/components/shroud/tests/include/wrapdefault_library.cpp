// wrapdefault_library.cpp
// This is generated code, do not edit
// wrapdefault_library.cpp
#include "wrapdefault_library.h"
#include "global_header.hpp"

extern "C" {
namespace one {
namespace two {

void DEF_function1()
{
// splicer begin function.function1
function1();
return;
// splicer end function.function1
}

// splicer begin additional_functions
// splicer end additional_functions

}  // namespace one
}  // namespace two
}  // extern "C"

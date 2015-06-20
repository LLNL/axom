// wrapUserLibrary.cpp
// This is generated code, do not edit
// blah blah
// yada yada
//
// wrapUserLibrary.cpp
#define EXAMPLE_WRAPPER_IMPL
#include "wrapUserLibrary.h"

extern "C" {
namespace example {
namespace nested {

void AA_local_function1()
{
// splicer begin function.local_function1
local_function1();
return;
// splicer end function.local_function1
}

bool AA_is_name_valid(const char * name)
{
// splicer begin function.isNameValid
return name != NULL;
// splicer end function.isNameValid
}

}  // namespace example
}  // namespace nested
}  // extern "C"

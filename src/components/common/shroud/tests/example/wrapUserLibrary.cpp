// wrapUserLibrary.cpp
// This is generated code, do not edit
// blah blah
// yada yada
//
// wrapUserLibrary.cpp
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
// splicer begin function.is_name_valid
return name != NULL;
// splicer end function.is_name_valid
}

bool AA_is_name_valid_bufferify(const char * name, int Lname)
{
// splicer begin function.is_name_valid_bufferify
return name != NULL;
// splicer end function.is_name_valid_bufferify
}

// splicer begin additional_functions
// splicer end additional_functions

}  // namespace example
}  // namespace nested
}  // extern "C"

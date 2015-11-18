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

void AA_test_names_empty(const char * name)
{
// splicer begin function.test_names_empty
test_names(name);
return;
// splicer end function.test_names_empty
}

void AA_test_names_bufferify(const char * name, int Lname)
{
// splicer begin function.test_names_bufferify
test_names(std::string(name, Lname));
return;
// splicer end function.test_names_bufferify
}

void AA_test_names_flag(const char * name, int flag)
{
// splicer begin function.test_names_flag
test_names(name, flag);
return;
// splicer end function.test_names_flag
}

void AA_test_names_bufferify(const char * name, int Lname, int flag)
{
// splicer begin function.test_names_bufferify
test_names(std::string(name, Lname), flag);
return;
// splicer end function.test_names_bufferify
}

// splicer begin additional_functions
// splicer end additional_functions

}  // namespace example
}  // namespace nested
}  // extern "C"

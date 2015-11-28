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

// void local_function1
// function_index=33
void AA_local_function1()
{
// splicer begin function.local_function1
local_function1();
return;
// splicer end function.local_function1
}

// bool isNameValid(const std::string& name)
// function_index=34
bool AA_is_name_valid(const char * name)
{
// splicer begin function.is_name_valid
return name != NULL;
// splicer end function.is_name_valid
}

// bool isNameValid(const std::string& name)
// function_index=38
bool AA_is_name_valid_bufferify(const char * name, int Lname)
{
// splicer begin function.is_name_valid_bufferify
return name != NULL;
// splicer end function.is_name_valid_bufferify
}

// void test_names(const std::string &name)
// function_index=35
void AA_test_names(const char * name)
{
// splicer begin function.test_names
test_names(name);
return;
// splicer end function.test_names
}

// void test_names(const std::string &name)
// function_index=39
void AA_test_names_bufferify(const char * name, int Lname)
{
// splicer begin function.test_names_bufferify
test_names(std::string(name, Lname));
return;
// splicer end function.test_names_bufferify
}

// void test_names(const std::string &name, int flag)
// function_index=36
void AA_test_names_flag(const char * name, int flag)
{
// splicer begin function.test_names_flag
test_names(name, flag);
return;
// splicer end function.test_names_flag
}

// void test_names(const std::string &name, int flag)
// function_index=40
void AA_test_names_flag_bufferify(const char * name, int Lname, int flag)
{
// splicer begin function.test_names_flag_bufferify
test_names(std::string(name, Lname), flag);
return;
// splicer end function.test_names_flag_bufferify
}

// void testoptional(int i = 1, long j=2)
// function_index=37
void AA_testoptional(int i, long j)
{
// splicer begin function.testoptional
testoptional(i, j);
return;
// splicer end function.testoptional
}

// splicer begin additional_functions
// splicer end additional_functions

}  // namespace example
}  // namespace nested
}  // extern "C"

// wrapUserLibrary.cpp
// This is generated code, do not edit
// blah blah
// yada yada
//
// wrapUserLibrary.cpp
#include "wrapUserLibrary.h"
#include <string>

extern "C" {
namespace example {
namespace nested {

// void local_function1()
// function_index=41
void AA_local_function1()
{
// splicer begin function.local_function1
local_function1();
return;
// splicer end function.local_function1
}

// bool isNameValid(const std::string & name+intent(in))
// function_index=42
bool AA_is_name_valid(const char * name)
{
// splicer begin function.is_name_valid
return name != NULL;
// splicer end function.is_name_valid
}

// bool isNameValid(const std::string & name+intent(in)+len_trim(Lname))
// function_index=50
bool AA_is_name_valid_bufferify(const char * name, int Lname)
{
// splicer begin function.is_name_valid_bufferify
return name != NULL;
// splicer end function.is_name_valid_bufferify
}

// bool isInitialized()
// function_index=43
bool AA_is_initialized()
{
// splicer begin function.is_initialized
bool rv = isInitialized();
return rv;
// splicer end function.is_initialized
}

// void test_names(const std::string & name+intent(in))
// function_index=44
void AA_test_names(const char * name)
{
// splicer begin function.test_names
std::string SH_name(name);
test_names(SH_name);
return;
// splicer end function.test_names
}

// void test_names(const std::string & name+intent(in)+len_trim(Lname))
// function_index=51
void AA_test_names_bufferify(const char * name, int Lname)
{
// splicer begin function.test_names_bufferify
std::string SH_name(name, Lname);
test_names(SH_name);
return;
// splicer end function.test_names_bufferify
}

// void test_names(const std::string & name+intent(in), int flag+intent(in)+value)
// function_index=45
void AA_test_names_flag(const char * name, int flag)
{
// splicer begin function.test_names_flag
std::string SH_name(name);
test_names(SH_name, flag);
return;
// splicer end function.test_names_flag
}

// void test_names(const std::string & name+intent(in)+len_trim(Lname), int flag+intent(in)+value)
// function_index=52
void AA_test_names_flag_bufferify(const char * name, int Lname, int flag)
{
// splicer begin function.test_names_flag_bufferify
std::string SH_name(name, Lname);
test_names(SH_name, flag);
return;
// splicer end function.test_names_flag_bufferify
}

// void testoptional()
// function_index=48
void AA_testoptional_0()
{
// splicer begin function.testoptional_0
testoptional();
return;
// splicer end function.testoptional_0
}

// void testoptional(int i+default(1)+intent(in)+value)
// function_index=49
void AA_testoptional_1(int i)
{
// splicer begin function.testoptional_1
testoptional(i);
return;
// splicer end function.testoptional_1
}

// void testoptional(int i+default(1)+intent(in)+value, long j+default(2)+intent(in)+value)
// function_index=46
void AA_testoptional_2(int i, long j)
{
// splicer begin function.testoptional_2
testoptional(i, j);
return;
// splicer end function.testoptional_2
}

// void testmpi(MPI_Comm comm+intent(in)+value)
// function_index=47
void AA_testmpi(MPI_Fint comm)
{
// splicer begin function.testmpi
testmpi(MPI_Comm_f2c(comm));
return;
// splicer end function.testmpi
}

// splicer begin additional_functions
// splicer end additional_functions

}  // namespace example
}  // namespace nested
}  // extern "C"

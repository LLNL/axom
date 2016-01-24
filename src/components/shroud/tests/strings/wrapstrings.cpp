// wrapstrings.cpp
// This is generated code, do not edit
// wrapstrings.cpp
#include "wrapstrings.h"
#include <string>
#include "shroudrt.hpp"
#include "strings.hpp"

extern "C" {

// const string & getString1()+pure
// function_index=0
/**
 * \brief return a string as character(*)
 *
 */
const char * STR_get_string1()
{
// splicer begin function.get_string1
const std::string & rv = getString1();
return rv.c_str();
// splicer end function.get_string1
}

// void getString1(string_result_as_arg & SH_F_rv+intent(out)+len(LSH_F_rv))+pure
// function_index=5
/**
 * \brief return a string as character(*)
 *
 */
void STR_get_string1_bufferify(char * SH_F_rv, int LSH_F_rv)
{
// splicer begin function.get_string1_bufferify
const std::string & rv = getString1();
asctoolkit::shroud::FccCopy(SH_F_rv, LSH_F_rv, rv.c_str());
return;
// splicer end function.get_string1_bufferify
}

// const string & getString2()
// function_index=1
/**
 * \brief return string with fixed size (len=30)
 *
 */
const char * STR_get_string2()
{
// splicer begin function.get_string2
const std::string & rv = getString2();
// check for error
if (rv.empty()) {
    return NULL;
}

return rv.c_str();
// splicer end function.get_string2
}

// void getString2(string_result_as_arg & SH_F_rv+intent(out)+len(LSH_F_rv))
// function_index=7
/**
 * \brief return string with fixed size (len=30)
 *
 */
void STR_get_string2_bufferify(char * SH_F_rv, int LSH_F_rv)
{
// splicer begin function.get_string2_bufferify
const std::string & rv = getString2();
// check for error
if (rv.empty()) {
    std::memset(SH_F_rv, ' ', LSH_F_rv);
    return;
}

asctoolkit::shroud::FccCopy(SH_F_rv, LSH_F_rv, rv.c_str());
return;
// splicer end function.get_string2_bufferify
}

// const string & getString3()
// function_index=2
/**
 * \brief return a string as argument
 *
 */
const char * STR_get_string3()
{
// splicer begin function.get_string3
const std::string & rv = getString3();
// check for error
if (rv.empty()) {
    return NULL;
}

return rv.c_str();
// splicer end function.get_string3
}

// void getString3(string_result_as_arg & output+intent(out)+len(Loutput))
// function_index=8
/**
 * \brief return a string as argument
 *
 */
void STR_get_string3_bufferify(char * output, int Loutput)
{
// splicer begin function.get_string3_bufferify
const std::string & rv = getString3();
// check for error
if (rv.empty()) {
    std::memset(output, ' ', Loutput);
    return;
}

asctoolkit::shroud::FccCopy(output, Loutput, rv.c_str());
return;
// splicer end function.get_string3_bufferify
}

// void acceptStringConstReference(const std::string & arg1+intent(in))
// function_index=3
/**
 * \brief Accept a const string reference
 *
 * Save contents of arg1.
 * arg1 is assumed to be intent(IN) since it is const
 * Will copy in.
 */
void STR_accept_string_const_reference(const char * arg1)
{
// splicer begin function.accept_string_const_reference
std::string SH_arg1(arg1);
acceptStringConstReference(SH_arg1);
return;
// splicer end function.accept_string_const_reference
}

// void acceptStringConstReference(const std::string & arg1+intent(in)+len_trim(Larg1))
// function_index=10
/**
 * \brief Accept a const string reference
 *
 * Save contents of arg1.
 * arg1 is assumed to be intent(IN) since it is const
 * Will copy in.
 */
void STR_accept_string_const_reference_bufferify(const char * arg1, int Larg1)
{
// splicer begin function.accept_string_const_reference_bufferify
std::string SH_arg1(arg1, Larg1);
acceptStringConstReference(SH_arg1);
return;
// splicer end function.accept_string_const_reference_bufferify
}

// void acceptStringReference(std::string & arg1+intent(inout))
// function_index=4
/**
 * \brief Accept a string reference
 *
 * Append "dog" to the end of arg1.
 * arg1 is assumed to be intent(INOUT)
 * Must copy in and copy out.
 */
void STR_accept_string_reference(char * arg1)
{
// splicer begin function.accept_string_reference
std::string SH_arg1(arg1);
acceptStringReference(SH_arg1);
asctoolkit::shroud::FccCopy(arg1, 0, SH_arg1.c_str());
return;
// splicer end function.accept_string_reference
}

// void acceptStringReference(std::string & arg1+intent(inout)+len(Narg1)+len_trim(Larg1))
// function_index=11
/**
 * \brief Accept a string reference
 *
 * Append "dog" to the end of arg1.
 * arg1 is assumed to be intent(INOUT)
 * Must copy in and copy out.
 */
void STR_accept_string_reference_bufferify(char * arg1, int Larg1, int Narg1)
{
// splicer begin function.accept_string_reference_bufferify
std::string SH_arg1(arg1, Larg1);
acceptStringReference(SH_arg1);
asctoolkit::shroud::FccCopy(arg1, Narg1, SH_arg1.c_str());
return;
// splicer end function.accept_string_reference_bufferify
}

// splicer begin additional_functions
// splicer end additional_functions

}  // extern "C"

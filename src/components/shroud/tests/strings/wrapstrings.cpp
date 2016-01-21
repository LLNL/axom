// wrapstrings.cpp
// This is generated code, do not edit
// wrapstrings.cpp
#include "wrapstrings.h"
#include <string>
#include "shroudrt.hpp"
#include "strings.hpp"

extern "C" {

// const string & getName()+pure
// function_index=0
/**
 * \brief return a string as argument
 *
 */
const char * STR_get_name()
{
// splicer begin function.get_name
const std::string & rv = getName();
// check for error
if (rv.empty()) {
    return NULL;
}

return rv.c_str();
// splicer end function.get_name
}

// void getName(string_result_as_arg & output+intent(out)+len(Loutput))+pure
// function_index=3
/**
 * \brief return a string as argument
 *
 */
void STR_get_name_bufferify(char * output, int Loutput)
{
// splicer begin function.get_name_bufferify
const std::string & rv = getName();
// check for error
if (rv.empty()) {
    std::memset(output, ' ', Loutput);
    return;
}

asctoolkit::shroud::FccCopy(output, Loutput, rv.c_str());
return;
// splicer end function.get_name_bufferify
}

// void acceptStringConstReference(const std::string & arg1+intent(in))
// function_index=1
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
// function_index=5
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

// void acceptStringReference(std::string & arg1+intent(in))
// function_index=2
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
return;
// splicer end function.accept_string_reference
}

// void acceptStringReference(std::string & arg1+len_trim(Larg1)+intent(in))
// function_index=7
/**
 * \brief Accept a string reference
 *
 * Append "dog" to the end of arg1.
 * arg1 is assumed to be intent(INOUT)
 * Must copy in and copy out.
 */
void STR_accept_string_reference_bufferify(char * arg1, int Larg1)
{
// splicer begin function.accept_string_reference_bufferify
std::string SH_arg1(arg1, Larg1);
acceptStringReference(SH_arg1);
return;
// splicer end function.accept_string_reference_bufferify
}

// splicer begin additional_functions
// splicer end additional_functions

}  // extern "C"

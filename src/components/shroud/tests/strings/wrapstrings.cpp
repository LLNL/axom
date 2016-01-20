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
// function_index=2
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
 */
void STR_accept_string_const_reference(const char * arg1)
{
// splicer begin function.accept_string_const_reference
acceptStringConstReference(arg1);
return;
// splicer end function.accept_string_const_reference
}

// void acceptStringConstReference(const std::string & arg1+intent(in)+len_trim(Larg1))
// function_index=4
/**
 * \brief Accept a const string reference
 *
 */
void STR_accept_string_const_reference_bufferify(const char * arg1, int Larg1)
{
// splicer begin function.accept_string_const_reference_bufferify
acceptStringConstReference(std::string(arg1, Larg1));
return;
// splicer end function.accept_string_const_reference_bufferify
}

// splicer begin additional_functions
// splicer end additional_functions

}  // extern "C"

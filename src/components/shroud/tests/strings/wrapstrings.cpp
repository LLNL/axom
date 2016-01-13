// wrapstrings.cpp
// This is generated code, do not edit
// wrapstrings.cpp
#include "wrapstrings.h"

extern "C" {

// const string & getName() const
// function_index=0
/**
 * \brief return a string as argument
 *
 */
const char * STR_get_name()
{
// splicer begin function.get_name
const std::string & rv = getName();
return rv.c_str();
// splicer end function.get_name
}

// splicer begin additional_functions
// splicer end additional_functions

}  // extern "C"

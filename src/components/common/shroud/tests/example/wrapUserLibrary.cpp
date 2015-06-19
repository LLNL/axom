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

// splicer begin AA_local_function1
local_function1();
return;
// splicer end AA_local_function1
}

bool AA_is_name_valid(const char * name)
{

// splicer begin AA_is_name_valid
return name != NULL;
// splicer end AA_is_name_valid
}

}  // namespace example
}  // namespace nested
}  // extern "C"

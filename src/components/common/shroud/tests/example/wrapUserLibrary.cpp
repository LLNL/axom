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
selfobj->local_function1();
return;
// splicer end AA_local_function1
}

}  // namespace example
}  // namespace nested
}  // extern "C"

// wrapClass2.cpp
// This is generated code, do not edit
// wrapClass2.cpp
#include "wrapClass2.h"
#include "global_header.hpp"

extern "C" {
namespace one {
namespace two {

void DEF_class2_method1(DEF_class2 * self)
{
Class2 *selfobj = static_cast<Class2 *>(static_cast<void *>(self));
// splicer begin class.Class2.method.method1
selfobj->method1();
return;
// splicer end class.Class2.method.method1
}

// splicer begin class.Class2.additional_functions
// splicer end class.Class2.additional_functions

}  // namespace one
}  // namespace two
}  // extern "C"

// wrapClass1.cpp
// This is generated code, do not edit
// wrapClass1.cpp
#define EXAMPLE_WRAPPER_IMPL
#include "wrapClass1.h"
#include "tutorial.hpp"

extern "C" {
namespace tutorial {

TUT_class1 * TUT_class1_new()
{
Class1 *selfobj = new Class1();
// splicer begin class.Class1.method.new
return (TUT_class1 *) selfobj;
// splicer end class.Class1.method.new
}

void TUT_class1_method1(TUT_class1 * self)
{
Class1 *selfobj = static_cast<Class1 *>(self);
// splicer begin class.Class1.method.method1
selfobj->Method1();
return;
// splicer end class.Class1.method.method1
}

// splicer begin class.Class1.additional_functions
// splicer end class.Class1.additional_functions

}  // namespace tutorial
}  // extern "C"

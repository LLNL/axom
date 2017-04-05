// wrapClass1.cpp
// This is generated code, do not edit
// wrapClass1.cpp
#include "wrapClass1.h"
#include "tutorial.hpp"

extern "C" {
namespace tutorial {

// Class1 * new()+constructor
// function_index=0
TUT_class1 * TUT_class1_new()
{
// splicer begin class.Class1.method.new
    Class1 * SH_rv = new Class1();
    return static_cast<TUT_class1 *>(static_cast<void *>(SH_rv));
// splicer end class.Class1.method.new
}

// void delete()+destructor
// function_index=1
void TUT_class1_delete(TUT_class1 * self)
{
Class1 *SH_this = static_cast<Class1 *>(static_cast<void *>(self));
// splicer begin class.Class1.method.delete
    delete SH_this;
// splicer end class.Class1.method.delete
}

// void Method1()
// function_index=2
void TUT_class1_method1(TUT_class1 * self)
{
Class1 *SH_this = static_cast<Class1 *>(static_cast<void *>(self));
// splicer begin class.Class1.method.method1
    SH_this->Method1();
    return;
// splicer end class.Class1.method.method1
}

// splicer begin class.Class1.additional_functions
// splicer end class.Class1.additional_functions

}  // namespace tutorial
}  // extern "C"

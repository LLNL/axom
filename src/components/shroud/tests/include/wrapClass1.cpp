// wrapClass1.cpp
// This is generated code, do not edit
// wrapClass1.cpp
#include "wrapClass1.h"
#include "class_header.hpp"
#include "type_header.hpp"

namespace three {


extern "C" {


void DEF_class1_method1(DEF_class1 * self, int arg1)
{
    Class1 *SH_this = static_cast<Class1 *>(static_cast<void *>(self));
    SH_this->method1(arg1);
    return;
}

}  // extern "C"

}  // namespace three

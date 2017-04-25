// foo.cpp
// This is generated code, do not edit
// foo.cpp
#include "foo.h"

// splicer begin class.Names.CXX_definitions
// splicer end class.Names.CXX_definitions

extern "C" {

// splicer begin class.Names.C_definitions
// splicer end class.Names.C_definitions

// void method1()
// function_index=0
void XXX_TES_names_method1(TES_names * self)
{
// splicer begin class.Names.method.method1
    Names *SH_this = static_cast<Names *>(static_cast<void *>(self));
    SH_this->method1();
    return;
// splicer end class.Names.method.method1
}

// void method2()
// function_index=1
void XXX_TES_names_method2(TES_names * self)
{
// splicer begin class.Names.method.method2
    Names *SH_this = static_cast<Names *>(static_cast<void *>(self));
    SH_this->method2();
    return;
// splicer end class.Names.method.method2
}

}  // extern "C"

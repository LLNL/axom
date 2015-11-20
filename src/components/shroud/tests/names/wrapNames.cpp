// wrapNames.cpp
// This is generated code, do not edit
// wrapNames.cpp
#include "wrapNames.h"

extern "C" {

void DEF_names_method1(DEF_names * self)
{
Names *selfobj = static_cast<Names *>(static_cast<void *>(self));
// splicer begin class.Names.method.method1
selfobj->method1();
return;
// splicer end class.Names.method.method1
}

void DEF_names_method2(DEF_names * self)
{
Names *selfobj = static_cast<Names *>(static_cast<void *>(self));
// splicer begin class.Names.method.method2
selfobj->method2();
return;
// splicer end class.Names.method.method2
}

// splicer begin class.Names.additional_functions
// splicer end class.Names.additional_functions

}  // extern "C"

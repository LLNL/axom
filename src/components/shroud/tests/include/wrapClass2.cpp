// wrapClass2.cpp
// This is generated code, do not edit
// wrapClass2.cpp
#include "wrapClass2.h"
#include "global_header.hpp"

extern "C" {

void DEF_class2_method1(DEF_class2 * self, MPI_Fint comm)
{
Class2 *selfobj = static_cast<Class2 *>(static_cast<void *>(self));
// splicer begin class.Class2.method.method1
    selfobj->method1(MPI_Comm_f2c(comm));
    return;
// splicer end class.Class2.method.method1
}

void DEF_class2_method2(DEF_class2 * self, DEF_class1 * c2)
{
Class2 *selfobj = static_cast<Class2 *>(static_cast<void *>(self));
// splicer begin class.Class2.method.method2
    selfobj->method2(static_cast<Class1 *>(static_cast<void *>(c2)));
    return;
// splicer end class.Class2.method.method2
}

// splicer begin class.Class2.additional_functions
// splicer end class.Class2.additional_functions

}  // extern "C"

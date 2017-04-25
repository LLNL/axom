// wrapClass2.cpp
// This is generated code, do not edit
// wrapClass2.cpp
#include "wrapClass2.h"
#include "global_header.hpp"

// splicer begin class.Class2.CXX_definitions
// splicer end class.Class2.CXX_definitions

extern "C" {

// splicer begin class.Class2.C_definitions
// splicer end class.Class2.C_definitions

void DEF_class2_method1(DEF_class2 * self, MPI_Fint comm)
{
// splicer begin class.Class2.method.method1
    Class2 *SH_this = static_cast<Class2 *>(static_cast<void *>(self));
    SH_this->method1(MPI_Comm_f2c(comm));
    return;
// splicer end class.Class2.method.method1
}

void DEF_class2_method2(DEF_class2 * self, DEF_class1 * c2)
{
// splicer begin class.Class2.method.method2
    Class2 *SH_this = static_cast<Class2 *>(static_cast<void *>(self));
    SH_this->method2(static_cast<Class1 *>(static_cast<void *>(c2)));
    return;
// splicer end class.Class2.method.method2
}

}  // extern "C"

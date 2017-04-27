// wrapClass2.cpp
// This is generated code, do not edit
// wrapClass2.cpp
#include "wrapClass2.h"
#include "global_header.hpp"


extern "C" {


void DEF_class2_method1(DEF_class2 * self, MPI_Fint comm)
{
    Class2 *SH_this = static_cast<Class2 *>(static_cast<void *>(self));
    SH_this->method1(MPI_Comm_f2c(comm));
    return;
}

void DEF_class2_method2(DEF_class2 * self, DEF_class1 * c2)
{
    Class2 *SH_this = static_cast<Class2 *>(static_cast<void *>(self));
    SH_this->method2(static_cast<Class1 *>(static_cast<void *>(c2)));
    return;
}

}  // extern "C"

// wrapSidre.cpp
// This is generated code, do not edit
//
// Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
//
// Produced at the Lawrence Livermore National Laboratory
//
// LLNL-CODE-741217
//
// All rights reserved.
//
// This file is part of Axom.
//
// For details about use and distribution, please read axom/LICENSE.
//
#include "wrapSidre.h"
#include <stdlib.h>
#include <string>
#include "axom/sidre/core/DataStore.hpp"
#include "typesSidre.h"

// splicer begin CXX_definitions
// splicer end CXX_definitions

extern "C" {

// splicer begin C_definitions
// equivalent to C_LOC
// called from Fortran
// https://gcc.gnu.org/bugzilla/show_bug.cgi?id=53945
// Work around a problem with gfortran 4.7 where C_LOC does not work
// with assumed shape array.  Passing the first element of the
// array to a function without an interface will force the compiler
// to use f77 semantics and pass the address of the data, essentially
// the same as C_LOC.
// XXX Pass the first element, not the entire array, to avoid getting
// XXX a copy of the array.
//
// The result must be an argument because some compilers (Intel)
// cannot return type(C_PTR)
void sidre_c_loc(void* addr, void** out)
{
  *out = addr;
}
void sidre_c_loc_(void* addr, void** out)
{
  *out = addr;
}

// splicer end C_definitions

bool SIDRE_name_is_valid(const char * name)
{
// splicer begin function.name_is_valid
    return name != NULL;
// splicer end function.name_is_valid
}

// Release C++ allocated memory.
void SIDRE_SHROUD_memory_destructor(SID_SHROUD_capsule_data *cap)
{
    void *ptr = cap->addr;
    switch (cap->idtor) {
    case 0:   // --none--
    {
        // Nothing to delete
        break;
    }
    case 1:   // axom::sidre::DataStore
    {
        axom::sidre::DataStore *cxx_ptr = reinterpret_cast<axom::sidre::DataStore *>(ptr);
        delete cxx_ptr;
        break;
    }
    default:
    {
        // Unexpected case in destructor
        break;
    }
    }
    cap->addr = NULL;
    cap->idtor = 0;  // avoid deleting again
}

}  // extern "C"

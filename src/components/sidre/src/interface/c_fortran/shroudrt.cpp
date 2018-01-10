// shroudrt.cpp
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


/* *INDENT-OFF* */
#ifdef __cplusplus
// Standard C++ headers
#include <cstring>
using namespace std;
extern "C" {
#else
#include <string.h>
#endif
/* *INDENT-ON* */

void shroud_FccCopy(char* a, int la, const char* s)
{
  int ls,nm;
  ls = strlen(s);
  nm = ls < la ? ls : la;
  memcpy(a,s,nm);
  if(la > nm)
  {
    memset(a+nm,' ',la-nm);
  }
}

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
void shroud_c_loc(void* addr, void** out)
{
  *out = addr;
}
void shroud_c_loc_(void* addr, void** out)
{
  *out = addr;
}


/* *INDENT-OFF* */
#ifdef __cplusplus
}  // extern "C"
#endif
/* *INDENT-ON* */

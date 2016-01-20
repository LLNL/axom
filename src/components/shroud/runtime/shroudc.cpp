/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */
//
// Support functions for shroud
//
// Standard C++ headers
#include <cstring>

// Other CS Toolkit headers
#include "common/FC.h"


namespace asctoolkit
{
namespace shroud
{

extern "C" {

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
void FC_GLOBAL(shroud_c_loc,SHROUD_C_LOC)(void * addr, void * * out)
{
  *out = addr;
}

}  // extern "C"

/*--------------------------------------------------------------------------*/

/* copy a C-string into a Fortran character variable
 * blank-fill result
 */
void FccCopy(char *a, int la, const char *s)
{
   int ls,nm;
   ls = std::strlen(s);
   nm = ls < la ? ls : la;
   memcpy(a,s,nm);
   if(la > nm) { memset(a+nm,' ',la-nm);}
}

/*--------------------------------------------------------------------------*/

} /* end namespace sidre */
} /* end namespace asctoolkit */

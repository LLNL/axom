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
// SidreAllocatable.cpp - Routines used by Fortran interface
//
#include <cstddef>
#include "common/CommonTypes.hpp"
#include "common/FC.h"

#include "sidre/SidreTypes.hpp"
#include "sidre/DataGroup.hpp"
#include "sidre/SidreAllocatable.hpp"

namespace asctoolkit
{
namespace sidre
{
class DataView;

extern "C" {

/*!
 * \brief Return DataView for a Fortran array.
 *
 * Create an external view using information derived from a Fortran array.
 * Called from Fortran.
 */
void * SIDRE_create_array_view(void * group, char * name, int lname,
                               void * addr, int type, int rank, SidreLength * extents)
{
  DataGroup * grp = static_cast<DataGroup *>(group);
#if 0
  DataView * view = grp->createExternalView(std::string(name, lname),
                                            addr,
                                            static_cast<TypeID>(type),
                                            static_cast<SidreLength>(extents[0]));
#endif
  DataView * view = grp->createExternalView(std::string(name, lname),
                                            addr,
                                            static_cast<TypeID>(type),
                                            rank, extents);
  return view;
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
void FC_GLOBAL(sidre_c_loc,SIDRE_C_LOC)(void * addr, void * * out)
{
  *out = addr;
}

}  // extern "C"


}  // namespace asctoolkit
}  // namespace sidre

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
  DataView * view = grp->createExternalView(std::string(name, lname),
                                            addr,
                                            static_cast<TypeID>(type),
                                            rank, extents);
  return view;
}


}  // extern "C"


}  // namespace asctoolkit
}  // namespace sidre

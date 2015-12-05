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
// Uses cog to insert some generated code into this file.
//
#include <cstddef>
#include "common/CommonTypes.hpp"
#include "common/FC.h"

#include "sidre/DataGroup.hpp"
#include "sidre/SidreAllocatable.hpp"

// import cog once
//[[[cog import cog;import genfsidresplicer as gen ]]]
//[[[end]]]

namespace asctoolkit
{
namespace sidre
{
class DataView;

extern "C" {
//[[[cog
//gen.print_lines(cog.outl, gen.print_atk_size_allocatable_header)
//]]]
#define SIZE_ALLOCATABLE_INT_SCALAR FC_GLOBAL(atk_size_allocatable_int_scalar,ATK_SIZE_ALLOCATABLE_INT_SCALAR)
size_t SIZE_ALLOCATABLE_INT_SCALAR(void * array);

#define SIZE_ALLOCATABLE_INT_1D FC_GLOBAL(atk_size_allocatable_int_1d,ATK_SIZE_ALLOCATABLE_INT_1D)
size_t SIZE_ALLOCATABLE_INT_1D(void * array);

#define SIZE_ALLOCATABLE_LONG_SCALAR FC_GLOBAL(atk_size_allocatable_long_scalar,ATK_SIZE_ALLOCATABLE_LONG_SCALAR)
size_t SIZE_ALLOCATABLE_LONG_SCALAR(void * array);

#define SIZE_ALLOCATABLE_LONG_1D FC_GLOBAL(atk_size_allocatable_long_1d,ATK_SIZE_ALLOCATABLE_LONG_1D)
size_t SIZE_ALLOCATABLE_LONG_1D(void * array);

#define SIZE_ALLOCATABLE_FLOAT_SCALAR FC_GLOBAL(atk_size_allocatable_float_scalar,ATK_SIZE_ALLOCATABLE_FLOAT_SCALAR)
size_t SIZE_ALLOCATABLE_FLOAT_SCALAR(void * array);

#define SIZE_ALLOCATABLE_FLOAT_1D FC_GLOBAL(atk_size_allocatable_float_1d,ATK_SIZE_ALLOCATABLE_FLOAT_1D)
size_t SIZE_ALLOCATABLE_FLOAT_1D(void * array);

#define SIZE_ALLOCATABLE_DOUBLE_SCALAR FC_GLOBAL(atk_size_allocatable_double_scalar,ATK_SIZE_ALLOCATABLE_DOUBLE_SCALAR)
size_t SIZE_ALLOCATABLE_DOUBLE_SCALAR(void * array);

#define SIZE_ALLOCATABLE_DOUBLE_1D FC_GLOBAL(atk_size_allocatable_double_1d,ATK_SIZE_ALLOCATABLE_DOUBLE_1D)
size_t SIZE_ALLOCATABLE_DOUBLE_1D(void * array);

//[[[end]]]

//[[[cog
//gen.print_lines(cog.outl, gen.print_atk_address_allocatable_header)
//]]]
#define ADDRESS_ALLOCATABLE_INT_SCALAR FC_GLOBAL(atk_address_allocatable_int_scalar,ATK_ADDRESS_ALLOCATABLE_INT_SCALAR)
void ADDRESS_ALLOCATABLE_INT_SCALAR(void * array, void * * addr);

#define ADDRESS_ALLOCATABLE_INT_1D FC_GLOBAL(atk_address_allocatable_int_1d,ATK_ADDRESS_ALLOCATABLE_INT_1D)
void ADDRESS_ALLOCATABLE_INT_1D(void * array, void * * addr);

#define ADDRESS_ALLOCATABLE_LONG_SCALAR FC_GLOBAL(atk_address_allocatable_long_scalar,ATK_ADDRESS_ALLOCATABLE_LONG_SCALAR)
void ADDRESS_ALLOCATABLE_LONG_SCALAR(void * array, void * * addr);

#define ADDRESS_ALLOCATABLE_LONG_1D FC_GLOBAL(atk_address_allocatable_long_1d,ATK_ADDRESS_ALLOCATABLE_LONG_1D)
void ADDRESS_ALLOCATABLE_LONG_1D(void * array, void * * addr);

#define ADDRESS_ALLOCATABLE_FLOAT_SCALAR FC_GLOBAL(atk_address_allocatable_float_scalar,ATK_ADDRESS_ALLOCATABLE_FLOAT_SCALAR)
void ADDRESS_ALLOCATABLE_FLOAT_SCALAR(void * array, void * * addr);

#define ADDRESS_ALLOCATABLE_FLOAT_1D FC_GLOBAL(atk_address_allocatable_float_1d,ATK_ADDRESS_ALLOCATABLE_FLOAT_1D)
void ADDRESS_ALLOCATABLE_FLOAT_1D(void * array, void * * addr);

#define ADDRESS_ALLOCATABLE_DOUBLE_SCALAR FC_GLOBAL(atk_address_allocatable_double_scalar,ATK_ADDRESS_ALLOCATABLE_DOUBLE_SCALAR)
void ADDRESS_ALLOCATABLE_DOUBLE_SCALAR(void * array, void * * addr);

#define ADDRESS_ALLOCATABLE_DOUBLE_1D FC_GLOBAL(atk_address_allocatable_double_1d,ATK_ADDRESS_ALLOCATABLE_DOUBLE_1D)
void ADDRESS_ALLOCATABLE_DOUBLE_1D(void * array, void * * addr);

//[[[end]]]

//[[[cog
//gen.print_lines(cog.outl, gen.print_atk_allocate_allocatable_header)
//]]]
#define ALLOCATE_ALLOCATABLE_INT_SCALAR FC_GLOBAL(atk_allocate_allocatable_int_scalar,ATK_ALLOCATE_ALLOCATABLE_INT_SCALAR)
void ALLOCATE_ALLOCATABLE_INT_SCALAR(void * array, long * nitems);

#define ALLOCATE_ALLOCATABLE_INT_1D FC_GLOBAL(atk_allocate_allocatable_int_1d,ATK_ALLOCATE_ALLOCATABLE_INT_1D)
void ALLOCATE_ALLOCATABLE_INT_1D(void * array, long * nitems);

#define ALLOCATE_ALLOCATABLE_LONG_SCALAR FC_GLOBAL(atk_allocate_allocatable_long_scalar,ATK_ALLOCATE_ALLOCATABLE_LONG_SCALAR)
void ALLOCATE_ALLOCATABLE_LONG_SCALAR(void * array, long * nitems);

#define ALLOCATE_ALLOCATABLE_LONG_1D FC_GLOBAL(atk_allocate_allocatable_long_1d,ATK_ALLOCATE_ALLOCATABLE_LONG_1D)
void ALLOCATE_ALLOCATABLE_LONG_1D(void * array, long * nitems);

#define ALLOCATE_ALLOCATABLE_FLOAT_SCALAR FC_GLOBAL(atk_allocate_allocatable_float_scalar,ATK_ALLOCATE_ALLOCATABLE_FLOAT_SCALAR)
void ALLOCATE_ALLOCATABLE_FLOAT_SCALAR(void * array, long * nitems);

#define ALLOCATE_ALLOCATABLE_FLOAT_1D FC_GLOBAL(atk_allocate_allocatable_float_1d,ATK_ALLOCATE_ALLOCATABLE_FLOAT_1D)
void ALLOCATE_ALLOCATABLE_FLOAT_1D(void * array, long * nitems);

#define ALLOCATE_ALLOCATABLE_DOUBLE_SCALAR FC_GLOBAL(atk_allocate_allocatable_double_scalar,ATK_ALLOCATE_ALLOCATABLE_DOUBLE_SCALAR)
void ALLOCATE_ALLOCATABLE_DOUBLE_SCALAR(void * array, long * nitems);

#define ALLOCATE_ALLOCATABLE_DOUBLE_1D FC_GLOBAL(atk_allocate_allocatable_double_1d,ATK_ALLOCATE_ALLOCATABLE_DOUBLE_1D)
void ALLOCATE_ALLOCATABLE_DOUBLE_1D(void * array, long * nitems);

//[[[end]]]

//[[[cog
//gen.print_lines(cog.outl, gen.print_atk_deallocate_allocatable_header)
//]]]
#define DEALLOCATE_ALLOCATABLE_INT_SCALAR FC_GLOBAL(atk_deallocate_allocatable_int_scalar,ATK_DEALLOCATE_ALLOCATABLE_INT_SCALAR)
void DEALLOCATE_ALLOCATABLE_INT_SCALAR(void * array);

#define DEALLOCATE_ALLOCATABLE_INT_1D FC_GLOBAL(atk_deallocate_allocatable_int_1d,ATK_DEALLOCATE_ALLOCATABLE_INT_1D)
void DEALLOCATE_ALLOCATABLE_INT_1D(void * array);

#define DEALLOCATE_ALLOCATABLE_LONG_SCALAR FC_GLOBAL(atk_deallocate_allocatable_long_scalar,ATK_DEALLOCATE_ALLOCATABLE_LONG_SCALAR)
void DEALLOCATE_ALLOCATABLE_LONG_SCALAR(void * array);

#define DEALLOCATE_ALLOCATABLE_LONG_1D FC_GLOBAL(atk_deallocate_allocatable_long_1d,ATK_DEALLOCATE_ALLOCATABLE_LONG_1D)
void DEALLOCATE_ALLOCATABLE_LONG_1D(void * array);

#define DEALLOCATE_ALLOCATABLE_FLOAT_SCALAR FC_GLOBAL(atk_deallocate_allocatable_float_scalar,ATK_DEALLOCATE_ALLOCATABLE_FLOAT_SCALAR)
void DEALLOCATE_ALLOCATABLE_FLOAT_SCALAR(void * array);

#define DEALLOCATE_ALLOCATABLE_FLOAT_1D FC_GLOBAL(atk_deallocate_allocatable_float_1d,ATK_DEALLOCATE_ALLOCATABLE_FLOAT_1D)
void DEALLOCATE_ALLOCATABLE_FLOAT_1D(void * array);

#define DEALLOCATE_ALLOCATABLE_DOUBLE_SCALAR FC_GLOBAL(atk_deallocate_allocatable_double_scalar,ATK_DEALLOCATE_ALLOCATABLE_DOUBLE_SCALAR)
void DEALLOCATE_ALLOCATABLE_DOUBLE_SCALAR(void * array);

#define DEALLOCATE_ALLOCATABLE_DOUBLE_1D FC_GLOBAL(atk_deallocate_allocatable_double_1d,ATK_DEALLOCATE_ALLOCATABLE_DOUBLE_1D)
void DEALLOCATE_ALLOCATABLE_DOUBLE_1D(void * array);

//[[[end]]]

//[[[cog
//gen.print_lines(cog.outl, gen.print_atk_reallocate_allocatable_header)
//]]]
#define REALLOCATE_ALLOCATABLE_INT_SCALAR FC_GLOBAL(atk_reallocate_allocatable_int_scalar,ATK_REALLOCATE_ALLOCATABLE_INT_SCALAR)
void REALLOCATE_ALLOCATABLE_INT_SCALAR(void * array, long nitems);

#define REALLOCATE_ALLOCATABLE_INT_1D FC_GLOBAL(atk_reallocate_allocatable_int_1d,ATK_REALLOCATE_ALLOCATABLE_INT_1D)
void REALLOCATE_ALLOCATABLE_INT_1D(void * array, long nitems);

#define REALLOCATE_ALLOCATABLE_LONG_SCALAR FC_GLOBAL(atk_reallocate_allocatable_long_scalar,ATK_REALLOCATE_ALLOCATABLE_LONG_SCALAR)
void REALLOCATE_ALLOCATABLE_LONG_SCALAR(void * array, long nitems);

#define REALLOCATE_ALLOCATABLE_LONG_1D FC_GLOBAL(atk_reallocate_allocatable_long_1d,ATK_REALLOCATE_ALLOCATABLE_LONG_1D)
void REALLOCATE_ALLOCATABLE_LONG_1D(void * array, long nitems);

#define REALLOCATE_ALLOCATABLE_FLOAT_SCALAR FC_GLOBAL(atk_reallocate_allocatable_float_scalar,ATK_REALLOCATE_ALLOCATABLE_FLOAT_SCALAR)
void REALLOCATE_ALLOCATABLE_FLOAT_SCALAR(void * array, long nitems);

#define REALLOCATE_ALLOCATABLE_FLOAT_1D FC_GLOBAL(atk_reallocate_allocatable_float_1d,ATK_REALLOCATE_ALLOCATABLE_FLOAT_1D)
void REALLOCATE_ALLOCATABLE_FLOAT_1D(void * array, long nitems);

#define REALLOCATE_ALLOCATABLE_DOUBLE_SCALAR FC_GLOBAL(atk_reallocate_allocatable_double_scalar,ATK_REALLOCATE_ALLOCATABLE_DOUBLE_SCALAR)
void REALLOCATE_ALLOCATABLE_DOUBLE_SCALAR(void * array, long nitems);

#define REALLOCATE_ALLOCATABLE_DOUBLE_1D FC_GLOBAL(atk_reallocate_allocatable_double_1d,ATK_REALLOCATE_ALLOCATABLE_DOUBLE_1D)
void REALLOCATE_ALLOCATABLE_DOUBLE_1D(void * array, long nitems);

//[[[end]]]

}

//----------------------------------------------------------------------

/*
 * Call a Fortran function to find the size of an allocatable.
 */
size_t SizeAllocatable(void * array, TypeID type, int rank)
//[[[cog
//gen.SizeAllocatable(cog.outl)
//]]]
{
  size_t nitems = 0;
  switch(type)
  {
  case CONDUIT_NATIVE_INT_DATATYPE_ID:
    switch(rank)
    {
    case 0:
      nitems = SIZE_ALLOCATABLE_INT_SCALAR(array);
      break;
    case 1:
      nitems = SIZE_ALLOCATABLE_INT_1D(array);
      break;
    default:
      break;
    }
    break;
  case CONDUIT_NATIVE_LONG_DATATYPE_ID:
    switch(rank)
    {
    case 0:
      nitems = SIZE_ALLOCATABLE_LONG_SCALAR(array);
      break;
    case 1:
      nitems = SIZE_ALLOCATABLE_LONG_1D(array);
      break;
    default:
      break;
    }
    break;
  case CONDUIT_NATIVE_FLOAT_DATATYPE_ID:
    switch(rank)
    {
    case 0:
      nitems = SIZE_ALLOCATABLE_FLOAT_SCALAR(array);
      break;
    case 1:
      nitems = SIZE_ALLOCATABLE_FLOAT_1D(array);
      break;
    default:
      break;
    }
    break;
  case CONDUIT_NATIVE_DOUBLE_DATATYPE_ID:
    switch(rank)
    {
    case 0:
      nitems = SIZE_ALLOCATABLE_DOUBLE_SCALAR(array);
      break;
    case 1:
      nitems = SIZE_ALLOCATABLE_DOUBLE_1D(array);
      break;
    default:
      break;
    }
    break;
  default:
    break;
  }
  return nitems;
}
//[[[end]]]

/*
 * Call a Fortran function to find the address of an allocatable.
 */
void * AddressAllocatable(void * array, TypeID type, int rank)
//[[[cog
//gen.AddressAllocatable(cog.outl)
//]]]
{
  void * addr = ATK_NULLPTR;
  switch(type)
  {
  case CONDUIT_NATIVE_INT_DATATYPE_ID:
    switch(rank)
    {
    case 0:
      ADDRESS_ALLOCATABLE_INT_SCALAR(array, &addr);
      break;
    case 1:
      ADDRESS_ALLOCATABLE_INT_1D(array, &addr);
      break;
    default:
      break;
    }
    break;
  case CONDUIT_NATIVE_LONG_DATATYPE_ID:
    switch(rank)
    {
    case 0:
      ADDRESS_ALLOCATABLE_LONG_SCALAR(array, &addr);
      break;
    case 1:
      ADDRESS_ALLOCATABLE_LONG_1D(array, &addr);
      break;
    default:
      break;
    }
    break;
  case CONDUIT_NATIVE_FLOAT_DATATYPE_ID:
    switch(rank)
    {
    case 0:
      ADDRESS_ALLOCATABLE_FLOAT_SCALAR(array, &addr);
      break;
    case 1:
      ADDRESS_ALLOCATABLE_FLOAT_1D(array, &addr);
      break;
    default:
      break;
    }
    break;
  case CONDUIT_NATIVE_DOUBLE_DATATYPE_ID:
    switch(rank)
    {
    case 0:
      ADDRESS_ALLOCATABLE_DOUBLE_SCALAR(array, &addr);
      break;
    case 1:
      ADDRESS_ALLOCATABLE_DOUBLE_1D(array, &addr);
      break;
    default:
      break;
    }
    break;
  default:
    break;
  }
  return addr;
}
//[[[end]]]

/*
 * Call a Fortran function to allocate an allocatable.
 */
void * AllocateAllocatable(void * array, TypeID type, int rank, SidreLength nitems)
//[[[cog
//gen.AllocateAllocatable(cog.outl)
//]]]
{
  void * addr = ATK_NULLPTR;
  switch(type)
  {
  case CONDUIT_NATIVE_INT_DATATYPE_ID:
    switch(rank)
    {
    case 0:
      ALLOCATE_ALLOCATABLE_INT_SCALAR(array, &nitems);
      ADDRESS_ALLOCATABLE_INT_SCALAR(array, &addr);
      break;
    case 1:
      ALLOCATE_ALLOCATABLE_INT_1D(array, &nitems);
      ADDRESS_ALLOCATABLE_INT_1D(array, &addr);
      break;
    default:
      break;
    }
    break;
  case CONDUIT_NATIVE_LONG_DATATYPE_ID:
    switch(rank)
    {
    case 0:
      ALLOCATE_ALLOCATABLE_LONG_SCALAR(array, &nitems);
      ADDRESS_ALLOCATABLE_LONG_SCALAR(array, &addr);
      break;
    case 1:
      ALLOCATE_ALLOCATABLE_LONG_1D(array, &nitems);
      ADDRESS_ALLOCATABLE_LONG_1D(array, &addr);
      break;
    default:
      break;
    }
    break;
  case CONDUIT_NATIVE_FLOAT_DATATYPE_ID:
    switch(rank)
    {
    case 0:
      ALLOCATE_ALLOCATABLE_FLOAT_SCALAR(array, &nitems);
      ADDRESS_ALLOCATABLE_FLOAT_SCALAR(array, &addr);
      break;
    case 1:
      ALLOCATE_ALLOCATABLE_FLOAT_1D(array, &nitems);
      ADDRESS_ALLOCATABLE_FLOAT_1D(array, &addr);
      break;
    default:
      break;
    }
    break;
  case CONDUIT_NATIVE_DOUBLE_DATATYPE_ID:
    switch(rank)
    {
    case 0:
      ALLOCATE_ALLOCATABLE_DOUBLE_SCALAR(array, &nitems);
      ADDRESS_ALLOCATABLE_DOUBLE_SCALAR(array, &addr);
      break;
    case 1:
      ALLOCATE_ALLOCATABLE_DOUBLE_1D(array, &nitems);
      ADDRESS_ALLOCATABLE_DOUBLE_1D(array, &addr);
      break;
    default:
      break;
    }
    break;
  default:
    break;
  }
  return addr;
}
//[[[end]]]

/*
 * Call a Fortran function to deallocate an allocatable.
 */
void DeallocateAllocatable(void * array, TypeID type, int rank)
//[[[cog
//gen.DeallocateAllocatable(cog.outl)
//]]]
{
  switch(type)
  {
  case CONDUIT_NATIVE_INT_DATATYPE_ID:
    switch(rank)
    {
    case 0:
      DEALLOCATE_ALLOCATABLE_INT_SCALAR(array);
      break;
    case 1:
      DEALLOCATE_ALLOCATABLE_INT_1D(array);
      break;
    default:
      break;
    }
    break;
  case CONDUIT_NATIVE_LONG_DATATYPE_ID:
    switch(rank)
    {
    case 0:
      DEALLOCATE_ALLOCATABLE_LONG_SCALAR(array);
      break;
    case 1:
      DEALLOCATE_ALLOCATABLE_LONG_1D(array);
      break;
    default:
      break;
    }
    break;
  case CONDUIT_NATIVE_FLOAT_DATATYPE_ID:
    switch(rank)
    {
    case 0:
      DEALLOCATE_ALLOCATABLE_FLOAT_SCALAR(array);
      break;
    case 1:
      DEALLOCATE_ALLOCATABLE_FLOAT_1D(array);
      break;
    default:
      break;
    }
    break;
  case CONDUIT_NATIVE_DOUBLE_DATATYPE_ID:
    switch(rank)
    {
    case 0:
      DEALLOCATE_ALLOCATABLE_DOUBLE_SCALAR(array);
      break;
    case 1:
      DEALLOCATE_ALLOCATABLE_DOUBLE_1D(array);
      break;
    default:
      break;
    }
    break;
  default:
    break;
  }
  return;
}
//[[[end]]]

extern "C" {

/*!
 * \brief Return DataView for a Fortran allocatable.
 *
 * The Fortran allocatable array is the buffer for the DataView.
 * Called from Fortran.
 */

void * ATK_create_fortran_allocatable_view(void * group,
                                           char * name, int lname,
                                           void * array, int type, int rank)
{
  DataGroup * grp = static_cast<DataGroup *>(group);
  DataView * view = grp->createFortranAllocatableView(std::string(name, lname),
                                                      array,
                                                      getTypeID(type), // XXX static_cast<TypeID>(type),
                                                      rank);
  return view;
}


/*!
 * \brief Return DataView for a Fortran array.
 *
 * Create an external view using information derived from a Fortran array.
 * Called from Fortran.
 */
void * ATK_create_array_view(void * group, char * name, int lname,
                             void * addr, int type, long nitems)
{
  DataGroup * grp = static_cast<DataGroup *>(group);
  DataView * view = grp->createExternalView(std::string(name, lname),
                                            addr,
                                            static_cast<TypeID>(type),
                                            static_cast<SidreLength>(nitems));
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
void FC_GLOBAL(atk_c_loc,ATK_C_LOC)(void * addr, void * * out)
{
  *out = addr;
}

/*
 *************************************************************************
 * These routines are called from Fortran with an interface
 * but without BIND(C)
 * since they need the address of the allocatable array,
 * not the address of the contents of the allocatable array.
 *
 * They are subroutines instead of functions since Intel has a different ABI
 * for returning derived types like type(CPTR). [They add an additional
 * argument at the beginning for the result]
 *
 * XXX - In the future it may be possible to replace them with one routine
 * with an interface like:
 *    type(*), allocatable :: array(..)
 *************************************************************************
 */
//[[[cog
//gen.print_lines(cog.outl, gen.print_atk_c_loc_allocatable)
//]]]

void FC_GLOBAL(atk_c_loc_allocatable_int_scalar,ATK_C_LOC_ALLOCATABLE_INT_SCALAR)
  (void * allocatable, void * * addr)
{
  *addr = allocatable;
}

void FC_GLOBAL(atk_c_loc_allocatable_int_1d,ATK_C_LOC_ALLOCATABLE_INT_1D)
  (void * allocatable, void * * addr)
{
  *addr = allocatable;
}

void FC_GLOBAL(atk_c_loc_allocatable_long_scalar,ATK_C_LOC_ALLOCATABLE_LONG_SCALAR)
  (void * allocatable, void * * addr)
{
  *addr = allocatable;
}

void FC_GLOBAL(atk_c_loc_allocatable_long_1d,ATK_C_LOC_ALLOCATABLE_LONG_1D)
  (void * allocatable, void * * addr)
{
  *addr = allocatable;
}

void FC_GLOBAL(atk_c_loc_allocatable_float_scalar,ATK_C_LOC_ALLOCATABLE_FLOAT_SCALAR)
  (void * allocatable, void * * addr)
{
  *addr = allocatable;
}

void FC_GLOBAL(atk_c_loc_allocatable_float_1d,ATK_C_LOC_ALLOCATABLE_FLOAT_1D)
  (void * allocatable, void * * addr)
{
  *addr = allocatable;
}

void FC_GLOBAL(atk_c_loc_allocatable_double_scalar,ATK_C_LOC_ALLOCATABLE_DOUBLE_SCALAR)
  (void * allocatable, void * * addr)
{
  *addr = allocatable;
}

void FC_GLOBAL(atk_c_loc_allocatable_double_1d,ATK_C_LOC_ALLOCATABLE_DOUBLE_1D)
  (void * allocatable, void * * addr)
{
  *addr = allocatable;
}
//[[[end]]]



}  // extern "C"


}  // namespace asctoolkit
}  // namespace sidre

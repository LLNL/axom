//
// SidreAllocatable.hpp - Routines used by Fortran interface
// Uses cog to insert some generated code into this file.
//
#ifndef SIDREALLOCATABLE_HPP_
#define SIDREALLOCATABLE_HPP_

extern "C" {
namespace asctoolkit
{
namespace sidre
{

size_t SizeAllocatable(void * array, TypeID type, int rank);

void * AddressAllocatable(void * array, TypeID type, int rank);

void * AllocateAllocatable(void * array, TypeID type, int rank, SidreLength nitems);

void DeallocateAllocatable(void * array, TypeID type, int rank);

void * ATK_create_fortran_allocatable_view(void * group,
                                           char * name, int lname,
                                           void * array, int type, int rank);

void * ATK_create_array_view(void * group, char * name, int lname,
                             void * addr, int type, long nitems);

}  // namespace asctoolkit
}  // namespace sidre
}  // extern "C"
#endif /* SIDREALLOCATABLE_HPP_ */

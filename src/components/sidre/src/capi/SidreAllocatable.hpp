//
// SidreAllocatable.hpp - Routines used by Fortran interface
// Uses cog to insert some generated code into this file.
//
#ifndef SIDREALLOCATABLE_HPP_
#define SIDREALLOCATABLE_HPP_

extern "C" {
namespace asctoolkit {
namespace sidre {

size_t SizeAllocatable(void * array, TypeID type, int rank);

void * AddressAllocatable(void * array, TypeID type, int rank);

void * AllocateAllocatable(void * array, TypeID type, int rank, SidreLength nitems);

void RegisterFortranAllocatableMetaBuffers(void);

}  // namespace asctoolkit
}  // namespace sidre
}  // extern "C"
#endif /* SIDREALLOCATABLE_HPP_ */


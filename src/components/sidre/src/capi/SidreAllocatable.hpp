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

void * SIDRE_create_array_view(void * group, char * name, int lname,
                               void * addr, int type, long nitems);

}  // namespace asctoolkit
}  // namespace sidre
}  // extern "C"
#endif /* SIDREALLOCATABLE_HPP_ */

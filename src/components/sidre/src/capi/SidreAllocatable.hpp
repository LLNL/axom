//
// SidreAllocatable.hpp - Routines used by Fortran interface
// Uses cog to insert some generated code into this file.
//
#ifndef SIDREALLOCATABLE_HPP_
#define SIDREALLOCATABLE_HPP_

extern "C" {
namespace asctoolkit {
namespace sidre {

//[[[cog
//import cog
//import genfsidresplicer as gen
//gen.print_lines(cog.outl, gen.print_atk_size_allocatable_header)
//]]]
size_t atk_size_allocatable_int_scalar_ptr_(void *array);
size_t atk_size_allocatable_int_1d_ptr_(void *array);
size_t atk_size_allocatable_long_scalar_ptr_(void *array);
size_t atk_size_allocatable_long_1d_ptr_(void *array);
size_t atk_size_allocatable_float_scalar_ptr_(void *array);
size_t atk_size_allocatable_float_1d_ptr_(void *array);
size_t atk_size_allocatable_double_scalar_ptr_(void *array);
size_t atk_size_allocatable_double_1d_ptr_(void *array);
//[[[end]]]

}  // namespace asctoolkit
}  // namespace sidre
}  // extern "C"
#endif /* SIDREALLOCATABLE_HPP_ */


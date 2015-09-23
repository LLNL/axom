//
// SidreAllocatable.cpp - Routines used by Fortran interface
// Uses cog to insert some generated code into this file.
//
#include <cstddef>
#include "common/CommonTypes.hpp"
#include "SidreWrapperHelpers.hpp"

extern "C" {
namespace asctoolkit {
namespace sidre {


//[[[cog
//import cog
//import genfsidresplicer as gen
//gen.print_lines(cog.outl, gen.print_atk_register_allocatable)
//]]]

// Fortran callable routine.
// Needed for each type-kind-rank to get address of allocatable array.
void *atk_register_allocatable_int_scalar_ptr_(
    DataGroup *group,
    char *name, int lname,
    void *array, int atk_type, int rank)
{
    return register_allocatable(group, std::string(name, lname), array, atk_type, rank); 
}

// Fortran callable routine.
// Needed for each type-kind-rank to get address of allocatable array.
void *atk_register_allocatable_int_1d_ptr_(
    DataGroup *group,
    char *name, int lname,
    void *array, int atk_type, int rank)
{
    return register_allocatable(group, std::string(name, lname), array, atk_type, rank); 
}

// Fortran callable routine.
// Needed for each type-kind-rank to get address of allocatable array.
void *atk_register_allocatable_long_scalar_ptr_(
    DataGroup *group,
    char *name, int lname,
    void *array, int atk_type, int rank)
{
    return register_allocatable(group, std::string(name, lname), array, atk_type, rank); 
}

// Fortran callable routine.
// Needed for each type-kind-rank to get address of allocatable array.
void *atk_register_allocatable_long_1d_ptr_(
    DataGroup *group,
    char *name, int lname,
    void *array, int atk_type, int rank)
{
    return register_allocatable(group, std::string(name, lname), array, atk_type, rank); 
}

// Fortran callable routine.
// Needed for each type-kind-rank to get address of allocatable array.
void *atk_register_allocatable_float_scalar_ptr_(
    DataGroup *group,
    char *name, int lname,
    void *array, int atk_type, int rank)
{
    return register_allocatable(group, std::string(name, lname), array, atk_type, rank); 
}

// Fortran callable routine.
// Needed for each type-kind-rank to get address of allocatable array.
void *atk_register_allocatable_float_1d_ptr_(
    DataGroup *group,
    char *name, int lname,
    void *array, int atk_type, int rank)
{
    return register_allocatable(group, std::string(name, lname), array, atk_type, rank); 
}

// Fortran callable routine.
// Needed for each type-kind-rank to get address of allocatable array.
void *atk_register_allocatable_double_scalar_ptr_(
    DataGroup *group,
    char *name, int lname,
    void *array, int atk_type, int rank)
{
    return register_allocatable(group, std::string(name, lname), array, atk_type, rank); 
}

// Fortran callable routine.
// Needed for each type-kind-rank to get address of allocatable array.
void *atk_register_allocatable_double_1d_ptr_(
    DataGroup *group,
    char *name, int lname,
    void *array, int atk_type, int rank)
{
    return register_allocatable(group, std::string(name, lname), array, atk_type, rank); 
}
//[[[end]]]


}  // namespace asctoolkit
}  // namespace sidre
}  // extern "C"


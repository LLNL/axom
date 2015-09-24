//
// SidreAllocatable.cpp - Routines used by Fortran interface
// Uses cog to insert some generated code into this file.
//
#include <cstddef>
#include "common/CommonTypes.hpp"
#include "SidreWrapperHelpers.hpp"

#include "sidre/DataGroup.hpp"
#include "sidre/MetaBuffer.hpp"

extern "C" {
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

//[[[cog
//gen.print_lines(cog.outl, gen.print_atk_address_allocatable_header)
//]]]
void *atk_address_allocatable_int_scalar_ptr_(void *array);
void *atk_address_allocatable_int_1d_ptr_(void *array);
void *atk_address_allocatable_long_scalar_ptr_(void *array);
void *atk_address_allocatable_long_1d_ptr_(void *array);
void *atk_address_allocatable_float_scalar_ptr_(void *array);
void *atk_address_allocatable_float_1d_ptr_(void *array);
void *atk_address_allocatable_double_scalar_ptr_(void *array);
void *atk_address_allocatable_double_1d_ptr_(void *array);
//[[[end]]]
}

namespace asctoolkit {
namespace sidre {

static int global_int;

// Holds pointers to Fortran functions since they cannot be in the
// Class AllocatableMetaBuffer directly.
struct Fptrs {
    size_t (*getNumberOfElements)(void *context);
    void *(*getDataPointer)(void *context);
};

class AllocatableMetaBuffer : public MetaBuffer
{
public:
  Fptrs *getFptrs() { return &m_callbacks; };

  virtual void *getDataPointer(void *context) const
    {
	return m_callbacks.getDataPointer(context);
    }

  virtual size_t getNumberOfElements(void *context) const
    {
	return m_callbacks.getNumberOfElements(context);
    }

private:
  Fptrs m_callbacks;

};

// XXX - temp code
AllocatableMetaBuffer *JUNK_DUMMY = NULL;
void RegisterFortranAllocatableMetaBuffers(void);

/*!
 * \brief Return DataView for a Fortran allocatable.
 *
 * The Fortran allocatable array is the buffer for the DataView.
 */
void *register_allocatable(DataGroup *group,
			   const std::string &name,
			   void *array, int atk_type, int rank)
{
  if (JUNK_DUMMY == NULL) {
      RegisterFortranAllocatableMetaBuffers();
  }
  global_int = atk_type;
  global_int = rank;
  return group->createViewWithMetaBuffer(name, array, JUNK_DUMMY);
}

extern "C" {

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

}  // extern "C"


// Create MetaBuffer classes to describe how to access Fortran allocatables.
// An instance is needed for each type-kind-rank since this depends on using
// the Fortran compiler to unpack the allocatable dope-vector.

void RegisterFortranAllocatableMetaBuffers(void)
{
    JUNK_DUMMY = new AllocatableMetaBuffer;
    Fptrs *callbacks = JUNK_DUMMY->getFptrs();

    callbacks->getNumberOfElements = atk_size_allocatable_int_1d_ptr_;
    callbacks->getDataPointer = atk_address_allocatable_int_1d_ptr_;
}


}  // namespace asctoolkit
}  // namespace sidre


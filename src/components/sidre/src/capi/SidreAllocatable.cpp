//
// SidreAllocatable.cpp - Routines used by Fortran interface
// Uses cog to insert some generated code into this file.
//
#include <cstddef>
#include "common/CommonTypes.hpp"
#include "SidreWrapperHelpers.hpp"

#include "sidre/DataGroup.hpp"
#include "sidre/MetaBuffer.hpp"
#include "sidre/SidreMalloc.hpp"  // for static variables

// import cog once
//[[[cog import cog;import genfsidresplicer as gen ]]]
//[[[end]]]

namespace asctoolkit {
namespace sidre {
class DataView;

extern "C" {
//[[[cog
//gen.print_lines(cog.outl, gen.print_atk_size_allocatable_header)
//]]]
size_t atk_size_allocatable_int_scalar_ptr_(void * array);
size_t atk_size_allocatable_int_1d_ptr_(void * array);
size_t atk_size_allocatable_long_scalar_ptr_(void * array);
size_t atk_size_allocatable_long_1d_ptr_(void * array);
size_t atk_size_allocatable_float_scalar_ptr_(void * array);
size_t atk_size_allocatable_float_1d_ptr_(void * array);
size_t atk_size_allocatable_double_scalar_ptr_(void * array);
size_t atk_size_allocatable_double_1d_ptr_(void * array);
//[[[end]]]

//[[[cog
//gen.print_lines(cog.outl, gen.print_atk_address_allocatable_header)
//]]]
void *atk_address_allocatable_int_scalar_ptr_(void * array);
void *atk_address_allocatable_int_1d_ptr_(void * array);
void *atk_address_allocatable_long_scalar_ptr_(void * array);
void *atk_address_allocatable_long_1d_ptr_(void * array);
void *atk_address_allocatable_float_scalar_ptr_(void * array);
void *atk_address_allocatable_float_1d_ptr_(void * array);
void *atk_address_allocatable_double_scalar_ptr_(void * array);
void *atk_address_allocatable_double_1d_ptr_(void * array);
//[[[end]]]

//[[[cog
//gen.print_lines(cog.outl, gen.print_atk_allocate_allocatable_header)
//]]]
void atk_allocate_allocatable_int_scalar_ptr_(void *array, long nitems);
void atk_allocate_allocatable_int_1d_ptr_(void *array, long nitems);
void atk_allocate_allocatable_long_scalar_ptr_(void *array, long nitems);
void atk_allocate_allocatable_long_1d_ptr_(void *array, long nitems);
void atk_allocate_allocatable_float_scalar_ptr_(void *array, long nitems);
void atk_allocate_allocatable_float_1d_ptr_(void *array, long nitems);
void atk_allocate_allocatable_double_scalar_ptr_(void *array, long nitems);
void atk_allocate_allocatable_double_1d_ptr_(void *array, long nitems);
//[[[end]]]

}

// Holds pointers to Fortran functions since they cannot be in the
// Class AllocatableMetaBuffer directly.
struct Fptrs {
    int rank;
    int type;
    size_t (*getNumberOfElements)(void *array);
    void *(*getDataPointer)(void *array);
    void (*allocate)(void *array, long nitems);
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

  virtual TypeID getTypeID(void *context) const
  {
      return static_cast<TypeID>(m_callbacks.type);
  }

    virtual void allocate(void *context,
			  TypeID type, SidreLength nitems) const
    {
    // XXX - type is fixed in the context, unused
    // XXX - check requested type vs context type
	m_callbacks.allocate(context, nitems);
	return;
    }


private:
  Fptrs m_callbacks;

};

// XXX - temp code
bool JUNK_DUMMY = false;
void RegisterFortranAllocatableMetaBuffers(void);

//[[[cog cog.outl('static AllocatableMetaBuffer *metabuffer_cache[%s];' % gen.num_metabuffers()) ]]]
static AllocatableMetaBuffer *metabuffer_cache[8];
//[[[end]]]

/*!
 * \brief Return DataView for a Fortran allocatable.
 *
 * The Fortran allocatable array is the buffer for the DataView.
 */
static void *register_allocatable(DataGroup *group,
				  const std::string &name,
				  void *context, int imetabuffer)
{
  if (! JUNK_DUMMY) {
      RegisterFortranAllocatableMetaBuffers();
      JUNK_DUMMY = true;
  }
  return group->createViewWithMetaBuffer(name, context,
					 metabuffer_cache[imetabuffer]);
}

extern "C" {

// called from Fortran
void *ATK_register_static(void *group, char *name, int lname, void *addr, int type, long nitems)
{
    return registerStaticNode(
			      static_cast<DataGroup *>(group),
			      std::string(name, lname),
			      addr,
			      static_cast<TypeID>(type),
			      static_cast<SidreLength>(nitems));
}


//[[[cog
//gen.print_lines(cog.outl, gen.print_atk_register_allocatable)
//]]]

// Fortran callable routine.
// Needed for each type-kind-rank to get address of allocatable array.
// array is address of allocatable, not the result of C_LOC(array)
void *atk_register_allocatable_int_scalar_ptr_(
    DataGroup *group,
    char *name, int lname,
    void *array, int indx)
{
    return register_allocatable(group, std::string(name, lname), array, indx); 
}

// Fortran callable routine.
// Needed for each type-kind-rank to get address of allocatable array.
// array is address of allocatable, not the result of C_LOC(array)
void *atk_register_allocatable_int_1d_ptr_(
    DataGroup *group,
    char *name, int lname,
    void *array, int indx)
{
    return register_allocatable(group, std::string(name, lname), array, indx); 
}

// Fortran callable routine.
// Needed for each type-kind-rank to get address of allocatable array.
// array is address of allocatable, not the result of C_LOC(array)
void *atk_register_allocatable_long_scalar_ptr_(
    DataGroup *group,
    char *name, int lname,
    void *array, int indx)
{
    return register_allocatable(group, std::string(name, lname), array, indx); 
}

// Fortran callable routine.
// Needed for each type-kind-rank to get address of allocatable array.
// array is address of allocatable, not the result of C_LOC(array)
void *atk_register_allocatable_long_1d_ptr_(
    DataGroup *group,
    char *name, int lname,
    void *array, int indx)
{
    return register_allocatable(group, std::string(name, lname), array, indx); 
}

// Fortran callable routine.
// Needed for each type-kind-rank to get address of allocatable array.
// array is address of allocatable, not the result of C_LOC(array)
void *atk_register_allocatable_float_scalar_ptr_(
    DataGroup *group,
    char *name, int lname,
    void *array, int indx)
{
    return register_allocatable(group, std::string(name, lname), array, indx); 
}

// Fortran callable routine.
// Needed for each type-kind-rank to get address of allocatable array.
// array is address of allocatable, not the result of C_LOC(array)
void *atk_register_allocatable_float_1d_ptr_(
    DataGroup *group,
    char *name, int lname,
    void *array, int indx)
{
    return register_allocatable(group, std::string(name, lname), array, indx); 
}

// Fortran callable routine.
// Needed for each type-kind-rank to get address of allocatable array.
// array is address of allocatable, not the result of C_LOC(array)
void *atk_register_allocatable_double_scalar_ptr_(
    DataGroup *group,
    char *name, int lname,
    void *array, int indx)
{
    return register_allocatable(group, std::string(name, lname), array, indx); 
}

// Fortran callable routine.
// Needed for each type-kind-rank to get address of allocatable array.
// array is address of allocatable, not the result of C_LOC(array)
void *atk_register_allocatable_double_1d_ptr_(
    DataGroup *group,
    char *name, int lname,
    void *array, int indx)
{
    return register_allocatable(group, std::string(name, lname), array, indx); 
}
//[[[end]]]

}  // extern "C"


// Create MetaBuffer classes to describe how to access Fortran allocatables.
// An instance is needed for each type-kind-rank since this depends on using
// the Fortran compiler to unpack the allocatable dope-vector.

void RegisterFortranAllocatableMetaBuffers(void)
{
  Fptrs *callbacks;


//[[[cog
//gen.print_lines(cog.out, gen.print_metabuffer)
//]]]

metabuffer_cache[0] = new AllocatableMetaBuffer;
callbacks = metabuffer_cache[0]->getFptrs();
callbacks->rank = 0;
callbacks->type = ATK_C_INT_T;
callbacks->getNumberOfElements = atk_size_allocatable_int_scalar_ptr_;
callbacks->getDataPointer = atk_address_allocatable_int_scalar_ptr_;
callbacks->allocate = atk_allocate_allocatable_int_scalar_ptr_;

metabuffer_cache[1] = new AllocatableMetaBuffer;
callbacks = metabuffer_cache[1]->getFptrs();
callbacks->rank = 1;
callbacks->type = ATK_C_INT_T;
callbacks->getNumberOfElements = atk_size_allocatable_int_1d_ptr_;
callbacks->getDataPointer = atk_address_allocatable_int_1d_ptr_;
callbacks->allocate = atk_allocate_allocatable_int_1d_ptr_;

metabuffer_cache[2] = new AllocatableMetaBuffer;
callbacks = metabuffer_cache[2]->getFptrs();
callbacks->rank = 0;
callbacks->type = ATK_C_LONG_T;
callbacks->getNumberOfElements = atk_size_allocatable_long_scalar_ptr_;
callbacks->getDataPointer = atk_address_allocatable_long_scalar_ptr_;
callbacks->allocate = atk_allocate_allocatable_long_scalar_ptr_;

metabuffer_cache[3] = new AllocatableMetaBuffer;
callbacks = metabuffer_cache[3]->getFptrs();
callbacks->rank = 1;
callbacks->type = ATK_C_LONG_T;
callbacks->getNumberOfElements = atk_size_allocatable_long_1d_ptr_;
callbacks->getDataPointer = atk_address_allocatable_long_1d_ptr_;
callbacks->allocate = atk_allocate_allocatable_long_1d_ptr_;

metabuffer_cache[4] = new AllocatableMetaBuffer;
callbacks = metabuffer_cache[4]->getFptrs();
callbacks->rank = 0;
callbacks->type = ATK_C_FLOAT_T;
callbacks->getNumberOfElements = atk_size_allocatable_float_scalar_ptr_;
callbacks->getDataPointer = atk_address_allocatable_float_scalar_ptr_;
callbacks->allocate = atk_allocate_allocatable_float_scalar_ptr_;

metabuffer_cache[5] = new AllocatableMetaBuffer;
callbacks = metabuffer_cache[5]->getFptrs();
callbacks->rank = 1;
callbacks->type = ATK_C_FLOAT_T;
callbacks->getNumberOfElements = atk_size_allocatable_float_1d_ptr_;
callbacks->getDataPointer = atk_address_allocatable_float_1d_ptr_;
callbacks->allocate = atk_allocate_allocatable_float_1d_ptr_;

metabuffer_cache[6] = new AllocatableMetaBuffer;
callbacks = metabuffer_cache[6]->getFptrs();
callbacks->rank = 0;
callbacks->type = ATK_C_DOUBLE_T;
callbacks->getNumberOfElements = atk_size_allocatable_double_scalar_ptr_;
callbacks->getDataPointer = atk_address_allocatable_double_scalar_ptr_;
callbacks->allocate = atk_allocate_allocatable_double_scalar_ptr_;

metabuffer_cache[7] = new AllocatableMetaBuffer;
callbacks = metabuffer_cache[7]->getFptrs();
callbacks->rank = 1;
callbacks->type = ATK_C_DOUBLE_T;
callbacks->getNumberOfElements = atk_size_allocatable_double_1d_ptr_;
callbacks->getDataPointer = atk_address_allocatable_double_1d_ptr_;
callbacks->allocate = atk_allocate_allocatable_double_1d_ptr_;
//[[[end]]]
}

}  // namespace asctoolkit
}  // namespace sidre


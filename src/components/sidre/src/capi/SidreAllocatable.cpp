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

//[[[cog
//gen.print_lines(cog.outl, gen.print_atk_deallocate_allocatable_header)
//]]]
void atk_deallocate_allocatable_int_scalar_ptr_(void *array);
void atk_deallocate_allocatable_int_1d_ptr_(void *array);
void atk_deallocate_allocatable_long_scalar_ptr_(void *array);
void atk_deallocate_allocatable_long_1d_ptr_(void *array);
void atk_deallocate_allocatable_float_scalar_ptr_(void *array);
void atk_deallocate_allocatable_float_1d_ptr_(void *array);
void atk_deallocate_allocatable_double_scalar_ptr_(void *array);
void atk_deallocate_allocatable_double_1d_ptr_(void *array);
//[[[end]]]

//[[[cog
//gen.print_lines(cog.outl, gen.print_atk_reallocate_allocatable_header)
//]]]
void atk_reallocate_allocatable_int_scalar_ptr_(void *array, long nitems);
void atk_reallocate_allocatable_int_1d_ptr_(void *array, long nitems);
void atk_reallocate_allocatable_long_scalar_ptr_(void *array, long nitems);
void atk_reallocate_allocatable_long_1d_ptr_(void *array, long nitems);
void atk_reallocate_allocatable_float_scalar_ptr_(void *array, long nitems);
void atk_reallocate_allocatable_float_1d_ptr_(void *array, long nitems);
void atk_reallocate_allocatable_double_scalar_ptr_(void *array, long nitems);
void atk_reallocate_allocatable_double_1d_ptr_(void *array, long nitems);
//[[[end]]]

}

//----------------------------------------------------------------------
// Holds pointers to Fortran functions since they cannot be in the
// Class AllocatableMetaBuffer directly.
struct Fptrs {
    int rank;
    int type;
    size_t (*getNumberOfElements)(void *array);
    void *(*getDataPointer)(void *array);
    void (*allocate)(void *array, long nitems);
    void (*deallocate)(void *array);
    void (*reallocate)(void *array, long nitems);
};

//[[[cog cog.outl('static struct Fptrs fptrs_cache[%s] = {' % gen.num_metabuffers()) ]]]
static struct Fptrs fptrs_cache[8] = {
//[[[end]]]

//[[[cog
//gen.print_lines(cog.out, gen.print_fptrs_cache)
//]]]

{
  0,   // rank
  ATK_C_INT_T,
  atk_size_allocatable_int_scalar_ptr_,
  atk_address_allocatable_int_scalar_ptr_,
  atk_allocate_allocatable_int_scalar_ptr_,
  atk_deallocate_allocatable_int_scalar_ptr_,
  atk_reallocate_allocatable_int_scalar_ptr_
},

{
  1,   // rank
  ATK_C_INT_T,
  atk_size_allocatable_int_1d_ptr_,
  atk_address_allocatable_int_1d_ptr_,
  atk_allocate_allocatable_int_1d_ptr_,
  atk_deallocate_allocatable_int_1d_ptr_,
  atk_reallocate_allocatable_int_1d_ptr_
},

{
  0,   // rank
  ATK_C_LONG_T,
  atk_size_allocatable_long_scalar_ptr_,
  atk_address_allocatable_long_scalar_ptr_,
  atk_allocate_allocatable_long_scalar_ptr_,
  atk_deallocate_allocatable_long_scalar_ptr_,
  atk_reallocate_allocatable_long_scalar_ptr_
},

{
  1,   // rank
  ATK_C_LONG_T,
  atk_size_allocatable_long_1d_ptr_,
  atk_address_allocatable_long_1d_ptr_,
  atk_allocate_allocatable_long_1d_ptr_,
  atk_deallocate_allocatable_long_1d_ptr_,
  atk_reallocate_allocatable_long_1d_ptr_
},

{
  0,   // rank
  ATK_C_FLOAT_T,
  atk_size_allocatable_float_scalar_ptr_,
  atk_address_allocatable_float_scalar_ptr_,
  atk_allocate_allocatable_float_scalar_ptr_,
  atk_deallocate_allocatable_float_scalar_ptr_,
  atk_reallocate_allocatable_float_scalar_ptr_
},

{
  1,   // rank
  ATK_C_FLOAT_T,
  atk_size_allocatable_float_1d_ptr_,
  atk_address_allocatable_float_1d_ptr_,
  atk_allocate_allocatable_float_1d_ptr_,
  atk_deallocate_allocatable_float_1d_ptr_,
  atk_reallocate_allocatable_float_1d_ptr_
},

{
  0,   // rank
  ATK_C_DOUBLE_T,
  atk_size_allocatable_double_scalar_ptr_,
  atk_address_allocatable_double_scalar_ptr_,
  atk_allocate_allocatable_double_scalar_ptr_,
  atk_deallocate_allocatable_double_scalar_ptr_,
  atk_reallocate_allocatable_double_scalar_ptr_
},

{
  1,   // rank
  ATK_C_DOUBLE_T,
  atk_size_allocatable_double_1d_ptr_,
  atk_address_allocatable_double_1d_ptr_,
  atk_allocate_allocatable_double_1d_ptr_,
  atk_deallocate_allocatable_double_1d_ptr_,
  atk_reallocate_allocatable_double_1d_ptr_
},
//[[[end]]]
};

//----------------------------------------------------------------------

class AllocatableMetaBuffer : public MetaBuffer
{
public:
  virtual void *getDataPointer() const
    {
	return m_callbacks->getDataPointer(m_context);
    }

  virtual size_t getNumberOfElements() const
    {
	return m_callbacks->getNumberOfElements(m_context);
    }

  virtual TypeID getTypeID() const
  {
      return static_cast<TypeID>(m_callbacks->type);
  }

    virtual void * allocate(TypeID type, SidreLength nitems) const
    {
    // XXX - type is fixed in the context, unused
    // XXX - check requested type vs context type
	m_callbacks->allocate(m_context, nitems);
	return m_callbacks->getDataPointer(m_context);
    }

  AllocatableMetaBuffer(void *context, const Fptrs * callbacks) :
    m_context(context),
    m_callbacks(callbacks)
  { }

  void setFptrs(const Fptrs * callbacks)
  {
    m_callbacks = callbacks;
  }


private:
  void * m_context;   // pointer to Fortran allocatable
  const Fptrs * m_callbacks;
};

/*!
 * \brief Return DataView for a Fortran allocatable.
 *
 * The Fortran allocatable array is the buffer for the DataView.
 */
static void *register_allocatable(DataGroup *group,
				  const std::string &name,
				  void *context, int imetabuffer)
{
  AllocatableMetaBuffer * metabuffer = new AllocatableMetaBuffer(context, fptrs_cache + imetabuffer);
  return group->createMetaBufferView(name, metabuffer);
}

extern "C" {

// called from Fortran
// return pointer to a DataView
void *ATK_register_static(void *group, char *name, int lname, void *addr, int type, long nitems)
{
  DataGroup *grp = static_cast<DataGroup *>(group);
  return grp->createExternalView( std::string(name, lname),
				  addr,
				  static_cast<TypeID>(type),
				  static_cast<SidreLength>(nitems));
}

// called from Fortran
// https://gcc.gnu.org/bugzilla/show_bug.cgi?id=53945
// Work around a problem with gfortran 4.7 where C_LOC does not work
// with assumed shape array.  Passing the first element of the
// array to a function without an interface will force the compiler
// to use f77 semantics and pass the address of the data, essentially
// the same as C_LOC.
// XXX Pass the first element, not the entire array, to avoid getting
// XXX a copy of the array.
void *atk_c_loc_(void *addr)
{
    return addr;
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


}  // namespace asctoolkit
}  // namespace sidre


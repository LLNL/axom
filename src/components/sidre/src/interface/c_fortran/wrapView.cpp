// wrapView.cpp
// This is generated code, do not edit
//
// Copyright (c) 2015, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory.
//
// All rights reserved.
//
// This source code cannot be distributed without permission and
// further review from Lawrence Livermore National Laboratory.
//
// wrapView.cpp
#include "wrapView.h"
#include <string>
#include "shroudrt.hpp"
#include "sidre/SidreTypes.hpp"
#include "sidre/View.hpp"

extern "C" {
namespace axom {
namespace sidre {

void SIDRE_view_allocate_simple(SIDRE_view * self)
{
View *selfobj = static_cast<View *>(static_cast<void *>(self));
// splicer begin class.View.method.allocate_simple
    selfobj->allocate();
    return;
// splicer end class.View.method.allocate_simple
}

void SIDRE_view_allocate_from_type(SIDRE_view * self, int type, SIDRE_SidreLength num_elems)
{
View *selfobj = static_cast<View *>(static_cast<void *>(self));
// splicer begin class.View.method.allocate_from_type
    selfobj->allocate(getTypeID(type), num_elems);
    return;
// splicer end class.View.method.allocate_from_type
}

void SIDRE_view_reallocate(SIDRE_view * self, SIDRE_SidreLength num_elems)
{
View *selfobj = static_cast<View *>(static_cast<void *>(self));
// splicer begin class.View.method.reallocate
    selfobj->reallocate(num_elems);
    return;
// splicer end class.View.method.reallocate
}

void SIDRE_view_apply_0(SIDRE_view * self)
{
View *selfobj = static_cast<View *>(static_cast<void *>(self));
// splicer begin class.View.method.apply_0
    selfobj->apply();
    return;
// splicer end class.View.method.apply_0
}

void SIDRE_view_attach_buffer_only(SIDRE_view * self, SIDRE_buffer * buff)
{
View *selfobj = static_cast<View *>(static_cast<void *>(self));
// splicer begin class.View.method.attach_buffer_only
    selfobj->attachBuffer(static_cast<Buffer *>(static_cast<void *>(buff)));
    return;
// splicer end class.View.method.attach_buffer_only
}

void SIDRE_view_attach_buffer_type(SIDRE_view * self, int type, SIDRE_SidreLength num_elems, SIDRE_buffer * buff)
{
View *selfobj = static_cast<View *>(static_cast<void *>(self));
// splicer begin class.View.method.attach_buffer_type
    selfobj->attachBuffer(getTypeID(type), num_elems, static_cast<Buffer *>(static_cast<void *>(buff)));
    return;
// splicer end class.View.method.attach_buffer_type
}

void SIDRE_view_attach_buffer_shape(SIDRE_view * self, int type, int ndims, SIDRE_SidreLength * shape, SIDRE_buffer * buff)
{
View *selfobj = static_cast<View *>(static_cast<void *>(self));
// splicer begin class.View.method.attach_buffer_shape
    selfobj->attachBuffer(getTypeID(type), ndims, shape, static_cast<Buffer *>(static_cast<void *>(buff)));
    return;
// splicer end class.View.method.attach_buffer_shape
}

void SIDRE_view_apply_nelems(SIDRE_view * self, SIDRE_SidreLength num_elems)
{
View *selfobj = static_cast<View *>(static_cast<void *>(self));
// splicer begin class.View.method.apply_nelems
    selfobj->apply(num_elems);
    return;
// splicer end class.View.method.apply_nelems
}

void SIDRE_view_apply_nelems_offset(SIDRE_view * self, SIDRE_SidreLength num_elems, SIDRE_SidreLength offset)
{
View *selfobj = static_cast<View *>(static_cast<void *>(self));
// splicer begin class.View.method.apply_nelems_offset
    selfobj->apply(num_elems, offset);
    return;
// splicer end class.View.method.apply_nelems_offset
}

void SIDRE_view_apply_nelems_offset_stride(SIDRE_view * self, SIDRE_SidreLength num_elems, SIDRE_SidreLength offset, SIDRE_SidreLength stride)
{
View *selfobj = static_cast<View *>(static_cast<void *>(self));
// splicer begin class.View.method.apply_nelems_offset_stride
    selfobj->apply(num_elems, offset, stride);
    return;
// splicer end class.View.method.apply_nelems_offset_stride
}

void SIDRE_view_apply_type_nelems(SIDRE_view * self, int type, SIDRE_SidreLength num_elems)
{
View *selfobj = static_cast<View *>(static_cast<void *>(self));
// splicer begin class.View.method.apply_type_nelems
    selfobj->apply(getTypeID(type), num_elems);
    return;
// splicer end class.View.method.apply_type_nelems
}

void SIDRE_view_apply_type_nelems_offset(SIDRE_view * self, int type, SIDRE_SidreLength num_elems, SIDRE_SidreLength offset)
{
View *selfobj = static_cast<View *>(static_cast<void *>(self));
// splicer begin class.View.method.apply_type_nelems_offset
    selfobj->apply(getTypeID(type), num_elems, offset);
    return;
// splicer end class.View.method.apply_type_nelems_offset
}

void SIDRE_view_apply_type_nelems_offset_stride(SIDRE_view * self, int type, SIDRE_SidreLength num_elems, SIDRE_SidreLength offset, SIDRE_SidreLength stride)
{
View *selfobj = static_cast<View *>(static_cast<void *>(self));
// splicer begin class.View.method.apply_type_nelems_offset_stride
    selfobj->apply(getTypeID(type), num_elems, offset, stride);
    return;
// splicer end class.View.method.apply_type_nelems_offset_stride
}

void SIDRE_view_apply_type_shape(SIDRE_view * self, int type, int ndims, SIDRE_SidreLength * shape)
{
View *selfobj = static_cast<View *>(static_cast<void *>(self));
// splicer begin class.View.method.apply_type_shape
    selfobj->apply(getTypeID(type), ndims, shape);
    return;
// splicer end class.View.method.apply_type_shape
}

bool SIDRE_view_has_buffer(const SIDRE_view * self)
{
const View *selfobj = static_cast<const View *>(static_cast<const void *>(self));
// splicer begin class.View.method.has_buffer
    bool rv = selfobj->hasBuffer();
    return rv;
// splicer end class.View.method.has_buffer
}

bool SIDRE_view_is_external(const SIDRE_view * self)
{
const View *selfobj = static_cast<const View *>(static_cast<const void *>(self));
// splicer begin class.View.method.is_external
    bool rv = selfobj->isExternal();
    return rv;
// splicer end class.View.method.is_external
}

bool SIDRE_view_is_allocated(SIDRE_view * self)
{
View *selfobj = static_cast<View *>(static_cast<void *>(self));
// splicer begin class.View.method.is_allocated
    bool rv = selfobj->isAllocated();
    return rv;
// splicer end class.View.method.is_allocated
}

bool SIDRE_view_is_applied(const SIDRE_view * self)
{
const View *selfobj = static_cast<const View *>(static_cast<const void *>(self));
// splicer begin class.View.method.is_applied
    bool rv = selfobj->isApplied();
    return rv;
// splicer end class.View.method.is_applied
}

bool SIDRE_view_is_described(const SIDRE_view * self)
{
const View *selfobj = static_cast<const View *>(static_cast<const void *>(self));
// splicer begin class.View.method.is_described
    bool rv = selfobj->isDescribed();
    return rv;
// splicer end class.View.method.is_described
}

bool SIDRE_view_is_empty(const SIDRE_view * self)
{
const View *selfobj = static_cast<const View *>(static_cast<const void *>(self));
// splicer begin class.View.method.is_empty
    bool rv = selfobj->isEmpty();
    return rv;
// splicer end class.View.method.is_empty
}

bool SIDRE_view_is_opaque(const SIDRE_view * self)
{
const View *selfobj = static_cast<const View *>(static_cast<const void *>(self));
// splicer begin class.View.method.is_opaque
    bool rv = selfobj->isOpaque();
    return rv;
// splicer end class.View.method.is_opaque
}

bool SIDRE_view_is_scalar(const SIDRE_view * self)
{
const View *selfobj = static_cast<const View *>(static_cast<const void *>(self));
// splicer begin class.View.method.is_scalar
    bool rv = selfobj->isScalar();
    return rv;
// splicer end class.View.method.is_scalar
}

bool SIDRE_view_is_string(const SIDRE_view * self)
{
const View *selfobj = static_cast<const View *>(static_cast<const void *>(self));
// splicer begin class.View.method.is_string
    bool rv = selfobj->isString();
    return rv;
// splicer end class.View.method.is_string
}

const char * SIDRE_view_get_name(const SIDRE_view * self)
{
const View *selfobj = static_cast<const View *>(static_cast<const void *>(self));
// splicer begin class.View.method.get_name
    const std::string & rv = selfobj->getName();
    return rv.c_str();
// splicer end class.View.method.get_name
}

void SIDRE_view_get_name_bufferify(const SIDRE_view * self, char * SH_F_rv, int LSH_F_rv)
{
const View *selfobj = static_cast<const View *>(static_cast<const void *>(self));
// splicer begin class.View.method.get_name_bufferify
    const std::string & rv = selfobj->getName();
    shroud::FccCopy(SH_F_rv, LSH_F_rv, rv.c_str());
    return;
// splicer end class.View.method.get_name_bufferify
}

SIDRE_buffer * SIDRE_view_get_buffer(SIDRE_view * self)
{
View *selfobj = static_cast<View *>(static_cast<void *>(self));
// splicer begin class.View.method.get_buffer
    Buffer * rv = selfobj->getBuffer();
    return static_cast<SIDRE_buffer *>(static_cast<void *>(rv));
// splicer end class.View.method.get_buffer
}

void * SIDRE_view_get_void_ptr(const SIDRE_view * self)
{
const View *selfobj = static_cast<const View *>(static_cast<const void *>(self));
// splicer begin class.View.method.get_void_ptr
    void * rv = selfobj->getVoidPtr();
    return rv;
// splicer end class.View.method.get_void_ptr
}

void SIDRE_view_set_scalar_int(SIDRE_view * self, int value)
{
View *selfobj = static_cast<View *>(static_cast<void *>(self));
// splicer begin class.View.method.set_scalar_int
    selfobj->setScalar<int>(value);
    return;
// splicer end class.View.method.set_scalar_int
}

void SIDRE_view_set_scalar_long(SIDRE_view * self, long value)
{
View *selfobj = static_cast<View *>(static_cast<void *>(self));
// splicer begin class.View.method.set_scalar_long
    selfobj->setScalar<long>(value);
    return;
// splicer end class.View.method.set_scalar_long
}

void SIDRE_view_set_scalar_float(SIDRE_view * self, float value)
{
View *selfobj = static_cast<View *>(static_cast<void *>(self));
// splicer begin class.View.method.set_scalar_float
    selfobj->setScalar<float>(value);
    return;
// splicer end class.View.method.set_scalar_float
}

void SIDRE_view_set_scalar_double(SIDRE_view * self, double value)
{
View *selfobj = static_cast<View *>(static_cast<void *>(self));
// splicer begin class.View.method.set_scalar_double
    selfobj->setScalar<double>(value);
    return;
// splicer end class.View.method.set_scalar_double
}

void SIDRE_view_set_external_data_ptr_only(SIDRE_view * self, void * external_ptr)
{
View *selfobj = static_cast<View *>(static_cast<void *>(self));
// splicer begin class.View.method.set_external_data_ptr_only
    selfobj->setExternalDataPtr(external_ptr);
    return;
// splicer end class.View.method.set_external_data_ptr_only
}

void SIDRE_view_set_external_data_ptr_type(SIDRE_view * self, int type, SIDRE_SidreLength num_elems, void * external_ptr)
{
View *selfobj = static_cast<View *>(static_cast<void *>(self));
// splicer begin class.View.method.set_external_data_ptr_type
    selfobj->setExternalDataPtr(getTypeID(type), num_elems, external_ptr);
    return;
// splicer end class.View.method.set_external_data_ptr_type
}

void SIDRE_view_set_string(SIDRE_view * self, const char * value)
{
View *selfobj = static_cast<View *>(static_cast<void *>(self));
// splicer begin class.View.method.set_string
    const std::string SH_value(value);
    selfobj->setString(SH_value);
    return;
// splicer end class.View.method.set_string
}

void SIDRE_view_set_string_bufferify(SIDRE_view * self, const char * value, int Lvalue)
{
View *selfobj = static_cast<View *>(static_cast<void *>(self));
// splicer begin class.View.method.set_string_bufferify
    const std::string SH_value(value, Lvalue);
    selfobj->setString(SH_value);
    return;
// splicer end class.View.method.set_string_bufferify
}

void SIDRE_view_set_external_data_ptr_shape(SIDRE_view * self, int type, int ndims, SIDRE_SidreLength * shape, void * external_ptr)
{
View *selfobj = static_cast<View *>(static_cast<void *>(self));
// splicer begin class.View.method.set_external_data_ptr_shape
    selfobj->setExternalDataPtr(getTypeID(type), ndims, shape, external_ptr);
    return;
// splicer end class.View.method.set_external_data_ptr_shape
}

const char * SIDRE_view_get_string(SIDRE_view * self)
{
View *selfobj = static_cast<View *>(static_cast<void *>(self));
// splicer begin class.View.method.get_string
    const char * rv = selfobj->getString();
    return rv;
// splicer end class.View.method.get_string
}

void SIDRE_view_get_string_bufferify(SIDRE_view * self, char * name, int Lname)
{
View *selfobj = static_cast<View *>(static_cast<void *>(self));
// splicer begin class.View.method.get_string_bufferify
    const char * rv = selfobj->getString();
    shroud::FccCopy(name, Lname, rv);
    return;
// splicer end class.View.method.get_string_bufferify
}

int SIDRE_view_get_data_int(SIDRE_view * self)
{
View *selfobj = static_cast<View *>(static_cast<void *>(self));
// splicer begin class.View.method.get_data_int
    int rv = selfobj->getData<int>();
    return rv;
// splicer end class.View.method.get_data_int
}

long SIDRE_view_get_data_long(SIDRE_view * self)
{
View *selfobj = static_cast<View *>(static_cast<void *>(self));
// splicer begin class.View.method.get_data_long
    long rv = selfobj->getData<long>();
    return rv;
// splicer end class.View.method.get_data_long
}

float SIDRE_view_get_data_float(SIDRE_view * self)
{
View *selfobj = static_cast<View *>(static_cast<void *>(self));
// splicer begin class.View.method.get_data_float
    float rv = selfobj->getData<float>();
    return rv;
// splicer end class.View.method.get_data_float
}

double SIDRE_view_get_data_double(SIDRE_view * self)
{
View *selfobj = static_cast<View *>(static_cast<void *>(self));
// splicer begin class.View.method.get_data_double
    double rv = selfobj->getData<double>();
    return rv;
// splicer end class.View.method.get_data_double
}

SIDRE_group * SIDRE_view_get_owning_group(SIDRE_view * self)
{
View *selfobj = static_cast<View *>(static_cast<void *>(self));
// splicer begin class.View.method.get_owning_group
    Group * rv = selfobj->getOwningGroup();
    return static_cast<SIDRE_group *>(static_cast<void *>(rv));
// splicer end class.View.method.get_owning_group
}

int SIDRE_view_get_type_id(const SIDRE_view * self)
{
const View *selfobj = static_cast<const View *>(static_cast<const void *>(self));
// splicer begin class.View.method.get_type_id
    TypeID rv = selfobj->getTypeID();
    return static_cast<int>(rv);
// splicer end class.View.method.get_type_id
}

size_t SIDRE_view_get_total_bytes(const SIDRE_view * self)
{
const View *selfobj = static_cast<const View *>(static_cast<const void *>(self));
// splicer begin class.View.method.get_total_bytes
    size_t rv = selfobj->getTotalBytes();
    return rv;
// splicer end class.View.method.get_total_bytes
}

size_t SIDRE_view_get_bytes_per_element(const SIDRE_view * self)
{
const View *selfobj = static_cast<const View *>(static_cast<const void *>(self));
// splicer begin class.View.method.get_bytes_per_element
    size_t rv = selfobj->getBytesPerElement();
    return rv;
// splicer end class.View.method.get_bytes_per_element
}

size_t SIDRE_view_get_num_elements(const SIDRE_view * self)
{
const View *selfobj = static_cast<const View *>(static_cast<const void *>(self));
// splicer begin class.View.method.get_num_elements
    size_t rv = selfobj->getNumElements();
    return rv;
// splicer end class.View.method.get_num_elements
}

size_t SIDRE_view_get_offset(const SIDRE_view * self)
{
const View *selfobj = static_cast<const View *>(static_cast<const void *>(self));
// splicer begin class.View.method.get_offset
    size_t rv = selfobj->getOffset();
    return rv;
// splicer end class.View.method.get_offset
}

size_t SIDRE_view_get_stride(const SIDRE_view * self)
{
const View *selfobj = static_cast<const View *>(static_cast<const void *>(self));
// splicer begin class.View.method.get_stride
    size_t rv = selfobj->getStride();
    return rv;
// splicer end class.View.method.get_stride
}

int SIDRE_view_get_num_dimensions(const SIDRE_view * self)
{
const View *selfobj = static_cast<const View *>(static_cast<const void *>(self));
// splicer begin class.View.method.get_num_dimensions
    int rv = selfobj->getNumDimensions();
    return rv;
// splicer end class.View.method.get_num_dimensions
}

int SIDRE_view_get_shape(const SIDRE_view * self, int ndims, SIDRE_SidreLength * shape)
{
const View *selfobj = static_cast<const View *>(static_cast<const void *>(self));
// splicer begin class.View.method.get_shape
    int rv = selfobj->getShape(ndims, shape);
    return rv;
// splicer end class.View.method.get_shape
}

bool SIDRE_view_rename(SIDRE_view * self, const char * new_name)
{
View *selfobj = static_cast<View *>(static_cast<void *>(self));
// splicer begin class.View.method.rename
    const std::string SH_new_name(new_name);
    bool rv = selfobj->rename(SH_new_name);
    return rv;
// splicer end class.View.method.rename
}

bool SIDRE_view_rename_bufferify(SIDRE_view * self, const char * new_name, int Lnew_name)
{
View *selfobj = static_cast<View *>(static_cast<void *>(self));
// splicer begin class.View.method.rename_bufferify
    const std::string SH_new_name(new_name, Lnew_name);
    bool rv = selfobj->rename(SH_new_name);
    return rv;
// splicer end class.View.method.rename_bufferify
}

void SIDRE_view_print(const SIDRE_view * self)
{
const View *selfobj = static_cast<const View *>(static_cast<const void *>(self));
// splicer begin class.View.method.print
    selfobj->print();
    return;
// splicer end class.View.method.print
}

// splicer begin class.View.additional_functions
// splicer end class.View.additional_functions

}  // namespace axom
}  // namespace sidre
}  // extern "C"

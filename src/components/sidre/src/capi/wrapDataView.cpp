// wrapDataView.cpp
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
// wrapDataView.cpp
#include "wrapDataView.h"
#include "sidre/DataView.hpp"
#include "sidre/SidreTypes.hpp"

extern "C" {
namespace asctoolkit {
namespace sidre {

void SIDRE_dataview_allocate_simple(SIDRE_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.allocate_simple
selfobj->allocate();
return;
// splicer end class.DataView.method.allocate_simple
}

void SIDRE_dataview_allocate_from_type(SIDRE_dataview * self, int type, SIDRE_SidreLength num_elems)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.allocate_from_type
selfobj->allocate(getTypeID(type), num_elems);
return;
// splicer end class.DataView.method.allocate_from_type
}

void SIDRE_dataview_reallocate(SIDRE_dataview * self, SIDRE_SidreLength num_elems)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.reallocate
selfobj->reallocate(num_elems);
return;
// splicer end class.DataView.method.reallocate
}

void SIDRE_dataview_apply_0(SIDRE_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.apply_0
selfobj->apply();
return;
// splicer end class.DataView.method.apply_0
}

void SIDRE_dataview_attach_buffer(SIDRE_dataview * self, SIDRE_databuffer * buff)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.attach_buffer
selfobj->attachBuffer(static_cast<DataBuffer *>(static_cast<void *>(buff)));
return;
// splicer end class.DataView.method.attach_buffer
}

void SIDRE_dataview_apply_nelems(SIDRE_dataview * self, SIDRE_SidreLength num_elems)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.apply_nelems
selfobj->apply(num_elems);
return;
// splicer end class.DataView.method.apply_nelems
}

void SIDRE_dataview_apply_nelems_offset(SIDRE_dataview * self, SIDRE_SidreLength num_elems, SIDRE_SidreLength offset)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.apply_nelems_offset
selfobj->apply(num_elems, offset);
return;
// splicer end class.DataView.method.apply_nelems_offset
}

void SIDRE_dataview_apply_nelems_offset_stride(SIDRE_dataview * self, SIDRE_SidreLength num_elems, SIDRE_SidreLength offset, SIDRE_SidreLength stride)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.apply_nelems_offset_stride
selfobj->apply(num_elems, offset, stride);
return;
// splicer end class.DataView.method.apply_nelems_offset_stride
}

void SIDRE_dataview_apply_type_nelems(SIDRE_dataview * self, int type, SIDRE_SidreLength num_elems)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.apply_type_nelems
selfobj->apply(getTypeID(type), num_elems);
return;
// splicer end class.DataView.method.apply_type_nelems
}

void SIDRE_dataview_apply_type_nelems_offset(SIDRE_dataview * self, int type, SIDRE_SidreLength num_elems, SIDRE_SidreLength offset)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.apply_type_nelems_offset
selfobj->apply(getTypeID(type), num_elems, offset);
return;
// splicer end class.DataView.method.apply_type_nelems_offset
}

void SIDRE_dataview_apply_type_nelems_offset_stride(SIDRE_dataview * self, int type, SIDRE_SidreLength num_elems, SIDRE_SidreLength offset, SIDRE_SidreLength stride)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.apply_type_nelems_offset_stride
selfobj->apply(getTypeID(type), num_elems, offset, stride);
return;
// splicer end class.DataView.method.apply_type_nelems_offset_stride
}

bool SIDRE_dataview_has_buffer(SIDRE_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.has_buffer
bool rv = selfobj->hasBuffer();
return rv;
// splicer end class.DataView.method.has_buffer
}

bool SIDRE_dataview_is_opaque(SIDRE_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.is_opaque
bool rv = selfobj->isOpaque();
return rv;
// splicer end class.DataView.method.is_opaque
}

const char * SIDRE_dataview_get_name(const SIDRE_dataview * self)
{
const DataView *selfobj = static_cast<const DataView *>(static_cast<const void *>(self));
// splicer begin class.DataView.method.get_name
const std::string & rv = selfobj->getName();
return rv.c_str();
// splicer end class.DataView.method.get_name
}

void * SIDRE_dataview_get_opaque(SIDRE_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.get_opaque
void * rv = selfobj->getOpaque();
return rv;
// splicer end class.DataView.method.get_opaque
}

SIDRE_databuffer * SIDRE_dataview_get_buffer(SIDRE_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.get_buffer
DataBuffer * rv = selfobj->getBuffer();
return static_cast<SIDRE_databuffer *>(static_cast<void *>(rv));
// splicer end class.DataView.method.get_buffer
}

void * SIDRE_dataview_get_data_pointer(SIDRE_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.get_data_pointer
void * rv = selfobj->getDataPointer();
return rv;
// splicer end class.DataView.method.get_data_pointer
}

void SIDRE_dataview_set_value_int(SIDRE_dataview * self, int value)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.set_value_int
selfobj->setValue<int>(value);
return;
// splicer end class.DataView.method.set_value_int
}

void SIDRE_dataview_set_value_long(SIDRE_dataview * self, long value)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.set_value_long
selfobj->setValue<long>(value);
return;
// splicer end class.DataView.method.set_value_long
}

void SIDRE_dataview_set_value_float(SIDRE_dataview * self, float value)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.set_value_float
selfobj->setValue<float>(value);
return;
// splicer end class.DataView.method.set_value_float
}

void SIDRE_dataview_set_value_double(SIDRE_dataview * self, double value)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.set_value_double
selfobj->setValue<double>(value);
return;
// splicer end class.DataView.method.set_value_double
}

int SIDRE_dataview_get_value_int(SIDRE_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.get_value_int
int rv = selfobj->getValue<int>();
return rv;
// splicer end class.DataView.method.get_value_int
}

long SIDRE_dataview_get_value_long(SIDRE_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.get_value_long
long rv = selfobj->getValue<long>();
return rv;
// splicer end class.DataView.method.get_value_long
}

float SIDRE_dataview_get_value_float(SIDRE_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.get_value_float
float rv = selfobj->getValue<float>();
return rv;
// splicer end class.DataView.method.get_value_float
}

double SIDRE_dataview_get_value_double(SIDRE_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.get_value_double
double rv = selfobj->getValue<double>();
return rv;
// splicer end class.DataView.method.get_value_double
}

SIDRE_dataview * SIDRE_dataview_set_opaque(SIDRE_dataview * self, void * opaque_ptr)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.set_opaque
DataView * rv = selfobj->setOpaque(opaque_ptr);
return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataView.method.set_opaque
}

SIDRE_datagroup * SIDRE_dataview_get_owning_group(SIDRE_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.get_owning_group
DataGroup * rv = selfobj->getOwningGroup();
return static_cast<SIDRE_datagroup *>(static_cast<void *>(rv));
// splicer end class.DataView.method.get_owning_group
}

int SIDRE_dataview_get_type_id(SIDRE_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.get_type_id
TypeID rv = selfobj->getTypeID();
return static_cast<int>(rv);
// splicer end class.DataView.method.get_type_id
}

size_t SIDRE_dataview_get_total_bytes(SIDRE_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.get_total_bytes
size_t rv = selfobj->getTotalBytes();
return rv;
// splicer end class.DataView.method.get_total_bytes
}

size_t SIDRE_dataview_get_num_elements(SIDRE_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.get_num_elements
size_t rv = selfobj->getNumElements();
return rv;
// splicer end class.DataView.method.get_num_elements
}

int SIDRE_dataview_get_num_dimensions(SIDRE_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.get_num_dimensions
int rv = selfobj->getNumDimensions();
return rv;
// splicer end class.DataView.method.get_num_dimensions
}

int SIDRE_dataview_get_shape(SIDRE_dataview * self, int ndims, SIDRE_SidreLength * shape)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.get_shape
int rv = selfobj->getShape(ndims, shape);
return rv;
// splicer end class.DataView.method.get_shape
}

void SIDRE_dataview_print(SIDRE_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.print
selfobj->print();
return;
// splicer end class.DataView.method.print
}

// splicer begin class.DataView.additional_functions
// splicer end class.DataView.additional_functions

}  // namespace asctoolkit
}  // namespace sidre
}  // extern "C"

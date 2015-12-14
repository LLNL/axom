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
#include "sidre/SidreWrapperHelpers.hpp"

extern "C" {
namespace asctoolkit {
namespace sidre {

void ATK_dataview_allocate_simple(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.allocate_simple
selfobj->allocate();
return;
// splicer end class.DataView.method.allocate_simple
}

void ATK_dataview_allocate_from_type(ATK_dataview * self, int type, ATK_SidreLength numelems)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.allocate_from_type
selfobj->allocate(getTypeID(type), numelems);
return;
// splicer end class.DataView.method.allocate_from_type
}

void ATK_dataview_reallocate(ATK_dataview * self, ATK_SidreLength numelems)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.reallocate
selfobj->reallocate(numelems);
return;
// splicer end class.DataView.method.reallocate
}

ATK_dataview * ATK_dataview_apply_simple(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.apply_simple
DataView * rv = selfobj->apply();
return static_cast<ATK_dataview *>(static_cast<void *>(rv));
// splicer end class.DataView.method.apply_simple
}

ATK_dataview * ATK_dataview_apply_nelems(ATK_dataview * self, ATK_SidreLength numelems)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.apply_nelems
DataView * rv = selfobj->apply(numelems);
return static_cast<ATK_dataview *>(static_cast<void *>(rv));
// splicer end class.DataView.method.apply_nelems
}

ATK_dataview * ATK_dataview_apply_nelems_offset(ATK_dataview * self, ATK_SidreLength numelems, ATK_SidreLength offset)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.apply_nelems_offset
DataView * rv = selfobj->apply(numelems, offset);
return static_cast<ATK_dataview *>(static_cast<void *>(rv));
// splicer end class.DataView.method.apply_nelems_offset
}

ATK_dataview * ATK_dataview_apply_nelems_offset_stride(ATK_dataview * self, ATK_SidreLength numelems, ATK_SidreLength offset, ATK_SidreLength stride)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.apply_nelems_offset_stride
DataView * rv = selfobj->apply(numelems, offset, stride);
return static_cast<ATK_dataview *>(static_cast<void *>(rv));
// splicer end class.DataView.method.apply_nelems_offset_stride
}

ATK_dataview * ATK_dataview_apply_type_nelems(ATK_dataview * self, int type, ATK_SidreLength numelems)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.apply_type_nelems
DataView * rv = selfobj->apply(getTypeID(type), numelems);
return static_cast<ATK_dataview *>(static_cast<void *>(rv));
// splicer end class.DataView.method.apply_type_nelems
}

ATK_dataview * ATK_dataview_apply_type_nelems_offset(ATK_dataview * self, int type, ATK_SidreLength numelems, ATK_SidreLength offset)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.apply_type_nelems_offset
DataView * rv = selfobj->apply(getTypeID(type), numelems, offset);
return static_cast<ATK_dataview *>(static_cast<void *>(rv));
// splicer end class.DataView.method.apply_type_nelems_offset
}

ATK_dataview * ATK_dataview_apply_type_nelems_offset_stride(ATK_dataview * self, int type, ATK_SidreLength numelems, ATK_SidreLength offset, ATK_SidreLength stride)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.apply_type_nelems_offset_stride
DataView * rv = selfobj->apply(getTypeID(type), numelems, offset, stride);
return static_cast<ATK_dataview *>(static_cast<void *>(rv));
// splicer end class.DataView.method.apply_type_nelems_offset_stride
}

bool ATK_dataview_has_buffer(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.has_buffer
bool rv = selfobj->hasBuffer();
return rv;
// splicer end class.DataView.method.has_buffer
}

bool ATK_dataview_is_opaque(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.is_opaque
bool rv = selfobj->isOpaque();
return rv;
// splicer end class.DataView.method.is_opaque
}

const char * ATK_dataview_get_name(const ATK_dataview * self)
{
const DataView *selfobj = static_cast<const DataView *>(static_cast<const void *>(self));
// splicer begin class.DataView.method.get_name
const std::string & rv = selfobj->getName();
return rv.c_str();
// splicer end class.DataView.method.get_name
}

ATK_databuffer * ATK_dataview_get_buffer(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.get_buffer
DataBuffer * rv = selfobj->getBuffer();
return static_cast<ATK_databuffer *>(static_cast<void *>(rv));
// splicer end class.DataView.method.get_buffer
}

void * ATK_dataview_get_void_ptr(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.get_void_ptr
void * rv = selfobj->getVoidPtr();
return rv;
// splicer end class.DataView.method.get_void_ptr
}

void ATK_dataview_set_scalar_int(ATK_dataview * self, int value)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.set_scalar_int
selfobj->setScalar<int>(value);
return;
// splicer end class.DataView.method.set_scalar_int
}

void ATK_dataview_set_scalar_long(ATK_dataview * self, long value)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.set_scalar_long
selfobj->setScalar<long>(value);
return;
// splicer end class.DataView.method.set_scalar_long
}

void ATK_dataview_set_scalar_float(ATK_dataview * self, float value)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.set_scalar_float
selfobj->setScalar<float>(value);
return;
// splicer end class.DataView.method.set_scalar_float
}

void ATK_dataview_set_scalar_double(ATK_dataview * self, double value)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.set_scalar_double
selfobj->setScalar<double>(value);
return;
// splicer end class.DataView.method.set_scalar_double
}

int ATK_dataview_get_data_int(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.get_data_int
int rv = selfobj->getData<int>();
return rv;
// splicer end class.DataView.method.get_data_int
}

long ATK_dataview_get_data_long(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.get_data_long
long rv = selfobj->getData<long>();
return rv;
// splicer end class.DataView.method.get_data_long
}

float ATK_dataview_get_data_float(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.get_data_float
float rv = selfobj->getData<float>();
return rv;
// splicer end class.DataView.method.get_data_float
}

double ATK_dataview_get_data_double(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.get_data_double
double rv = selfobj->getData<double>();
return rv;
// splicer end class.DataView.method.get_data_double
}

ATK_datagroup * ATK_dataview_get_owning_group(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.get_owning_group
DataGroup * rv = selfobj->getOwningGroup();
return static_cast<ATK_datagroup *>(static_cast<void *>(rv));
// splicer end class.DataView.method.get_owning_group
}

int ATK_dataview_get_type_id(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.get_type_id
TypeID rv = selfobj->getTypeID();
return rv;
// splicer end class.DataView.method.get_type_id
}

size_t ATK_dataview_get_total_bytes(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.get_total_bytes
size_t rv = selfobj->getTotalBytes();
return rv;
// splicer end class.DataView.method.get_total_bytes
}

size_t ATK_dataview_get_num_elements(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.get_num_elements
size_t rv = selfobj->getNumElements();
return rv;
// splicer end class.DataView.method.get_num_elements
}

void ATK_dataview_print(ATK_dataview * self)
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

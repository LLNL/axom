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
#define EXAMPLE_WRAPPER_IMPL
#include "wrapDataView.h"
#include "sidre/DataView.hpp"
#include "sidre/SidreWrapperHelpers.hpp"

extern "C" {
namespace asctoolkit {
namespace sidre {

void ATK_dataview_declare(ATK_dataview * self, int type, ATK_SidreLength len)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin class.DataView.method.declare
selfobj->declare(getTypeID(type), len);
return;
// splicer end class.DataView.method.declare
}

void ATK_dataview_allocate_simple(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin class.DataView.method.allocate_simple
selfobj->allocate();
return;
// splicer end class.DataView.method.allocate_simple
}

void ATK_dataview_allocate_from_type(ATK_dataview * self, int type, ATK_SidreLength len)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin class.DataView.method.allocate_from_type
selfobj->allocate(getTypeID(type), len);
return;
// splicer end class.DataView.method.allocate_from_type
}

void ATK_dataview_reallocate(ATK_dataview * self, int type, ATK_SidreLength len)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin class.DataView.method.reallocate
selfobj->reallocate(getTypeID(type), len);
return;
// splicer end class.DataView.method.reallocate
}

bool ATK_dataview_has_buffer(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin class.DataView.method.has_buffer
bool rv = selfobj->hasBuffer();
return rv;
// splicer end class.DataView.method.has_buffer
}

bool ATK_dataview_is_opaque(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin class.DataView.method.is_opaque
bool rv = selfobj->isOpaque();
return rv;
// splicer end class.DataView.method.is_opaque
}

const char * ATK_dataview_get_name(const ATK_dataview * self)
{
const DataView *selfobj = static_cast<const DataView *>(self);
// splicer begin class.DataView.method.get_name
const std::string & rv = selfobj->getName();
return rv.c_str();
// splicer end class.DataView.method.get_name
}

void * ATK_dataview_get_opaque(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin class.DataView.method.get_opaque
void * rv = selfobj->getOpaque();
return rv;
// splicer end class.DataView.method.get_opaque
}

ATK_databuffer * ATK_dataview_get_buffer(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin class.DataView.method.get_buffer
DataBuffer * rv = selfobj->getBuffer();
return rv;
// splicer end class.DataView.method.get_buffer
}

void * ATK_dataview_get_data_pointer(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin class.DataView.method.get_data_pointer
void * rv = selfobj->getDataPointer();
return rv;
// splicer end class.DataView.method.get_data_pointer
}

ATK_datagroup * ATK_dataview_get_owning_group(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin class.DataView.method.get_owning_group
DataGroup * rv = selfobj->getOwningGroup();
return rv;
// splicer end class.DataView.method.get_owning_group
}

int ATK_dataview_get_type_id(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin class.DataView.method.get_type_id
TypeID rv = selfobj->getTypeID();
return rv;
// splicer end class.DataView.method.get_type_id
}

size_t ATK_dataview_get_total_bytes(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin class.DataView.method.get_total_bytes
size_t rv = selfobj->getTotalBytes();
return rv;
// splicer end class.DataView.method.get_total_bytes
}

size_t ATK_dataview_get_number_of_elements(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin class.DataView.method.get_number_of_elements
size_t rv = selfobj->getNumberOfElements();
return rv;
// splicer end class.DataView.method.get_number_of_elements
}

void ATK_dataview_print(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin class.DataView.method.print
selfobj->print();
return;
// splicer end class.DataView.method.print
}

void ATK_dataview_set_value_int(ATK_dataview * self, int value)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin class.DataView.method.set_value_int
selfobj->setValue(value);
return;
// splicer end class.DataView.method.set_value_int
}

void ATK_dataview_set_value_long(ATK_dataview * self, long value)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin class.DataView.method.set_value_long
selfobj->setValue(value);
return;
// splicer end class.DataView.method.set_value_long
}

void ATK_dataview_set_value_float(ATK_dataview * self, float value)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin class.DataView.method.set_value_float
selfobj->setValue(value);
return;
// splicer end class.DataView.method.set_value_float
}

void ATK_dataview_set_value_double(ATK_dataview * self, double value)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin class.DataView.method.set_value_double
selfobj->setValue(value);
return;
// splicer end class.DataView.method.set_value_double
}

int ATK_dataview_get_value_int(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin class.DataView.method.get_value_int
int rv = selfobj->getValue<int>();
return rv;
// splicer end class.DataView.method.get_value_int
}

long ATK_dataview_get_value_long(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin class.DataView.method.get_value_long
long rv = selfobj->getValue<long>();
return rv;
// splicer end class.DataView.method.get_value_long
}

float ATK_dataview_get_value_float(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin class.DataView.method.get_value_float
float rv = selfobj->getValue<float>();
return rv;
// splicer end class.DataView.method.get_value_float
}

double ATK_dataview_get_value_double(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin class.DataView.method.get_value_double
double rv = selfobj->getValue<double>();
return rv;
// splicer end class.DataView.method.get_value_double
}

// splicer begin class.DataView.additional_functions
// splicer end class.DataView.additional_functions

}  // namespace asctoolkit
}  // namespace sidre
}  // extern "C"

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

void ATK_dataview_allocate(ATK_dataview * self, int type, ATK_SidreLength len)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin class.DataView.method.allocate
selfobj->allocate(getTypeID(type), len);
return;
// splicer end class.DataView.method.allocate
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
// splicer begin class.DataView.method.hasBuffer
bool rv = selfobj->hasBuffer();
return rv;
// splicer end class.DataView.method.hasBuffer
}

bool ATK_dataview_is_opaque(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin class.DataView.method.isOpaque
bool rv = selfobj->isOpaque();
return rv;
// splicer end class.DataView.method.isOpaque
}

const char * ATK_dataview_get_name(const ATK_dataview * self)
{
const DataView *selfobj = static_cast<const DataView *>(self);
// splicer begin class.DataView.method.getName
const std::string & rv = selfobj->getName();
return isNameValid(rv) ? rv.c_str() : ATK_InvalidName;
// splicer end class.DataView.method.getName
}

void * ATK_dataview_get_opaque(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin class.DataView.method.getOpaque
void * rv = selfobj->getOpaque();
return rv;
// splicer end class.DataView.method.getOpaque
}

ATK_databuffer * ATK_dataview_get_buffer(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin class.DataView.method.getBuffer
DataBuffer * rv = selfobj->getBuffer();
return rv;
// splicer end class.DataView.method.getBuffer
}

void * ATK_dataview_get_data_pointer(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin class.DataView.method.getDataPointer
void * rv = selfobj->getDataPointer();
return rv;
// splicer end class.DataView.method.getDataPointer
}

ATK_datagroup * ATK_dataview_get_owning_group(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin class.DataView.method.getOwningGroup
DataGroup * rv = selfobj->getOwningGroup();
return rv;
// splicer end class.DataView.method.getOwningGroup
}

int ATK_dataview_get_type_id(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin class.DataView.method.getTypeID
TypeID rv = selfobj->getTypeID();
return rv;
// splicer end class.DataView.method.getTypeID
}

size_t ATK_dataview_get_total_bytes(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin class.DataView.method.getTotalBytes
size_t rv = selfobj->getTotalBytes();
return rv;
// splicer end class.DataView.method.getTotalBytes
}

size_t ATK_dataview_get_number_of_elements(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin class.DataView.method.getNumberOfElements
size_t rv = selfobj->getNumberOfElements();
return rv;
// splicer end class.DataView.method.getNumberOfElements
}

// splicer begin class.DataView.additional_functions
// splicer end class.DataView.additional_functions

}  // namespace asctoolkit
}  // namespace sidre
}  // extern "C"

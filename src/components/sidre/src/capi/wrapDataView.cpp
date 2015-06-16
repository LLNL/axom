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

extern "C" {
namespace asctoolkit {
namespace sidre {
// splicer push class.DataView.method

void ATK_dataview_declare(ATK_dataview * self, int type, ATK_SidreLength len)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin ATK_dataview_declare
selfobj->declare(getTypeID(type), len);
return;
// splicer end ATK_dataview_declare
}

void ATK_dataview_allocate(ATK_dataview * self, int type, ATK_SidreLength len)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin ATK_dataview_allocate
selfobj->allocate(getTypeID(type), len);
return;
// splicer end ATK_dataview_allocate
}

void ATK_dataview_reallocate(ATK_dataview * self, int type, ATK_SidreLength len)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin ATK_dataview_reallocate
selfobj->reallocate(getTypeID(type), len);
return;
// splicer end ATK_dataview_reallocate
}

bool ATK_dataview_has_buffer(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin ATK_dataview_has_buffer
bool rv = selfobj->hasBuffer();
return rv;
// splicer end ATK_dataview_has_buffer
}

bool ATK_dataview_is_opaque(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin ATK_dataview_is_opaque
bool rv = selfobj->isOpaque();
return rv;
// splicer end ATK_dataview_is_opaque
}

const char * ATK_dataview_get_name(const ATK_dataview * self)
{
const DataView *selfobj = static_cast<const DataView *>(self);
// splicer begin ATK_dataview_get_name
const std::string & rv = selfobj->getName();
return rv.c_str();
// splicer end ATK_dataview_get_name
}

void * ATK_dataview_get_opaque(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin ATK_dataview_get_opaque
void * rv = selfobj->getOpaque();
return rv;
// splicer end ATK_dataview_get_opaque
}

ATK_databuffer * ATK_dataview_get_buffer(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin ATK_dataview_get_buffer
DataBuffer * rv = selfobj->getBuffer();
return rv;
// splicer end ATK_dataview_get_buffer
}

void * ATK_dataview_get_data_pointer(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin ATK_dataview_get_data_pointer
void * rv = selfobj->getDataPointer();
return rv;
// splicer end ATK_dataview_get_data_pointer
}

ATK_datagroup * ATK_dataview_get_owning_group(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin ATK_dataview_get_owning_group
DataGroup * rv = selfobj->getOwningGroup();
return rv;
// splicer end ATK_dataview_get_owning_group
}

int ATK_dataview_get_type_id(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin ATK_dataview_get_type_id
TypeID rv = selfobj->getTypeID();
return rv;
// splicer end ATK_dataview_get_type_id
}

size_t ATK_dataview_get_total_bytes(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin ATK_dataview_get_total_bytes
size_t rv = selfobj->getTotalBytes();
return rv;
// splicer end ATK_dataview_get_total_bytes
}

size_t ATK_dataview_get_number_of_elements(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin ATK_dataview_get_number_of_elements
size_t rv = selfobj->getNumberOfElements();
return rv;
// splicer end ATK_dataview_get_number_of_elements
}

// splicer pop.class.DataView method

}  // namespace asctoolkit
}  // namespace sidre
}  // extern "C"

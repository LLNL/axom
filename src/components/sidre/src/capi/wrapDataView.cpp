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

void ATK_dataview_declare(ATK_dataview * self, int type, ATK_SidreLength len)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin
selfobj->declare(getTypeID(type), len);
return;
// splicer end
}

void ATK_dataview_allocate(ATK_dataview * self, int type, ATK_SidreLength len)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin
selfobj->allocate(getTypeID(type), len);
return;
// splicer end
}

void ATK_dataview_reallocate(ATK_dataview * self, int type, ATK_SidreLength len)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin
selfobj->reallocate(getTypeID(type), len);
return;
// splicer end
}

bool ATK_dataview_has_buffer(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin
bool rv = selfobj->hasBuffer();
return rv;
// splicer end
}

bool ATK_dataview_is_opaque(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin
bool rv = selfobj->isOpaque();
return rv;
// splicer end
}

const char * ATK_dataview_get_name(const ATK_dataview * self)
{
const DataView *selfobj = static_cast<const DataView *>(self);
// splicer begin
const std::string & rv = selfobj->getName();
return rv.c_str();
// splicer end
}

void * ATK_dataview_get_data_in_buffer(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin
void * rv = selfobj->getDataInBuffer();
return rv;
// splicer end
}

void * ATK_dataview_get_opaque(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin
void * rv = selfobj->getOpaque();
return rv;
// splicer end
}

ATK_databuffer * ATK_dataview_get_buffer(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin
DataBuffer * rv = selfobj->getBuffer();
return rv;
// splicer end
}

void * ATK_dataview_get_data(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin
void * rv = selfobj->getData();
return rv;
// splicer end
}

ATK_datagroup * ATK_dataview_get_owning_group(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin
DataGroup * rv = selfobj->getOwningGroup();
return rv;
// splicer end
}

size_t ATK_dataview_get_total_bytes(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin
size_t rv = selfobj->getTotalBytes();
return rv;
// splicer end
}

size_t ATK_dataview_get_number_of_elements(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin
size_t rv = selfobj->getNumberOfElements();
return rv;
// splicer end
}

}  // namespace asctoolkit
}  // namespace sidre
}  // extern "C"

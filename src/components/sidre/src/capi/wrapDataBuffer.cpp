// wrapDataBuffer.cpp
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
// wrapDataBuffer.cpp
#define EXAMPLE_WRAPPER_IMPL
#include "wrapDataBuffer.h"
#include "sidre/DataBuffer.hpp"
#include "sidre/SidreWrapperHelpers.hpp"

extern "C" {
namespace asctoolkit {
namespace sidre {

ATK_IndexType ATK_databuffer_get_index(ATK_databuffer * self)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(self);
// splicer begin class.DataBuffer.method.getIndex
IndexType rv = selfobj->getIndex();
return rv;
// splicer end class.DataBuffer.method.getIndex
}

size_t ATK_databuffer_get_num_views(ATK_databuffer * self)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(self);
// splicer begin class.DataBuffer.method.getNumViews
size_t rv = selfobj->getNumViews();
return rv;
// splicer end class.DataBuffer.method.getNumViews
}

void ATK_databuffer_declare(ATK_databuffer * self, int type, ATK_SidreLength len)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(self);
// splicer begin class.DataBuffer.method.declare
selfobj->declare(getTypeID(type), len);
return;
// splicer end class.DataBuffer.method.declare
}

void ATK_databuffer_declare_external(ATK_databuffer * self, void * external_data, int type, ATK_SidreLength len)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(self);
// splicer begin class.DataBuffer.method.declareExternal
selfobj->declareExternal(external_data, getTypeID(type), len);
return;
// splicer end class.DataBuffer.method.declareExternal
}

void ATK_databuffer_allocate_existing(ATK_databuffer * self)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(self);
// splicer begin class.DataBuffer.method.allocate
selfobj->allocate();
return;
// splicer end class.DataBuffer.method.allocate
}

void ATK_databuffer_allocate_from_type(ATK_databuffer * self, int type, ATK_SidreLength len)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(self);
// splicer begin class.DataBuffer.method.allocate
selfobj->allocate(getTypeID(type), len);
return;
// splicer end class.DataBuffer.method.allocate
}

void ATK_databuffer_reallocate(ATK_databuffer * self, int type, ATK_SidreLength len)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(self);
// splicer begin class.DataBuffer.method.reallocate
selfobj->reallocate(getTypeID(type), len);
return;
// splicer end class.DataBuffer.method.reallocate
}

bool ATK_databuffer_is_external(ATK_databuffer * self)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(self);
// splicer begin class.DataBuffer.method.isExternal
bool rv = selfobj->isExternal();
return rv;
// splicer end class.DataBuffer.method.isExternal
}

void * ATK_databuffer_get_data(ATK_databuffer * self)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(self);
// splicer begin class.DataBuffer.method.getData
void * rv = selfobj->getData();
return rv;
// splicer end class.DataBuffer.method.getData
}

size_t ATK_databuffer_get_total_bytes(ATK_databuffer * self)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(self);
// splicer begin class.DataBuffer.method.getTotalBytes
size_t rv = selfobj->getTotalBytes();
return rv;
// splicer end class.DataBuffer.method.getTotalBytes
}

// splicer begin class.DataBuffer.additional_functions
// splicer end class.DataBuffer.additional_functions

}  // namespace asctoolkit
}  // namespace sidre
}  // extern "C"

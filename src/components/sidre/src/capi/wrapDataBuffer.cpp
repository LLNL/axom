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
#include "wrapDataBuffer.h"
#include "sidre/DataBuffer.hpp"
#include "sidre/SidreWrapperHelpers.hpp"

extern "C" {
namespace asctoolkit {
namespace sidre {

ATK_IndexType ATK_databuffer_get_index(ATK_databuffer * self)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(static_cast<void *>(self));
// splicer begin class.DataBuffer.method.get_index
IndexType rv = selfobj->getIndex();
return rv;
// splicer end class.DataBuffer.method.get_index
}

size_t ATK_databuffer_get_num_views(ATK_databuffer * self)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(static_cast<void *>(self));
// splicer begin class.DataBuffer.method.get_num_views
size_t rv = selfobj->getNumViews();
return rv;
// splicer end class.DataBuffer.method.get_num_views
}

void ATK_databuffer_declare(ATK_databuffer * self, int type, ATK_SidreLength len)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(static_cast<void *>(self));
// splicer begin class.DataBuffer.method.declare
selfobj->declare(getTypeID(type), len);
return;
// splicer end class.DataBuffer.method.declare
}

void ATK_databuffer_allocate_existing(ATK_databuffer * self)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(static_cast<void *>(self));
// splicer begin class.DataBuffer.method.allocate_existing
selfobj->allocate();
return;
// splicer end class.DataBuffer.method.allocate_existing
}

void ATK_databuffer_allocate_from_type(ATK_databuffer * self, int type, ATK_SidreLength len)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(static_cast<void *>(self));
// splicer begin class.DataBuffer.method.allocate_from_type
selfobj->allocate(getTypeID(type), len);
return;
// splicer end class.DataBuffer.method.allocate_from_type
}

void ATK_databuffer_reallocate(ATK_databuffer * self, ATK_SidreLength len)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(static_cast<void *>(self));
// splicer begin class.DataBuffer.method.reallocate
selfobj->reallocate(len);
return;
// splicer end class.DataBuffer.method.reallocate
}

void ATK_databuffer_set_external_data(ATK_databuffer * self, void * external_data)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(static_cast<void *>(self));
// splicer begin class.DataBuffer.method.set_external_data
selfobj->setExternalData(external_data);
return;
// splicer end class.DataBuffer.method.set_external_data
}

bool ATK_databuffer_is_external(ATK_databuffer * self)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(static_cast<void *>(self));
// splicer begin class.DataBuffer.method.is_external
bool rv = selfobj->isExternal();
return rv;
// splicer end class.DataBuffer.method.is_external
}

void * ATK_databuffer_get_data(ATK_databuffer * self)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(static_cast<void *>(self));
// splicer begin class.DataBuffer.method.get_data
void * rv = selfobj->getData();
return rv;
// splicer end class.DataBuffer.method.get_data
}

int ATK_databuffer_get_type_id(ATK_databuffer * self)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(static_cast<void *>(self));
// splicer begin class.DataBuffer.method.get_type_id
TypeID rv = selfobj->getTypeID();
return rv;
// splicer end class.DataBuffer.method.get_type_id
}

size_t ATK_databuffer_get_num_elements(ATK_databuffer * self)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(static_cast<void *>(self));
// splicer begin class.DataBuffer.method.get_num_elements
size_t rv = selfobj->getNumElements();
return rv;
// splicer end class.DataBuffer.method.get_num_elements
}

size_t ATK_databuffer_get_total_bytes(ATK_databuffer * self)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(static_cast<void *>(self));
// splicer begin class.DataBuffer.method.get_total_bytes
size_t rv = selfobj->getTotalBytes();
return rv;
// splicer end class.DataBuffer.method.get_total_bytes
}

void ATK_databuffer_print(ATK_databuffer * self)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(static_cast<void *>(self));
// splicer begin class.DataBuffer.method.print
selfobj->print();
return;
// splicer end class.DataBuffer.method.print
}

// splicer begin class.DataBuffer.additional_functions
// splicer end class.DataBuffer.additional_functions

}  // namespace asctoolkit
}  // namespace sidre
}  // extern "C"

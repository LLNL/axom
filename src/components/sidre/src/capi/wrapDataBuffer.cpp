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
// splicer push class.DataBuffer.method

ATK_IndexType ATK_databuffer_get_index(ATK_databuffer * self)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(self);
// splicer begin ATK_databuffer_get_index
IndexType rv = selfobj->getIndex();
return rv;
// splicer end ATK_databuffer_get_index
}

size_t ATK_databuffer_get_num_views(ATK_databuffer * self)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(self);
// splicer begin ATK_databuffer_get_num_views
size_t rv = selfobj->getNumViews();
return rv;
// splicer end ATK_databuffer_get_num_views
}

void ATK_databuffer_declare(ATK_databuffer * self, int type, ATK_SidreLength len)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(self);
// splicer begin ATK_databuffer_declare
selfobj->declare(getTypeID(type), len);
return;
// splicer end ATK_databuffer_declare
}

void ATK_databuffer_declare_external(ATK_databuffer * self, void * external_data, int type, ATK_SidreLength len)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(self);
// splicer begin ATK_databuffer_declare_external
selfobj->declareExternal(external_data, getTypeID(type), len);
return;
// splicer end ATK_databuffer_declare_external
}

void ATK_databuffer_allocate_existing(ATK_databuffer * self)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(self);
// splicer begin ATK_databuffer_allocate_existing
selfobj->allocate();
return;
// splicer end ATK_databuffer_allocate_existing
}

void ATK_databuffer_allocate_from_type(ATK_databuffer * self, int type, ATK_SidreLength len)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(self);
// splicer begin ATK_databuffer_allocate_from_type
selfobj->allocate(getTypeID(type), len);
return;
// splicer end ATK_databuffer_allocate_from_type
}

void ATK_databuffer_reallocate(ATK_databuffer * self, int type, ATK_SidreLength len)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(self);
// splicer begin ATK_databuffer_reallocate
selfobj->reallocate(getTypeID(type), len);
return;
// splicer end ATK_databuffer_reallocate
}

bool ATK_databuffer_is_external(ATK_databuffer * self)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(self);
// splicer begin ATK_databuffer_is_external
bool rv = selfobj->isExternal();
return rv;
// splicer end ATK_databuffer_is_external
}

void * ATK_databuffer_get_data(ATK_databuffer * self)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(self);
// splicer begin ATK_databuffer_get_data
void * rv = selfobj->getData();
return rv;
// splicer end ATK_databuffer_get_data
}

size_t ATK_databuffer_get_total_bytes(ATK_databuffer * self)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(self);
// splicer begin ATK_databuffer_get_total_bytes
size_t rv = selfobj->getTotalBytes();
return rv;
// splicer end ATK_databuffer_get_total_bytes
}

// splicer pop.class.DataBuffer method

}  // namespace asctoolkit
}  // namespace sidre
}  // extern "C"

//
// Copyright (c) 2015, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory.
//
// All rights reserved.
//
// This source code cannot be distributed without permission and
// further review from Lawrence Livermore National Laboratory.
//
// wrapDataStore.cpp
#define EXAMPLE_WRAPPER_IMPL
#include "wrapDataStore.h"
#include "sidre/DataStore.hpp"

extern "C" {
namespace asctoolkit {
namespace sidre {

ATK_datastore * ATK_datastore_new()
{
DataStore *selfobj = new DataStore();
// splicer begin
return (ATK_datastore *) selfobj;
// splicer end
}

void ATK_datastore_delete(ATK_datastore * self)
{
DataStore *selfobj = static_cast<DataStore *>(self);
// splicer begin
delete selfobj;
// splicer end
}

ATK_databuffer * ATK_datastore_create_buffer(ATK_datastore * self)
{
DataStore *selfobj = static_cast<DataStore *>(self);
// splicer begin
DataBuffer * rv = selfobj->createBuffer();
return rv;
// splicer end
}

void ATK_datastore_destroy_buffer(ATK_datastore * self, const ATK_IndexType id)
{
DataStore *selfobj = static_cast<DataStore *>(self);
// splicer begin
selfobj->destroyBuffer(id);
return;
// splicer end
}

ATK_datagroup * ATK_datastore_get_root(ATK_datastore * self)
{
DataStore *selfobj = static_cast<DataStore *>(self);
// splicer begin
DataGroup * rv = selfobj->getRoot();
return rv;
// splicer end
}

ATK_databuffer * ATK_datastore_get_buffer(ATK_datastore * self, const ATK_IndexType id)
{
DataStore *selfobj = static_cast<DataStore *>(self);
// splicer begin
DataBuffer * rv = selfobj->getBuffer(id);
return rv;
// splicer end
}

size_t ATK_datastore_get_num_buffers(ATK_datastore * self)
{
DataStore *selfobj = static_cast<DataStore *>(self);
// splicer begin
size_t rv = selfobj->getNumBuffers();
return rv;
// splicer end
}

}  // namespace asctoolkit
}  // namespace sidre
}  // extern "C"

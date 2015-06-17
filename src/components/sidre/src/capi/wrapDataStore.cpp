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
// splicer push class.DataStore.method

ATK_datastore * ATK_datastore_new()
{
DataStore *selfobj = new DataStore();
// splicer begin ATK_datastore_new
return (ATK_datastore *) selfobj;
// splicer end ATK_datastore_new
}

void ATK_datastore_delete(ATK_datastore * self)
{
DataStore *selfobj = static_cast<DataStore *>(self);
// splicer begin ATK_datastore_delete
delete selfobj;
// splicer end ATK_datastore_delete
}

ATK_datagroup * ATK_datastore_get_root(ATK_datastore * self)
{
DataStore *selfobj = static_cast<DataStore *>(self);
// splicer begin ATK_datastore_get_root
DataGroup * rv = selfobj->getRoot();
return rv;
// splicer end ATK_datastore_get_root
}

ATK_databuffer * ATK_datastore_get_buffer(ATK_datastore * self, ATK_IndexType idx)
{
DataStore *selfobj = static_cast<DataStore *>(self);
// splicer begin ATK_datastore_get_buffer
DataBuffer * rv = selfobj->getBuffer(idx);
return rv;
// splicer end ATK_datastore_get_buffer
}

ATK_databuffer * ATK_datastore_create_buffer(ATK_datastore * self)
{
DataStore *selfobj = static_cast<DataStore *>(self);
// splicer begin ATK_datastore_create_buffer
DataBuffer * rv = selfobj->createBuffer();
return rv;
// splicer end ATK_datastore_create_buffer
}

void ATK_datastore_destroy_buffer(ATK_datastore * self, ATK_IndexType id)
{
DataStore *selfobj = static_cast<DataStore *>(self);
// splicer begin ATK_datastore_destroy_buffer
selfobj->destroyBuffer(id);
return;
// splicer end ATK_datastore_destroy_buffer
}

size_t ATK_datastore_get_num_buffers(ATK_datastore * self)
{
DataStore *selfobj = static_cast<DataStore *>(self);
// splicer begin ATK_datastore_get_num_buffers
size_t rv = selfobj->getNumBuffers();
return rv;
// splicer end ATK_datastore_get_num_buffers
}

void ATK_datastore_print(ATK_datastore * self)
{
DataStore *selfobj = static_cast<DataStore *>(self);
// splicer begin ATK_datastore_print
selfobj->print();
return;
// splicer end ATK_datastore_print
}

// splicer pop.class.DataStore method

}  // namespace asctoolkit
}  // namespace sidre
}  // extern "C"

// wrapDataStore.cpp
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
// splicer begin class.DataStore.method.new
return (ATK_datastore *) selfobj;
// splicer end class.DataStore.method.new
}

void ATK_datastore_delete(ATK_datastore * self)
{
DataStore *selfobj = static_cast<DataStore *>(self);
// splicer begin class.DataStore.method.delete
delete selfobj;
// splicer end class.DataStore.method.delete
}

ATK_datagroup * ATK_datastore_get_root(ATK_datastore * self)
{
DataStore *selfobj = static_cast<DataStore *>(self);
// splicer begin class.DataStore.method.getRoot
DataGroup * rv = selfobj->getRoot();
return rv;
// splicer end class.DataStore.method.getRoot
}

ATK_databuffer * ATK_datastore_get_buffer(ATK_datastore * self, ATK_IndexType idx)
{
DataStore *selfobj = static_cast<DataStore *>(self);
// splicer begin class.DataStore.method.getBuffer
DataBuffer * rv = selfobj->getBuffer(idx);
return rv;
// splicer end class.DataStore.method.getBuffer
}

ATK_databuffer * ATK_datastore_create_buffer(ATK_datastore * self)
{
DataStore *selfobj = static_cast<DataStore *>(self);
// splicer begin class.DataStore.method.createBuffer
DataBuffer * rv = selfobj->createBuffer();
return rv;
// splicer end class.DataStore.method.createBuffer
}

void ATK_datastore_destroy_buffer(ATK_datastore * self, ATK_IndexType id)
{
DataStore *selfobj = static_cast<DataStore *>(self);
// splicer begin class.DataStore.method.destroyBuffer
selfobj->destroyBuffer(id);
return;
// splicer end class.DataStore.method.destroyBuffer
}

size_t ATK_datastore_get_num_buffers(ATK_datastore * self)
{
DataStore *selfobj = static_cast<DataStore *>(self);
// splicer begin class.DataStore.method.getNumBuffers
size_t rv = selfobj->getNumBuffers();
return rv;
// splicer end class.DataStore.method.getNumBuffers
}

void ATK_datastore_print(ATK_datastore * self)
{
DataStore *selfobj = static_cast<DataStore *>(self);
// splicer begin class.DataStore.method.print
selfobj->print();
return;
// splicer end class.DataStore.method.print
}

// splicer begin class.DataStore.additional_functions
// splicer end class.DataStore.additional_functions

}  // namespace asctoolkit
}  // namespace sidre
}  // extern "C"

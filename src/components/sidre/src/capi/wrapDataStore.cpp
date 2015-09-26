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
#include "wrapDataStore.h"
#include "sidre/DataStore.hpp"

extern "C" {
namespace asctoolkit {
namespace sidre {

ATK_datastore * ATK_datastore_new()
{
DataStore *selfobj = new DataStore();
// splicer begin class.DataStore.method.new
return static_cast<ATK_datastore *>(static_cast<void *>(selfobj));
// splicer end class.DataStore.method.new
}

void ATK_datastore_delete(ATK_datastore * self)
{
DataStore *selfobj = static_cast<DataStore *>(static_cast<void *>(self));
// splicer begin class.DataStore.method.delete
delete selfobj;
// splicer end class.DataStore.method.delete
}

ATK_datagroup * ATK_datastore_get_root(ATK_datastore * self)
{
DataStore *selfobj = static_cast<DataStore *>(static_cast<void *>(self));
// splicer begin class.DataStore.method.get_root
DataGroup * rv = selfobj->getRoot();
return static_cast<ATK_datagroup *>(static_cast<void *>(rv));
// splicer end class.DataStore.method.get_root
}

ATK_databuffer * ATK_datastore_get_buffer(ATK_datastore * self, ATK_IndexType idx)
{
DataStore *selfobj = static_cast<DataStore *>(static_cast<void *>(self));
// splicer begin class.DataStore.method.get_buffer
DataBuffer * rv = selfobj->getBuffer(idx);
return static_cast<ATK_databuffer *>(static_cast<void *>(rv));
// splicer end class.DataStore.method.get_buffer
}

ATK_databuffer * ATK_datastore_create_buffer(ATK_datastore * self)
{
DataStore *selfobj = static_cast<DataStore *>(static_cast<void *>(self));
// splicer begin class.DataStore.method.create_buffer
DataBuffer * rv = selfobj->createBuffer();
return static_cast<ATK_databuffer *>(static_cast<void *>(rv));
// splicer end class.DataStore.method.create_buffer
}

void ATK_datastore_destroy_buffer(ATK_datastore * self, ATK_IndexType id)
{
DataStore *selfobj = static_cast<DataStore *>(static_cast<void *>(self));
// splicer begin class.DataStore.method.destroy_buffer
selfobj->destroyBuffer(id);
return;
// splicer end class.DataStore.method.destroy_buffer
}

size_t ATK_datastore_get_num_buffers(ATK_datastore * self)
{
DataStore *selfobj = static_cast<DataStore *>(static_cast<void *>(self));
// splicer begin class.DataStore.method.get_num_buffers
size_t rv = selfobj->getNumBuffers();
return rv;
// splicer end class.DataStore.method.get_num_buffers
}

void ATK_datastore_print(ATK_datastore * self)
{
DataStore *selfobj = static_cast<DataStore *>(static_cast<void *>(self));
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

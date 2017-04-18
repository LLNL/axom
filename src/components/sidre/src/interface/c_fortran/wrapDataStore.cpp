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
#include "sidre/SidreTypes.hpp"

extern "C" {
namespace axom
{
namespace sidre
{

SIDRE_datastore * SIDRE_datastore_new()
{
// splicer begin class.DataStore.method.new
  DataStore * SH_rv = new DataStore();
  return static_cast<SIDRE_datastore *>(static_cast<void *>(SH_rv));
// splicer end class.DataStore.method.new
}

void SIDRE_datastore_delete(SIDRE_datastore * self)
{
// splicer begin class.DataStore.method.delete
  DataStore * SH_this = static_cast<DataStore *>(static_cast<void *>(self));
  delete SH_this;
  return;
// splicer end class.DataStore.method.delete
}

SIDRE_group * SIDRE_datastore_get_root(SIDRE_datastore * self)
{
// splicer begin class.DataStore.method.get_root
  DataStore * SH_this = static_cast<DataStore *>(static_cast<void *>(self));
  Group * SH_rv = SH_this->getRoot();
  return static_cast<SIDRE_group *>(static_cast<void *>(SH_rv));
// splicer end class.DataStore.method.get_root
}

SIDRE_buffer * SIDRE_datastore_get_buffer(SIDRE_datastore * self,
                                          SIDRE_IndexType idx)
{
// splicer begin class.DataStore.method.get_buffer
  DataStore * SH_this = static_cast<DataStore *>(static_cast<void *>(self));
  Buffer * SH_rv = SH_this->getBuffer(idx);
  return static_cast<SIDRE_buffer *>(static_cast<void *>(SH_rv));
// splicer end class.DataStore.method.get_buffer
}

SIDRE_buffer * SIDRE_datastore_create_buffer_empty(SIDRE_datastore * self)
{
// splicer begin class.DataStore.method.create_buffer_empty
  DataStore * SH_this = static_cast<DataStore *>(static_cast<void *>(self));
  Buffer * SH_rv = SH_this->createBuffer();
  return static_cast<SIDRE_buffer *>(static_cast<void *>(SH_rv));
// splicer end class.DataStore.method.create_buffer_empty
}

SIDRE_buffer * SIDRE_datastore_create_buffer_from_type(SIDRE_datastore * self,
                                                       int type,
                                                       SIDRE_SidreLength num_elems)
{
// splicer begin class.DataStore.method.create_buffer_from_type
  DataStore * SH_this = static_cast<DataStore *>(static_cast<void *>(self));
  Buffer * SH_rv = SH_this->createBuffer(getTypeID(type), num_elems);
  return static_cast<SIDRE_buffer *>(static_cast<void *>(SH_rv));
// splicer end class.DataStore.method.create_buffer_from_type
}

void SIDRE_datastore_destroy_buffer(SIDRE_datastore * self, SIDRE_IndexType id)
{
// splicer begin class.DataStore.method.destroy_buffer
  DataStore * SH_this = static_cast<DataStore *>(static_cast<void *>(self));
  SH_this->destroyBuffer(id);
  return;
// splicer end class.DataStore.method.destroy_buffer
}

size_t SIDRE_datastore_get_num_buffers(const SIDRE_datastore * self)
{
// splicer begin class.DataStore.method.get_num_buffers
  const DataStore * SH_this =
    static_cast<const DataStore *>(static_cast<const void *>(self));
  size_t SH_rv = SH_this->getNumBuffers();
  return SH_rv;
// splicer end class.DataStore.method.get_num_buffers
}

void SIDRE_datastore_print(const SIDRE_datastore * self)
{
// splicer begin class.DataStore.method.print
  const DataStore * SH_this =
    static_cast<const DataStore *>(static_cast<const void *>(self));
  SH_this->print();
  return;
// splicer end class.DataStore.method.print
}

// splicer begin class.DataStore.additional_functions
// splicer end class.DataStore.additional_functions

}  // namespace axom
}  // namespace sidre
}  // extern "C"

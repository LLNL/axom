// wrapDataStore.cpp
// This is generated code, do not edit
//
// Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
//
// Produced at the Lawrence Livermore National Laboratory
//
// LLNL-CODE-741217
//
// All rights reserved.
//
// This file is part of Axom.
//
// For details about use and distribution, please read axom/LICENSE.
//
#include "wrapDataStore.h"
#include "sidre/DataStore.hpp"
#include "sidre/SidreTypes.hpp"

namespace axom
{
namespace sidre
{

// splicer begin class.DataStore.CXX_definitions
// splicer end class.DataStore.CXX_definitions

extern "C" {

// splicer begin class.DataStore.C_definitions
// splicer end class.DataStore.C_definitions

SIDRE_datastore* SIDRE_datastore_new()
{
// splicer begin class.DataStore.method.new
  DataStore* SHCXX_rv = new DataStore();
  SIDRE_datastore* SHC_rv =
    static_cast<SIDRE_datastore*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.DataStore.method.new
}

void SIDRE_datastore_delete(SIDRE_datastore* self)
{
// splicer begin class.DataStore.method.delete
  DataStore* SH_this = static_cast<DataStore*>(static_cast<void*>(self));
  delete SH_this;
  return;
// splicer end class.DataStore.method.delete
}

SIDRE_group* SIDRE_datastore_get_root(SIDRE_datastore* self)
{
// splicer begin class.DataStore.method.get_root
  DataStore* SH_this = static_cast<DataStore*>(static_cast<void*>(self));
  Group* SHCXX_rv = SH_this->getRoot();
  SIDRE_group* SHC_rv = static_cast<SIDRE_group*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.DataStore.method.get_root
}

SIDRE_buffer* SIDRE_datastore_get_buffer(SIDRE_datastore* self,
                                         SIDRE_IndexType idx)
{
// splicer begin class.DataStore.method.get_buffer
  DataStore* SH_this = static_cast<DataStore*>(static_cast<void*>(self));
  Buffer* SHCXX_rv = SH_this->getBuffer(idx);
  SIDRE_buffer* SHC_rv =
    static_cast<SIDRE_buffer*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.DataStore.method.get_buffer
}

SIDRE_buffer* SIDRE_datastore_create_buffer_empty(SIDRE_datastore* self)
{
// splicer begin class.DataStore.method.create_buffer_empty
  DataStore* SH_this = static_cast<DataStore*>(static_cast<void*>(self));
  Buffer* SHCXX_rv = SH_this->createBuffer();
  SIDRE_buffer* SHC_rv =
    static_cast<SIDRE_buffer*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.DataStore.method.create_buffer_empty
}

SIDRE_buffer* SIDRE_datastore_create_buffer_from_type(SIDRE_datastore* self,
                                                      int type,
                                                      SIDRE_SidreLength num_elems)
{
// splicer begin class.DataStore.method.create_buffer_from_type
  DataStore* SH_this = static_cast<DataStore*>(static_cast<void*>(self));
  TypeID SHCXX_type = getTypeID(type);
  Buffer* SHCXX_rv = SH_this->createBuffer(SHCXX_type, num_elems);
  SIDRE_buffer* SHC_rv =
    static_cast<SIDRE_buffer*>(static_cast<void*>(SHCXX_rv));
  return SHC_rv;
// splicer end class.DataStore.method.create_buffer_from_type
}

void SIDRE_datastore_destroy_buffer(SIDRE_datastore* self, SIDRE_IndexType id)
{
// splicer begin class.DataStore.method.destroy_buffer
  DataStore* SH_this = static_cast<DataStore*>(static_cast<void*>(self));
  SH_this->destroyBuffer(id);
  return;
// splicer end class.DataStore.method.destroy_buffer
}

size_t SIDRE_datastore_get_num_buffers(const SIDRE_datastore* self)
{
// splicer begin class.DataStore.method.get_num_buffers
  const DataStore* SH_this =
    static_cast<const DataStore*>(static_cast<const void*>(self));
  size_t SHC_rv = SH_this->getNumBuffers();
  return SHC_rv;
// splicer end class.DataStore.method.get_num_buffers
}

void SIDRE_datastore_print(const SIDRE_datastore* self)
{
// splicer begin class.DataStore.method.print
  const DataStore* SH_this =
    static_cast<const DataStore*>(static_cast<const void*>(self));
  SH_this->print();
  return;
// splicer end class.DataStore.method.print
}

}  // extern "C"

}  // namespace sidre
}  // namespace axom

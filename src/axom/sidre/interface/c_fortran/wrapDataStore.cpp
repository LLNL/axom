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
#include <stdlib.h>
#include "axom/sidre/core/Buffer.hpp"
#include "axom/sidre/core/DataStore.hpp"
#include "axom/sidre/core/Group.hpp"
#include "axom/sidre/core/SidreTypes.hpp"

// splicer begin class.DataStore.CXX_definitions
// splicer end class.DataStore.CXX_definitions

extern "C" {

// splicer begin class.DataStore.C_definitions
// splicer end class.DataStore.C_definitions

SIDRE_datastore* SIDRE_datastore_new(SIDRE_datastore* SHC_rv)
{
// splicer begin class.DataStore.method.new
  axom::sidre::DataStore* SHCXX_rv = new axom::sidre::DataStore();
  SHC_rv->addr = static_cast<void*>(SHCXX_rv);
  SHC_rv->idtor = 0;
  return SHC_rv;
// splicer end class.DataStore.method.new
}

void SIDRE_datastore_delete(SIDRE_datastore* self)
{
// splicer begin class.DataStore.method.delete
  axom::sidre::DataStore* SH_this =
    static_cast<axom::sidre::DataStore*>(self->addr);
  delete SH_this;
  self->addr = NULL;
  return;
// splicer end class.DataStore.method.delete
}

SIDRE_group* SIDRE_datastore_get_root(SIDRE_datastore* self,
                                      SIDRE_group* SHC_rv)
{
// splicer begin class.DataStore.method.get_root
  axom::sidre::DataStore* SH_this =
    static_cast<axom::sidre::DataStore*>(self->addr);
  axom::sidre::Group* SHCXX_rv = SH_this->getRoot();
  SHC_rv->addr = static_cast<void*>(SHCXX_rv);
  SHC_rv->idtor = 0;
  return SHC_rv;
// splicer end class.DataStore.method.get_root
}

size_t SIDRE_datastore_get_num_buffers(const SIDRE_datastore* self)
{
// splicer begin class.DataStore.method.get_num_buffers
  const axom::sidre::DataStore* SH_this =
    static_cast<const axom::sidre::DataStore*>(self->addr);
  size_t SHC_rv = SH_this->getNumBuffers();
  return SHC_rv;
// splicer end class.DataStore.method.get_num_buffers
}

SIDRE_buffer* SIDRE_datastore_get_buffer(SIDRE_datastore* self,
                                         SIDRE_IndexType idx,
                                         SIDRE_buffer* SHC_rv)
{
// splicer begin class.DataStore.method.get_buffer
  axom::sidre::DataStore* SH_this =
    static_cast<axom::sidre::DataStore*>(self->addr);
  axom::sidre::Buffer* SHCXX_rv = SH_this->getBuffer(idx);
  // C_error_pattern
  if (SHCXX_rv == nullptr)
  {
    SHC_rv->addr = NULL;
    SHC_rv->idtor = 0;
    return NULL;
  }

  SHC_rv->addr = static_cast<void*>(SHCXX_rv);
  SHC_rv->idtor = 0;
  return SHC_rv;
// splicer end class.DataStore.method.get_buffer
}

SIDRE_buffer* SIDRE_datastore_create_buffer_empty(SIDRE_datastore* self,
                                                  SIDRE_buffer* SHC_rv)
{
// splicer begin class.DataStore.method.create_buffer_empty
  axom::sidre::DataStore* SH_this =
    static_cast<axom::sidre::DataStore*>(self->addr);
  axom::sidre::Buffer* SHCXX_rv = SH_this->createBuffer();
  SHC_rv->addr = static_cast<void*>(SHCXX_rv);
  SHC_rv->idtor = 0;
  return SHC_rv;
// splicer end class.DataStore.method.create_buffer_empty
}

SIDRE_buffer* SIDRE_datastore_create_buffer_from_type(SIDRE_datastore* self,
                                                      int type,
                                                      SIDRE_SidreLength num_elems,
                                                      SIDRE_buffer* SHC_rv)
{
// splicer begin class.DataStore.method.create_buffer_from_type
  axom::sidre::DataStore* SH_this =
    static_cast<axom::sidre::DataStore*>(self->addr);
  axom::sidre::TypeID SHCXX_type = axom::sidre::getTypeID(type);
  axom::sidre::Buffer* SHCXX_rv = SH_this->createBuffer(SHCXX_type, num_elems);
  // C_error_pattern
  if (SHCXX_rv == nullptr)
  {
    SHC_rv->addr = NULL;
    SHC_rv->idtor = 0;
    return NULL;
  }

  SHC_rv->addr = static_cast<void*>(SHCXX_rv);
  SHC_rv->idtor = 0;
  return SHC_rv;
// splicer end class.DataStore.method.create_buffer_from_type
}

void SIDRE_datastore_destroy_buffer(SIDRE_datastore* self, SIDRE_IndexType id)
{
// splicer begin class.DataStore.method.destroy_buffer
  axom::sidre::DataStore* SH_this =
    static_cast<axom::sidre::DataStore*>(self->addr);
  SH_this->destroyBuffer(id);
  return;
// splicer end class.DataStore.method.destroy_buffer
}

void SIDRE_datastore_print(const SIDRE_datastore* self)
{
// splicer begin class.DataStore.method.print
  const axom::sidre::DataStore* SH_this =
    static_cast<const axom::sidre::DataStore*>(self->addr);
  SH_this->print();
  return;
// splicer end class.DataStore.method.print
}

}  // extern "C"

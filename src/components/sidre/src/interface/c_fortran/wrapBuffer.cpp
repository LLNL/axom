// wrapBuffer.cpp
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
#include "wrapBuffer.h"
#include "sidre/Buffer.hpp"
#include "sidre/SidreTypes.hpp"

// splicer begin class.Buffer.CXX_definitions
// splicer end class.Buffer.CXX_definitions

extern "C" {

// splicer begin class.Buffer.C_definitions
// splicer end class.Buffer.C_definitions

SIDRE_IndexType SIDRE_buffer_get_index(const SIDRE_buffer* self)
{
// splicer begin class.Buffer.method.get_index
  const axom::sidre::Buffer* SH_this =
    static_cast<const axom::sidre::Buffer*>(static_cast<const void*>(self));
  axom::sidre::IndexType SHC_rv = SH_this->getIndex();
  return SHC_rv;
// splicer end class.Buffer.method.get_index
}

size_t SIDRE_buffer_get_num_views(const SIDRE_buffer* self)
{
// splicer begin class.Buffer.method.get_num_views
  const axom::sidre::Buffer* SH_this =
    static_cast<const axom::sidre::Buffer*>(static_cast<const void*>(self));
  size_t SHC_rv = SH_this->getNumViews();
  return SHC_rv;
// splicer end class.Buffer.method.get_num_views
}

void* SIDRE_buffer_get_void_ptr(SIDRE_buffer* self)
{
// splicer begin class.Buffer.method.get_void_ptr
  axom::sidre::Buffer* SH_this =
    static_cast<axom::sidre::Buffer*>(static_cast<void*>(self));
  void* SHC_rv = SH_this->getVoidPtr();
  return SHC_rv;
// splicer end class.Buffer.method.get_void_ptr
}

int SIDRE_buffer_get_type_id(const SIDRE_buffer* self)
{
// splicer begin class.Buffer.method.get_type_id
  const axom::sidre::Buffer* SH_this =
    static_cast<const axom::sidre::Buffer*>(static_cast<const void*>(self));
  axom::sidre::TypeID SHCXX_rv = SH_this->getTypeID();
  int SHC_rv = static_cast<int>(SHCXX_rv);
  return SHC_rv;
// splicer end class.Buffer.method.get_type_id
}

size_t SIDRE_buffer_get_num_elements(const SIDRE_buffer* self)
{
// splicer begin class.Buffer.method.get_num_elements
  const axom::sidre::Buffer* SH_this =
    static_cast<const axom::sidre::Buffer*>(static_cast<const void*>(self));
  size_t SHC_rv = SH_this->getNumElements();
  return SHC_rv;
// splicer end class.Buffer.method.get_num_elements
}

size_t SIDRE_buffer_get_total_bytes(const SIDRE_buffer* self)
{
// splicer begin class.Buffer.method.get_total_bytes
  const axom::sidre::Buffer* SH_this =
    static_cast<const axom::sidre::Buffer*>(static_cast<const void*>(self));
  size_t SHC_rv = SH_this->getTotalBytes();
  return SHC_rv;
// splicer end class.Buffer.method.get_total_bytes
}

size_t SIDRE_buffer_get_bytes_per_element(const SIDRE_buffer* self)
{
// splicer begin class.Buffer.method.get_bytes_per_element
  const axom::sidre::Buffer* SH_this =
    static_cast<const axom::sidre::Buffer*>(static_cast<const void*>(self));
  size_t SHC_rv = SH_this->getBytesPerElement();
  return SHC_rv;
// splicer end class.Buffer.method.get_bytes_per_element
}

void SIDRE_buffer_describe(SIDRE_buffer* self, int type,
                           SIDRE_SidreLength num_elems)
{
// splicer begin class.Buffer.method.describe
  axom::sidre::Buffer* SH_this =
    static_cast<axom::sidre::Buffer*>(static_cast<void*>(self));
  axom::sidre::TypeID SHCXX_type = axom::sidre::getTypeID(type);
  SH_this->describe(SHCXX_type, num_elems);
  return;
// splicer end class.Buffer.method.describe
}

void SIDRE_buffer_allocate_existing(SIDRE_buffer* self)
{
// splicer begin class.Buffer.method.allocate_existing
  axom::sidre::Buffer* SH_this =
    static_cast<axom::sidre::Buffer*>(static_cast<void*>(self));
  SH_this->allocate();
  return;
// splicer end class.Buffer.method.allocate_existing
}

void SIDRE_buffer_allocate_from_type(SIDRE_buffer* self, int type,
                                     SIDRE_SidreLength num_elems)
{
// splicer begin class.Buffer.method.allocate_from_type
  axom::sidre::Buffer* SH_this =
    static_cast<axom::sidre::Buffer*>(static_cast<void*>(self));
  axom::sidre::TypeID SHCXX_type = axom::sidre::getTypeID(type);
  SH_this->allocate(SHCXX_type, num_elems);
  return;
// splicer end class.Buffer.method.allocate_from_type
}

void SIDRE_buffer_reallocate(SIDRE_buffer* self, SIDRE_SidreLength num_elems)
{
// splicer begin class.Buffer.method.reallocate
  axom::sidre::Buffer* SH_this =
    static_cast<axom::sidre::Buffer*>(static_cast<void*>(self));
  SH_this->reallocate(num_elems);
  return;
// splicer end class.Buffer.method.reallocate
}

void SIDRE_buffer_print(const SIDRE_buffer* self)
{
// splicer begin class.Buffer.method.print
  const axom::sidre::Buffer* SH_this =
    static_cast<const axom::sidre::Buffer*>(static_cast<const void*>(self));
  SH_this->print();
  return;
// splicer end class.Buffer.method.print
}

}  // extern "C"

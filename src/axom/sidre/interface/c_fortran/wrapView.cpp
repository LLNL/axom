// wrapView.cpp
// This is generated code, do not edit
//
// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "wrapView.h"
#include <cstring>
#include <stdlib.h>
#include <string>
#include "axom/sidre/core/Buffer.hpp"
#include "axom/sidre/core/Group.hpp"
#include "axom/sidre/core/SidreTypes.hpp"
#include "axom/sidre/core/View.hpp"

// splicer begin class.View.CXX_definitions
// splicer end class.View.CXX_definitions

extern "C" {


// helper function
// Copy src into dest, blank fill to ndest characters
// Truncate if dest is too short.
// dest will not be NULL terminated.
static void ShroudStrCopy(char* dest, int ndest, const char* src, int nsrc)
{
  int nm = nsrc < ndest ? nsrc : ndest;
  std::memcpy(dest,src,nm);
  if(ndest > nm)
    std::memset(dest+nm,' ',ndest-nm);
}
// splicer begin class.View.C_definitions
// splicer end class.View.C_definitions

SIDRE_IndexType SIDRE_view_get_index(SIDRE_view* self)
{
// splicer begin class.View.method.get_index
  axom::sidre::View* SH_this = static_cast<axom::sidre::View*>(self->addr);
  axom::sidre::IndexType SHC_rv = SH_this->getIndex();
  return SHC_rv;
// splicer end class.View.method.get_index
}

const char* SIDRE_view_get_name(const SIDRE_view* self)
{
// splicer begin class.View.method.get_name
  const axom::sidre::View* SH_this =
    static_cast<const axom::sidre::View*>(self->addr);
  const std::string & SHCXX_rv = SH_this->getName();
  const char* SHC_rv = SHCXX_rv.c_str();
  return SHC_rv;
// splicer end class.View.method.get_name
}

void SIDRE_view_get_name_bufferify(const SIDRE_view* self, char* SHF_rv,
                                   int NSHF_rv)
{
// splicer begin class.View.method.get_name_bufferify
  const axom::sidre::View* SH_this =
    static_cast<const axom::sidre::View*>(self->addr);
  const std::string & SHCXX_rv = SH_this->getName();
  if (SHCXX_rv.empty())
  {
    std::memset(SHF_rv, ' ', NSHF_rv);
  }
  else
  {
    ShroudStrCopy(SHF_rv, NSHF_rv, SHCXX_rv.data(), SHCXX_rv.size());
  }
  return;
// splicer end class.View.method.get_name_bufferify
}

void SIDRE_view_get_path_bufferify(const SIDRE_view* self, char* SHF_rv,
                                   int NSHF_rv)
{
// splicer begin class.View.method.get_path_bufferify
  const axom::sidre::View* SH_this =
    static_cast<const axom::sidre::View*>(self->addr);
  std::string SHCXX_rv = SH_this->getPath();
  if (SHCXX_rv.empty())
  {
    std::memset(SHF_rv, ' ', NSHF_rv);
  }
  else
  {
    ShroudStrCopy(SHF_rv, NSHF_rv, SHCXX_rv.data(), SHCXX_rv.size());
  }
  return;
// splicer end class.View.method.get_path_bufferify
}

void SIDRE_view_get_path_name_bufferify(const SIDRE_view* self, char* SHF_rv,
                                        int NSHF_rv)
{
// splicer begin class.View.method.get_path_name_bufferify
  const axom::sidre::View* SH_this =
    static_cast<const axom::sidre::View*>(self->addr);
  std::string SHCXX_rv = SH_this->getPathName();
  if (SHCXX_rv.empty())
  {
    std::memset(SHF_rv, ' ', NSHF_rv);
  }
  else
  {
    ShroudStrCopy(SHF_rv, NSHF_rv, SHCXX_rv.data(), SHCXX_rv.size());
  }
  return;
// splicer end class.View.method.get_path_name_bufferify
}

SIDRE_group* SIDRE_view_get_owning_group(SIDRE_view* self, SIDRE_group* SHC_rv)
{
// splicer begin class.View.method.get_owning_group
  axom::sidre::View* SH_this = static_cast<axom::sidre::View*>(self->addr);
  axom::sidre::Group* SHCXX_rv = SH_this->getOwningGroup();
  SHC_rv->addr = static_cast<void*>(SHCXX_rv);
  SHC_rv->idtor = 0;
  return SHC_rv;
// splicer end class.View.method.get_owning_group
}

bool SIDRE_view_has_buffer(const SIDRE_view* self)
{
// splicer begin class.View.method.has_buffer
  const axom::sidre::View* SH_this =
    static_cast<const axom::sidre::View*>(self->addr);
  bool SHC_rv = SH_this->hasBuffer();
  return SHC_rv;
// splicer end class.View.method.has_buffer
}

SIDRE_buffer* SIDRE_view_get_buffer(SIDRE_view* self, SIDRE_buffer* SHC_rv)
{
// splicer begin class.View.method.get_buffer
  axom::sidre::View* SH_this = static_cast<axom::sidre::View*>(self->addr);
  axom::sidre::Buffer* SHCXX_rv = SH_this->getBuffer();
  SHC_rv->addr = static_cast<void*>(SHCXX_rv);
  SHC_rv->idtor = 0;
  return SHC_rv;
// splicer end class.View.method.get_buffer
}

bool SIDRE_view_is_external(const SIDRE_view* self)
{
// splicer begin class.View.method.is_external
  const axom::sidre::View* SH_this =
    static_cast<const axom::sidre::View*>(self->addr);
  bool SHC_rv = SH_this->isExternal();
  return SHC_rv;
// splicer end class.View.method.is_external
}

bool SIDRE_view_is_allocated(SIDRE_view* self)
{
// splicer begin class.View.method.is_allocated
  axom::sidre::View* SH_this = static_cast<axom::sidre::View*>(self->addr);
  bool SHC_rv = SH_this->isAllocated();
  return SHC_rv;
// splicer end class.View.method.is_allocated
}

bool SIDRE_view_is_applied(const SIDRE_view* self)
{
// splicer begin class.View.method.is_applied
  const axom::sidre::View* SH_this =
    static_cast<const axom::sidre::View*>(self->addr);
  bool SHC_rv = SH_this->isApplied();
  return SHC_rv;
// splicer end class.View.method.is_applied
}

bool SIDRE_view_is_described(const SIDRE_view* self)
{
// splicer begin class.View.method.is_described
  const axom::sidre::View* SH_this =
    static_cast<const axom::sidre::View*>(self->addr);
  bool SHC_rv = SH_this->isDescribed();
  return SHC_rv;
// splicer end class.View.method.is_described
}

bool SIDRE_view_is_empty(const SIDRE_view* self)
{
// splicer begin class.View.method.is_empty
  const axom::sidre::View* SH_this =
    static_cast<const axom::sidre::View*>(self->addr);
  bool SHC_rv = SH_this->isEmpty();
  return SHC_rv;
// splicer end class.View.method.is_empty
}

bool SIDRE_view_is_opaque(const SIDRE_view* self)
{
// splicer begin class.View.method.is_opaque
  const axom::sidre::View* SH_this =
    static_cast<const axom::sidre::View*>(self->addr);
  bool SHC_rv = SH_this->isOpaque();
  return SHC_rv;
// splicer end class.View.method.is_opaque
}

bool SIDRE_view_is_scalar(const SIDRE_view* self)
{
// splicer begin class.View.method.is_scalar
  const axom::sidre::View* SH_this =
    static_cast<const axom::sidre::View*>(self->addr);
  bool SHC_rv = SH_this->isScalar();
  return SHC_rv;
// splicer end class.View.method.is_scalar
}

bool SIDRE_view_is_string(const SIDRE_view* self)
{
// splicer begin class.View.method.is_string
  const axom::sidre::View* SH_this =
    static_cast<const axom::sidre::View*>(self->addr);
  bool SHC_rv = SH_this->isString();
  return SHC_rv;
// splicer end class.View.method.is_string
}

int SIDRE_view_get_type_id(const SIDRE_view* self)
{
// splicer begin class.View.method.get_type_id
  const axom::sidre::View* SH_this =
    static_cast<const axom::sidre::View*>(self->addr);
  axom::sidre::TypeID SHCXX_rv = SH_this->getTypeID();
  int SHC_rv = static_cast<int>(SHCXX_rv);
  return SHC_rv;
// splicer end class.View.method.get_type_id
}

size_t SIDRE_view_get_total_bytes(const SIDRE_view* self)
{
// splicer begin class.View.method.get_total_bytes
  const axom::sidre::View* SH_this =
    static_cast<const axom::sidre::View*>(self->addr);
  size_t SHC_rv = SH_this->getTotalBytes();
  return SHC_rv;
// splicer end class.View.method.get_total_bytes
}

size_t SIDRE_view_get_num_elements(const SIDRE_view* self)
{
// splicer begin class.View.method.get_num_elements
  const axom::sidre::View* SH_this =
    static_cast<const axom::sidre::View*>(self->addr);
  size_t SHC_rv = SH_this->getNumElements();
  return SHC_rv;
// splicer end class.View.method.get_num_elements
}

size_t SIDRE_view_get_bytes_per_element(const SIDRE_view* self)
{
// splicer begin class.View.method.get_bytes_per_element
  const axom::sidre::View* SH_this =
    static_cast<const axom::sidre::View*>(self->addr);
  size_t SHC_rv = SH_this->getBytesPerElement();
  return SHC_rv;
// splicer end class.View.method.get_bytes_per_element
}

size_t SIDRE_view_get_offset(const SIDRE_view* self)
{
// splicer begin class.View.method.get_offset
  const axom::sidre::View* SH_this =
    static_cast<const axom::sidre::View*>(self->addr);
  size_t SHC_rv = SH_this->getOffset();
  return SHC_rv;
// splicer end class.View.method.get_offset
}

size_t SIDRE_view_get_stride(const SIDRE_view* self)
{
// splicer begin class.View.method.get_stride
  const axom::sidre::View* SH_this =
    static_cast<const axom::sidre::View*>(self->addr);
  size_t SHC_rv = SH_this->getStride();
  return SHC_rv;
// splicer end class.View.method.get_stride
}

int SIDRE_view_get_num_dimensions(const SIDRE_view* self)
{
// splicer begin class.View.method.get_num_dimensions
  const axom::sidre::View* SH_this =
    static_cast<const axom::sidre::View*>(self->addr);
  int SHC_rv = SH_this->getNumDimensions();
  return SHC_rv;
// splicer end class.View.method.get_num_dimensions
}

int SIDRE_view_get_shape(const SIDRE_view* self, int ndims,
                         SIDRE_IndexType* shape)
{
// splicer begin class.View.method.get_shape
  const axom::sidre::View* SH_this =
    static_cast<const axom::sidre::View*>(self->addr);
  int SHC_rv = SH_this->getShape(ndims, shape);
  return SHC_rv;
// splicer end class.View.method.get_shape
}

void SIDRE_view_allocate_simple(SIDRE_view* self)
{
// splicer begin class.View.method.allocate_simple
  axom::sidre::View* SH_this = static_cast<axom::sidre::View*>(self->addr);
  SH_this->allocate();
  return;
// splicer end class.View.method.allocate_simple
}

void SIDRE_view_allocate_from_type(SIDRE_view* self, int type,
                                   SIDRE_IndexType num_elems)
{
// splicer begin class.View.method.allocate_from_type
  axom::sidre::View* SH_this = static_cast<axom::sidre::View*>(self->addr);
  axom::sidre::TypeID SHCXX_type = axom::sidre::getTypeID(type);
  SH_this->allocate(SHCXX_type, num_elems);
  return;
// splicer end class.View.method.allocate_from_type
}

void SIDRE_view_reallocate(SIDRE_view* self, SIDRE_IndexType num_elems)
{
// splicer begin class.View.method.reallocate
  axom::sidre::View* SH_this = static_cast<axom::sidre::View*>(self->addr);
  SH_this->reallocate(num_elems);
  return;
// splicer end class.View.method.reallocate
}

void SIDRE_view_attach_buffer_only(SIDRE_view* self, SIDRE_buffer* buff)
{
// splicer begin class.View.method.attach_buffer_only
  axom::sidre::View* SH_this = static_cast<axom::sidre::View*>(self->addr);
  axom::sidre::Buffer* SHCXX_buff =
    static_cast<axom::sidre::Buffer*>(buff->addr);
  SH_this->attachBuffer(SHCXX_buff);
  return;
// splicer end class.View.method.attach_buffer_only
}

void SIDRE_view_attach_buffer_type(SIDRE_view* self, int type,
                                   SIDRE_IndexType num_elems,
                                   SIDRE_buffer* buff)
{
// splicer begin class.View.method.attach_buffer_type
  axom::sidre::View* SH_this = static_cast<axom::sidre::View*>(self->addr);
  axom::sidre::TypeID SHCXX_type = axom::sidre::getTypeID(type);
  axom::sidre::Buffer* SHCXX_buff =
    static_cast<axom::sidre::Buffer*>(buff->addr);
  SH_this->attachBuffer(SHCXX_type, num_elems, SHCXX_buff);
  return;
// splicer end class.View.method.attach_buffer_type
}

void SIDRE_view_attach_buffer_shape(SIDRE_view* self, int type, int ndims,
                                    SIDRE_IndexType* shape, SIDRE_buffer* buff)
{
// splicer begin class.View.method.attach_buffer_shape
  axom::sidre::View* SH_this = static_cast<axom::sidre::View*>(self->addr);
  axom::sidre::TypeID SHCXX_type = axom::sidre::getTypeID(type);
  axom::sidre::Buffer* SHCXX_buff =
    static_cast<axom::sidre::Buffer*>(buff->addr);
  SH_this->attachBuffer(SHCXX_type, ndims, shape, SHCXX_buff);
  return;
// splicer end class.View.method.attach_buffer_shape
}

void SIDRE_view_apply_0(SIDRE_view* self)
{
// splicer begin class.View.method.apply_0
  axom::sidre::View* SH_this = static_cast<axom::sidre::View*>(self->addr);
  SH_this->apply();
  return;
// splicer end class.View.method.apply_0
}

void SIDRE_view_apply_nelems(SIDRE_view* self, SIDRE_IndexType num_elems)
{
// splicer begin class.View.method.apply_nelems
  axom::sidre::View* SH_this = static_cast<axom::sidre::View*>(self->addr);
  SH_this->apply(num_elems);
  return;
// splicer end class.View.method.apply_nelems
}

void SIDRE_view_apply_nelems_offset(SIDRE_view* self, SIDRE_IndexType num_elems,
                                    SIDRE_IndexType offset)
{
// splicer begin class.View.method.apply_nelems_offset
  axom::sidre::View* SH_this = static_cast<axom::sidre::View*>(self->addr);
  SH_this->apply(num_elems, offset);
  return;
// splicer end class.View.method.apply_nelems_offset
}

void SIDRE_view_apply_nelems_offset_stride(SIDRE_view* self,
                                           SIDRE_IndexType num_elems,
                                           SIDRE_IndexType offset,
                                           SIDRE_IndexType stride)
{
// splicer begin class.View.method.apply_nelems_offset_stride
  axom::sidre::View* SH_this = static_cast<axom::sidre::View*>(self->addr);
  SH_this->apply(num_elems, offset, stride);
  return;
// splicer end class.View.method.apply_nelems_offset_stride
}

void SIDRE_view_apply_type_nelems(SIDRE_view* self, int type,
                                  SIDRE_IndexType num_elems)
{
// splicer begin class.View.method.apply_type_nelems
  axom::sidre::View* SH_this = static_cast<axom::sidre::View*>(self->addr);
  axom::sidre::TypeID SHCXX_type = axom::sidre::getTypeID(type);
  SH_this->apply(SHCXX_type, num_elems);
  return;
// splicer end class.View.method.apply_type_nelems
}

void SIDRE_view_apply_type_nelems_offset(SIDRE_view* self, int type,
                                         SIDRE_IndexType num_elems,
                                         SIDRE_IndexType offset)
{
// splicer begin class.View.method.apply_type_nelems_offset
  axom::sidre::View* SH_this = static_cast<axom::sidre::View*>(self->addr);
  axom::sidre::TypeID SHCXX_type = axom::sidre::getTypeID(type);
  SH_this->apply(SHCXX_type, num_elems, offset);
  return;
// splicer end class.View.method.apply_type_nelems_offset
}

void SIDRE_view_apply_type_nelems_offset_stride(SIDRE_view* self, int type,
                                                SIDRE_IndexType num_elems,
                                                SIDRE_IndexType offset,
                                                SIDRE_IndexType stride)
{
// splicer begin class.View.method.apply_type_nelems_offset_stride
  axom::sidre::View* SH_this = static_cast<axom::sidre::View*>(self->addr);
  axom::sidre::TypeID SHCXX_type = axom::sidre::getTypeID(type);
  SH_this->apply(SHCXX_type, num_elems, offset, stride);
  return;
// splicer end class.View.method.apply_type_nelems_offset_stride
}

void SIDRE_view_apply_type_shape(SIDRE_view* self, int type, int ndims,
                                 SIDRE_IndexType* shape)
{
// splicer begin class.View.method.apply_type_shape
  axom::sidre::View* SH_this = static_cast<axom::sidre::View*>(self->addr);
  axom::sidre::TypeID SHCXX_type = axom::sidre::getTypeID(type);
  SH_this->apply(SHCXX_type, ndims, shape);
  return;
// splicer end class.View.method.apply_type_shape
}

void SIDRE_view_set_scalar_int(SIDRE_view* self, int value)
{
// splicer begin class.View.method.set_scalar_int
  axom::sidre::View* SH_this = static_cast<axom::sidre::View*>(self->addr);
  SH_this->setScalar<int>(value);
  return;
// splicer end class.View.method.set_scalar_int
}

void SIDRE_view_set_scalar_long(SIDRE_view* self, long value)
{
// splicer begin class.View.method.set_scalar_long
  axom::sidre::View* SH_this = static_cast<axom::sidre::View*>(self->addr);
  SH_this->setScalar<long>(value);
  return;
// splicer end class.View.method.set_scalar_long
}

void SIDRE_view_set_scalar_float(SIDRE_view* self, float value)
{
// splicer begin class.View.method.set_scalar_float
  axom::sidre::View* SH_this = static_cast<axom::sidre::View*>(self->addr);
  SH_this->setScalar<float>(value);
  return;
// splicer end class.View.method.set_scalar_float
}

void SIDRE_view_set_scalar_double(SIDRE_view* self, double value)
{
// splicer begin class.View.method.set_scalar_double
  axom::sidre::View* SH_this = static_cast<axom::sidre::View*>(self->addr);
  SH_this->setScalar<double>(value);
  return;
// splicer end class.View.method.set_scalar_double
}

void SIDRE_view_set_string(SIDRE_view* self, const char* value)
{
// splicer begin class.View.method.set_string
  axom::sidre::View* SH_this = static_cast<axom::sidre::View*>(self->addr);
  const std::string SH_value(value);
  SH_this->setString(SH_value);
  return;
// splicer end class.View.method.set_string
}

void SIDRE_view_set_string_bufferify(SIDRE_view* self, const char* value,
                                     int Lvalue)
{
// splicer begin class.View.method.set_string_bufferify
  axom::sidre::View* SH_this = static_cast<axom::sidre::View*>(self->addr);
  const std::string SH_value(value, Lvalue);
  SH_this->setString(SH_value);
  return;
// splicer end class.View.method.set_string_bufferify
}

void SIDRE_view_set_external_data_ptr_only(SIDRE_view* self, void* external_ptr)
{
// splicer begin class.View.method.set_external_data_ptr_only
  axom::sidre::View* SH_this = static_cast<axom::sidre::View*>(self->addr);
  SH_this->setExternalDataPtr(external_ptr);
  return;
// splicer end class.View.method.set_external_data_ptr_only
}

void SIDRE_view_set_external_data_ptr_type(SIDRE_view* self, int type,
                                           SIDRE_IndexType num_elems,
                                           void* external_ptr)
{
// splicer begin class.View.method.set_external_data_ptr_type
  axom::sidre::View* SH_this = static_cast<axom::sidre::View*>(self->addr);
  axom::sidre::TypeID SHCXX_type = axom::sidre::getTypeID(type);
  SH_this->setExternalDataPtr(SHCXX_type, num_elems, external_ptr);
  return;
// splicer end class.View.method.set_external_data_ptr_type
}

void SIDRE_view_set_external_data_ptr_shape(SIDRE_view* self, int type,
                                            int ndims, SIDRE_IndexType* shape,
                                            void* external_ptr)
{
// splicer begin class.View.method.set_external_data_ptr_shape
  axom::sidre::View* SH_this = static_cast<axom::sidre::View*>(self->addr);
  axom::sidre::TypeID SHCXX_type = axom::sidre::getTypeID(type);
  SH_this->setExternalDataPtr(SHCXX_type, ndims, shape, external_ptr);
  return;
// splicer end class.View.method.set_external_data_ptr_shape
}

const char* SIDRE_view_get_string(SIDRE_view* self)
{
// splicer begin class.View.method.get_string
  axom::sidre::View* SH_this = static_cast<axom::sidre::View*>(self->addr);
  const char* SHC_rv = SH_this->getString();
  return SHC_rv;
// splicer end class.View.method.get_string
}

void SIDRE_view_get_string_bufferify(SIDRE_view* self, char* name, int Nname)
{
// splicer begin class.View.method.get_string_bufferify
  axom::sidre::View* SH_this = static_cast<axom::sidre::View*>(self->addr);
  const char* SHC_rv = SH_this->getString();
  if (SHC_rv == NULL)
  {
    std::memset(name, ' ', Nname);
  }
  else
  {
    ShroudStrCopy(name, Nname, SHC_rv, std::strlen(SHC_rv));
  }
  return;
// splicer end class.View.method.get_string_bufferify
}

int SIDRE_view_get_data_int(SIDRE_view* self)
{
// splicer begin class.View.method.get_data_int
  axom::sidre::View* SH_this = static_cast<axom::sidre::View*>(self->addr);
  int SHC_rv = SH_this->getData<int>();
  return SHC_rv;
// splicer end class.View.method.get_data_int
}

long SIDRE_view_get_data_long(SIDRE_view* self)
{
// splicer begin class.View.method.get_data_long
  axom::sidre::View* SH_this = static_cast<axom::sidre::View*>(self->addr);
  long SHC_rv = SH_this->getData<long>();
  return SHC_rv;
// splicer end class.View.method.get_data_long
}

float SIDRE_view_get_data_float(SIDRE_view* self)
{
// splicer begin class.View.method.get_data_float
  axom::sidre::View* SH_this = static_cast<axom::sidre::View*>(self->addr);
  float SHC_rv = SH_this->getData<float>();
  return SHC_rv;
// splicer end class.View.method.get_data_float
}

double SIDRE_view_get_data_double(SIDRE_view* self)
{
// splicer begin class.View.method.get_data_double
  axom::sidre::View* SH_this = static_cast<axom::sidre::View*>(self->addr);
  double SHC_rv = SH_this->getData<double>();
  return SHC_rv;
// splicer end class.View.method.get_data_double
}

void* SIDRE_view_get_void_ptr(const SIDRE_view* self)
{
// splicer begin class.View.method.get_void_ptr
  const axom::sidre::View* SH_this =
    static_cast<const axom::sidre::View*>(self->addr);
  void* SHC_rv = SH_this->getVoidPtr();
  return SHC_rv;
// splicer end class.View.method.get_void_ptr
}

void SIDRE_view_print(const SIDRE_view* self)
{
// splicer begin class.View.method.print
  const axom::sidre::View* SH_this =
    static_cast<const axom::sidre::View*>(self->addr);
  SH_this->print();
  return;
// splicer end class.View.method.print
}

bool SIDRE_view_rename(SIDRE_view* self, const char* new_name)
{
// splicer begin class.View.method.rename
  axom::sidre::View* SH_this = static_cast<axom::sidre::View*>(self->addr);
  const std::string SH_new_name(new_name);
  bool SHC_rv = SH_this->rename(SH_new_name);
  return SHC_rv;
// splicer end class.View.method.rename
}

bool SIDRE_view_rename_bufferify(SIDRE_view* self, const char* new_name,
                                 int Lnew_name)
{
// splicer begin class.View.method.rename_bufferify
  axom::sidre::View* SH_this = static_cast<axom::sidre::View*>(self->addr);
  const std::string SH_new_name(new_name, Lnew_name);
  bool SHC_rv = SH_this->rename(SH_new_name);
  return SHC_rv;
// splicer end class.View.method.rename_bufferify
}

}  // extern "C"

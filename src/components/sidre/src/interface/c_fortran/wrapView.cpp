// wrapView.cpp
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
// wrapView.cpp
#include "wrapView.h"
#include <cstring>
#include <string>
#include "shroudrt.hpp"
#include "sidre/SidreTypes.hpp"
#include "sidre/View.hpp"

namespace axom
{
namespace sidre
{

// splicer begin class.View.CXX_definitions
// splicer end class.View.CXX_definitions

extern "C" {

// splicer begin class.View.C_definitions
// splicer end class.View.C_definitions

void SIDRE_view_allocate_simple(SIDRE_view* self)
{
// splicer begin class.View.method.allocate_simple
  View* SH_this = static_cast<View*>(static_cast<void*>(self));
  SH_this->allocate();
  return;
// splicer end class.View.method.allocate_simple
}

void SIDRE_view_allocate_from_type(SIDRE_view* self, int type,
                                   SIDRE_SidreLength num_elems)
{
// splicer begin class.View.method.allocate_from_type
  View* SH_this = static_cast<View*>(static_cast<void*>(self));
  SH_this->allocate(getTypeID(type), num_elems);
  return;
// splicer end class.View.method.allocate_from_type
}

void SIDRE_view_reallocate(SIDRE_view* self, SIDRE_SidreLength num_elems)
{
// splicer begin class.View.method.reallocate
  View* SH_this = static_cast<View*>(static_cast<void*>(self));
  SH_this->reallocate(num_elems);
  return;
// splicer end class.View.method.reallocate
}

void SIDRE_view_apply_0(SIDRE_view* self)
{
// splicer begin class.View.method.apply_0
  View* SH_this = static_cast<View*>(static_cast<void*>(self));
  SH_this->apply();
  return;
// splicer end class.View.method.apply_0
}

void SIDRE_view_attach_buffer_only(SIDRE_view* self, SIDRE_buffer* buff)
{
// splicer begin class.View.method.attach_buffer_only
  View* SH_this = static_cast<View*>(static_cast<void*>(self));
  SH_this->attachBuffer(static_cast<Buffer*>(static_cast<void*>(buff)));
  return;
// splicer end class.View.method.attach_buffer_only
}

void SIDRE_view_attach_buffer_type(SIDRE_view* self, int type,
                                   SIDRE_SidreLength num_elems,
                                   SIDRE_buffer* buff)
{
// splicer begin class.View.method.attach_buffer_type
  View* SH_this = static_cast<View*>(static_cast<void*>(self));
  SH_this->attachBuffer(getTypeID(type), num_elems,
                        static_cast<Buffer*>(static_cast<void*>(buff)));
  return;
// splicer end class.View.method.attach_buffer_type
}

void SIDRE_view_attach_buffer_shape(SIDRE_view* self, int type, int ndims,
                                    SIDRE_SidreLength* shape,
                                    SIDRE_buffer* buff)
{
// splicer begin class.View.method.attach_buffer_shape
  View* SH_this = static_cast<View*>(static_cast<void*>(self));
  SH_this->attachBuffer(getTypeID(type), ndims, shape,
                        static_cast<Buffer*>(static_cast<void*>(buff)));
  return;
// splicer end class.View.method.attach_buffer_shape
}

void SIDRE_view_apply_nelems(SIDRE_view* self, SIDRE_SidreLength num_elems)
{
// splicer begin class.View.method.apply_nelems
  View* SH_this = static_cast<View*>(static_cast<void*>(self));
  SH_this->apply(num_elems);
  return;
// splicer end class.View.method.apply_nelems
}

void SIDRE_view_apply_nelems_offset(SIDRE_view* self,
                                    SIDRE_SidreLength num_elems,
                                    SIDRE_SidreLength offset)
{
// splicer begin class.View.method.apply_nelems_offset
  View* SH_this = static_cast<View*>(static_cast<void*>(self));
  SH_this->apply(num_elems, offset);
  return;
// splicer end class.View.method.apply_nelems_offset
}

void SIDRE_view_apply_nelems_offset_stride(SIDRE_view* self,
                                           SIDRE_SidreLength num_elems,
                                           SIDRE_SidreLength offset,
                                           SIDRE_SidreLength stride)
{
// splicer begin class.View.method.apply_nelems_offset_stride
  View* SH_this = static_cast<View*>(static_cast<void*>(self));
  SH_this->apply(num_elems, offset, stride);
  return;
// splicer end class.View.method.apply_nelems_offset_stride
}

void SIDRE_view_apply_type_nelems(SIDRE_view* self, int type,
                                  SIDRE_SidreLength num_elems)
{
// splicer begin class.View.method.apply_type_nelems
  View* SH_this = static_cast<View*>(static_cast<void*>(self));
  SH_this->apply(getTypeID(type), num_elems);
  return;
// splicer end class.View.method.apply_type_nelems
}

void SIDRE_view_apply_type_nelems_offset(SIDRE_view* self, int type,
                                         SIDRE_SidreLength num_elems,
                                         SIDRE_SidreLength offset)
{
// splicer begin class.View.method.apply_type_nelems_offset
  View* SH_this = static_cast<View*>(static_cast<void*>(self));
  SH_this->apply(getTypeID(type), num_elems, offset);
  return;
// splicer end class.View.method.apply_type_nelems_offset
}

void SIDRE_view_apply_type_nelems_offset_stride(SIDRE_view* self, int type,
                                                SIDRE_SidreLength num_elems,
                                                SIDRE_SidreLength offset,
                                                SIDRE_SidreLength stride)
{
// splicer begin class.View.method.apply_type_nelems_offset_stride
  View* SH_this = static_cast<View*>(static_cast<void*>(self));
  SH_this->apply(getTypeID(type), num_elems, offset, stride);
  return;
// splicer end class.View.method.apply_type_nelems_offset_stride
}

void SIDRE_view_apply_type_shape(SIDRE_view* self, int type, int ndims,
                                 SIDRE_SidreLength* shape)
{
// splicer begin class.View.method.apply_type_shape
  View* SH_this = static_cast<View*>(static_cast<void*>(self));
  SH_this->apply(getTypeID(type), ndims, shape);
  return;
// splicer end class.View.method.apply_type_shape
}

bool SIDRE_view_has_buffer(const SIDRE_view* self)
{
// splicer begin class.View.method.has_buffer
  const View* SH_this =
    static_cast<const View*>(static_cast<const void*>(self));
  bool SH_rv = SH_this->hasBuffer();
  return SH_rv;
// splicer end class.View.method.has_buffer
}

bool SIDRE_view_is_external(const SIDRE_view* self)
{
// splicer begin class.View.method.is_external
  const View* SH_this =
    static_cast<const View*>(static_cast<const void*>(self));
  bool SH_rv = SH_this->isExternal();
  return SH_rv;
// splicer end class.View.method.is_external
}

bool SIDRE_view_is_allocated(SIDRE_view* self)
{
// splicer begin class.View.method.is_allocated
  View* SH_this = static_cast<View*>(static_cast<void*>(self));
  bool SH_rv = SH_this->isAllocated();
  return SH_rv;
// splicer end class.View.method.is_allocated
}

bool SIDRE_view_is_applied(const SIDRE_view* self)
{
// splicer begin class.View.method.is_applied
  const View* SH_this =
    static_cast<const View*>(static_cast<const void*>(self));
  bool SH_rv = SH_this->isApplied();
  return SH_rv;
// splicer end class.View.method.is_applied
}

bool SIDRE_view_is_described(const SIDRE_view* self)
{
// splicer begin class.View.method.is_described
  const View* SH_this =
    static_cast<const View*>(static_cast<const void*>(self));
  bool SH_rv = SH_this->isDescribed();
  return SH_rv;
// splicer end class.View.method.is_described
}

bool SIDRE_view_is_empty(const SIDRE_view* self)
{
// splicer begin class.View.method.is_empty
  const View* SH_this =
    static_cast<const View*>(static_cast<const void*>(self));
  bool SH_rv = SH_this->isEmpty();
  return SH_rv;
// splicer end class.View.method.is_empty
}

bool SIDRE_view_is_opaque(const SIDRE_view* self)
{
// splicer begin class.View.method.is_opaque
  const View* SH_this =
    static_cast<const View*>(static_cast<const void*>(self));
  bool SH_rv = SH_this->isOpaque();
  return SH_rv;
// splicer end class.View.method.is_opaque
}

bool SIDRE_view_is_scalar(const SIDRE_view* self)
{
// splicer begin class.View.method.is_scalar
  const View* SH_this =
    static_cast<const View*>(static_cast<const void*>(self));
  bool SH_rv = SH_this->isScalar();
  return SH_rv;
// splicer end class.View.method.is_scalar
}

bool SIDRE_view_is_string(const SIDRE_view* self)
{
// splicer begin class.View.method.is_string
  const View* SH_this =
    static_cast<const View*>(static_cast<const void*>(self));
  bool SH_rv = SH_this->isString();
  return SH_rv;
// splicer end class.View.method.is_string
}

SIDRE_IndexType SIDRE_view_get_index(SIDRE_view* self)
{
// splicer begin class.View.method.get_index
  View* SH_this = static_cast<View*>(static_cast<void*>(self));
  IndexType SH_rv = SH_this->getIndex();
  return SH_rv;
// splicer end class.View.method.get_index
}

const char* SIDRE_view_get_name(const SIDRE_view* self)
{
// splicer begin class.View.method.get_name
  const View* SH_this =
    static_cast<const View*>(static_cast<const void*>(self));
  const std::string & SH_rv = SH_this->getName();
  const char* XSH_rv = SH_rv.c_str();
  return XSH_rv;
// splicer end class.View.method.get_name
}

void SIDRE_view_get_name_bufferify(const SIDRE_view* self, char* SH_F_rv,
                                   int NSH_F_rv)
{
// splicer begin class.View.method.get_name_bufferify
  const View* SH_this =
    static_cast<const View*>(static_cast<const void*>(self));
  const std::string & SH_rv = SH_this->getName();
  if (SH_rv.empty())
  {
    std::memset(SH_F_rv, ' ', NSH_F_rv);
  }
  else
  {
    shroud_FccCopy(SH_F_rv, NSH_F_rv, SH_rv.c_str());
  }
  return;
// splicer end class.View.method.get_name_bufferify
}

void SIDRE_view_get_path_bufferify(const SIDRE_view* self, char* SH_F_rv,
                                   int NSH_F_rv)
{
// splicer begin class.View.method.get_path_bufferify
  const View* SH_this =
    static_cast<const View*>(static_cast<const void*>(self));
  std::string SH_rv = SH_this->getPath();
  if (SH_rv.empty())
  {
    std::memset(SH_F_rv, ' ', NSH_F_rv);
  }
  else
  {
    shroud_FccCopy(SH_F_rv, NSH_F_rv, SH_rv.c_str());
  }
  return;
// splicer end class.View.method.get_path_bufferify
}

void SIDRE_view_get_path_name_bufferify(const SIDRE_view* self, char* SH_F_rv,
                                        int NSH_F_rv)
{
// splicer begin class.View.method.get_path_name_bufferify
  const View* SH_this =
    static_cast<const View*>(static_cast<const void*>(self));
  std::string SH_rv = SH_this->getPathName();
  if (SH_rv.empty())
  {
    std::memset(SH_F_rv, ' ', NSH_F_rv);
  }
  else
  {
    shroud_FccCopy(SH_F_rv, NSH_F_rv, SH_rv.c_str());
  }
  return;
// splicer end class.View.method.get_path_name_bufferify
}

SIDRE_buffer* SIDRE_view_get_buffer(SIDRE_view* self)
{
// splicer begin class.View.method.get_buffer
  View* SH_this = static_cast<View*>(static_cast<void*>(self));
  Buffer* SH_rv = SH_this->getBuffer();
  SIDRE_buffer* XSH_rv =
    static_cast<SIDRE_buffer*>(static_cast<void*>(SH_rv));
  return XSH_rv;
// splicer end class.View.method.get_buffer
}

void* SIDRE_view_get_void_ptr(const SIDRE_view* self)
{
// splicer begin class.View.method.get_void_ptr
  const View* SH_this =
    static_cast<const View*>(static_cast<const void*>(self));
  void* SH_rv = SH_this->getVoidPtr();
  return SH_rv;
// splicer end class.View.method.get_void_ptr
}

void SIDRE_view_set_scalar_int(SIDRE_view* self, int value)
{
// splicer begin class.View.method.set_scalar_int
  View* SH_this = static_cast<View*>(static_cast<void*>(self));
  SH_this->setScalar<int>(value);
  return;
// splicer end class.View.method.set_scalar_int
}

void SIDRE_view_set_scalar_long(SIDRE_view* self, long value)
{
// splicer begin class.View.method.set_scalar_long
  View* SH_this = static_cast<View*>(static_cast<void*>(self));
  SH_this->setScalar<long>(value);
  return;
// splicer end class.View.method.set_scalar_long
}

void SIDRE_view_set_scalar_float(SIDRE_view* self, float value)
{
// splicer begin class.View.method.set_scalar_float
  View* SH_this = static_cast<View*>(static_cast<void*>(self));
  SH_this->setScalar<float>(value);
  return;
// splicer end class.View.method.set_scalar_float
}

void SIDRE_view_set_scalar_double(SIDRE_view* self, double value)
{
// splicer begin class.View.method.set_scalar_double
  View* SH_this = static_cast<View*>(static_cast<void*>(self));
  SH_this->setScalar<double>(value);
  return;
// splicer end class.View.method.set_scalar_double
}

void SIDRE_view_set_external_data_ptr_only(SIDRE_view* self,
                                           void* external_ptr)
{
// splicer begin class.View.method.set_external_data_ptr_only
  View* SH_this = static_cast<View*>(static_cast<void*>(self));
  SH_this->setExternalDataPtr(external_ptr);
  return;
// splicer end class.View.method.set_external_data_ptr_only
}

void SIDRE_view_set_external_data_ptr_type(SIDRE_view* self, int type,
                                           SIDRE_SidreLength num_elems,
                                           void* external_ptr)
{
// splicer begin class.View.method.set_external_data_ptr_type
  View* SH_this = static_cast<View*>(static_cast<void*>(self));
  SH_this->setExternalDataPtr(getTypeID(type), num_elems, external_ptr);
  return;
// splicer end class.View.method.set_external_data_ptr_type
}

void SIDRE_view_set_string(SIDRE_view* self, const char* value)
{
// splicer begin class.View.method.set_string
  View* SH_this = static_cast<View*>(static_cast<void*>(self));
  const std::string SH_value(value);
  SH_this->setString(SH_value);
  return;
// splicer end class.View.method.set_string
}

void SIDRE_view_set_string_bufferify(SIDRE_view* self, const char* value,
                                     int Lvalue)
{
// splicer begin class.View.method.set_string_bufferify
  View* SH_this = static_cast<View*>(static_cast<void*>(self));
  const std::string SH_value(value, Lvalue);
  SH_this->setString(SH_value);
  return;
// splicer end class.View.method.set_string_bufferify
}

void SIDRE_view_set_external_data_ptr_shape(SIDRE_view* self, int type,
                                            int ndims,
                                            SIDRE_SidreLength* shape,
                                            void* external_ptr)
{
// splicer begin class.View.method.set_external_data_ptr_shape
  View* SH_this = static_cast<View*>(static_cast<void*>(self));
  SH_this->setExternalDataPtr(getTypeID(type), ndims, shape, external_ptr);
  return;
// splicer end class.View.method.set_external_data_ptr_shape
}

const char* SIDRE_view_get_string(SIDRE_view* self)
{
// splicer begin class.View.method.get_string
  View* SH_this = static_cast<View*>(static_cast<void*>(self));
  const char* SH_rv = SH_this->getString();
  return SH_rv;
// splicer end class.View.method.get_string
}

void SIDRE_view_get_string_bufferify(SIDRE_view* self, char* name, int Nname)
{
// splicer begin class.View.method.get_string_bufferify
  View* SH_this = static_cast<View*>(static_cast<void*>(self));
  const char* SH_rv = SH_this->getString();
  if (SH_rv == NULL)
  {
    std::memset(name, ' ', Nname);
  }
  else
  {
    shroud_FccCopy(name, Nname, SH_rv);
  }
  return;
// splicer end class.View.method.get_string_bufferify
}

int SIDRE_view_get_data_int(SIDRE_view* self)
{
// splicer begin class.View.method.get_data_int
  View* SH_this = static_cast<View*>(static_cast<void*>(self));
  int SH_rv = SH_this->getData<int>();
  return SH_rv;
// splicer end class.View.method.get_data_int
}

long SIDRE_view_get_data_long(SIDRE_view* self)
{
// splicer begin class.View.method.get_data_long
  View* SH_this = static_cast<View*>(static_cast<void*>(self));
  long SH_rv = SH_this->getData<long>();
  return SH_rv;
// splicer end class.View.method.get_data_long
}

float SIDRE_view_get_data_float(SIDRE_view* self)
{
// splicer begin class.View.method.get_data_float
  View* SH_this = static_cast<View*>(static_cast<void*>(self));
  float SH_rv = SH_this->getData<float>();
  return SH_rv;
// splicer end class.View.method.get_data_float
}

double SIDRE_view_get_data_double(SIDRE_view* self)
{
// splicer begin class.View.method.get_data_double
  View* SH_this = static_cast<View*>(static_cast<void*>(self));
  double SH_rv = SH_this->getData<double>();
  return SH_rv;
// splicer end class.View.method.get_data_double
}

SIDRE_group* SIDRE_view_get_owning_group(SIDRE_view* self)
{
// splicer begin class.View.method.get_owning_group
  View* SH_this = static_cast<View*>(static_cast<void*>(self));
  Group* SH_rv = SH_this->getOwningGroup();
  SIDRE_group* XSH_rv = static_cast<SIDRE_group*>(static_cast<void*>(SH_rv));
  return XSH_rv;
// splicer end class.View.method.get_owning_group
}

int SIDRE_view_get_type_id(const SIDRE_view* self)
{
// splicer begin class.View.method.get_type_id
  const View* SH_this =
    static_cast<const View*>(static_cast<const void*>(self));
  TypeID SH_rv = SH_this->getTypeID();
  int XSH_rv = static_cast<int>(SH_rv);
  return XSH_rv;
// splicer end class.View.method.get_type_id
}

size_t SIDRE_view_get_total_bytes(const SIDRE_view* self)
{
// splicer begin class.View.method.get_total_bytes
  const View* SH_this =
    static_cast<const View*>(static_cast<const void*>(self));
  size_t SH_rv = SH_this->getTotalBytes();
  return SH_rv;
// splicer end class.View.method.get_total_bytes
}

size_t SIDRE_view_get_bytes_per_element(const SIDRE_view* self)
{
// splicer begin class.View.method.get_bytes_per_element
  const View* SH_this =
    static_cast<const View*>(static_cast<const void*>(self));
  size_t SH_rv = SH_this->getBytesPerElement();
  return SH_rv;
// splicer end class.View.method.get_bytes_per_element
}

size_t SIDRE_view_get_num_elements(const SIDRE_view* self)
{
// splicer begin class.View.method.get_num_elements
  const View* SH_this =
    static_cast<const View*>(static_cast<const void*>(self));
  size_t SH_rv = SH_this->getNumElements();
  return SH_rv;
// splicer end class.View.method.get_num_elements
}

size_t SIDRE_view_get_offset(const SIDRE_view* self)
{
// splicer begin class.View.method.get_offset
  const View* SH_this =
    static_cast<const View*>(static_cast<const void*>(self));
  size_t SH_rv = SH_this->getOffset();
  return SH_rv;
// splicer end class.View.method.get_offset
}

size_t SIDRE_view_get_stride(const SIDRE_view* self)
{
// splicer begin class.View.method.get_stride
  const View* SH_this =
    static_cast<const View*>(static_cast<const void*>(self));
  size_t SH_rv = SH_this->getStride();
  return SH_rv;
// splicer end class.View.method.get_stride
}

int SIDRE_view_get_num_dimensions(const SIDRE_view* self)
{
// splicer begin class.View.method.get_num_dimensions
  const View* SH_this =
    static_cast<const View*>(static_cast<const void*>(self));
  int SH_rv = SH_this->getNumDimensions();
  return SH_rv;
// splicer end class.View.method.get_num_dimensions
}

int SIDRE_view_get_shape(const SIDRE_view* self, int ndims,
                         SIDRE_SidreLength* shape)
{
// splicer begin class.View.method.get_shape
  const View* SH_this =
    static_cast<const View*>(static_cast<const void*>(self));
  int SH_rv = SH_this->getShape(ndims, shape);
  return SH_rv;
// splicer end class.View.method.get_shape
}

bool SIDRE_view_rename(SIDRE_view* self, const char* new_name)
{
// splicer begin class.View.method.rename
  View* SH_this = static_cast<View*>(static_cast<void*>(self));
  const std::string SH_new_name(new_name);
  bool SH_rv = SH_this->rename(SH_new_name);
  return SH_rv;
// splicer end class.View.method.rename
}

bool SIDRE_view_rename_bufferify(SIDRE_view* self, const char* new_name,
                                 int Lnew_name)
{
// splicer begin class.View.method.rename_bufferify
  View* SH_this = static_cast<View*>(static_cast<void*>(self));
  const std::string SH_new_name(new_name, Lnew_name);
  bool SH_rv = SH_this->rename(SH_new_name);
  return SH_rv;
// splicer end class.View.method.rename_bufferify
}

void SIDRE_view_print(const SIDRE_view* self)
{
// splicer begin class.View.method.print
  const View* SH_this =
    static_cast<const View*>(static_cast<const void*>(self));
  SH_this->print();
  return;
// splicer end class.View.method.print
}

}  // extern "C"

}  // namespace sidre
}  // namespace axom

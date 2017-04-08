// wrapDataView.cpp
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
// wrapDataView.cpp
#include "wrapDataView.h"
#include <string>
#include "shroudrt.hpp"
#include "sidre/DataView.hpp"
#include "sidre/SidreTypes.hpp"

extern "C" {
namespace axom
{
namespace sidre
{

void SIDRE_dataview_allocate_simple(SIDRE_dataview * self)
{
  DataView * SH_this = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.allocate_simple
  SH_this->allocate();
  return;
// splicer end class.DataView.method.allocate_simple
}

void SIDRE_dataview_allocate_from_type(SIDRE_dataview * self, int type,
                                       SIDRE_SidreLength num_elems)
{
  DataView * SH_this = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.allocate_from_type
  SH_this->allocate(getTypeID(type), num_elems);
  return;
// splicer end class.DataView.method.allocate_from_type
}

void SIDRE_dataview_reallocate(SIDRE_dataview * self,
                               SIDRE_SidreLength num_elems)
{
  DataView * SH_this = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.reallocate
  SH_this->reallocate(num_elems);
  return;
// splicer end class.DataView.method.reallocate
}

void SIDRE_dataview_apply_0(SIDRE_dataview * self)
{
  DataView * SH_this = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.apply_0
  SH_this->apply();
  return;
// splicer end class.DataView.method.apply_0
}

void SIDRE_dataview_attach_buffer_only(SIDRE_dataview * self,
                                       SIDRE_databuffer * buff)
{
  DataView * SH_this = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.attach_buffer_only
  SH_this->attachBuffer(static_cast<DataBuffer *>(static_cast<void *>(buff)));
  return;
// splicer end class.DataView.method.attach_buffer_only
}

void SIDRE_dataview_attach_buffer_type(SIDRE_dataview * self, int type,
                                       SIDRE_SidreLength num_elems,
                                       SIDRE_databuffer * buff)
{
  DataView * SH_this = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.attach_buffer_type
  SH_this->attachBuffer(getTypeID(type), num_elems,
                        static_cast<DataBuffer *>(static_cast<void *>(buff)));
  return;
// splicer end class.DataView.method.attach_buffer_type
}

void SIDRE_dataview_attach_buffer_shape(SIDRE_dataview * self, int type,
                                        int ndims, SIDRE_SidreLength * shape,
                                        SIDRE_databuffer * buff)
{
  DataView * SH_this = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.attach_buffer_shape
  SH_this->attachBuffer(getTypeID(
                          type), ndims, shape,
                        static_cast<DataBuffer *>(static_cast<void *>(buff)));
  return;
// splicer end class.DataView.method.attach_buffer_shape
}

void SIDRE_dataview_apply_nelems(SIDRE_dataview * self,
                                 SIDRE_SidreLength num_elems)
{
  DataView * SH_this = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.apply_nelems
  SH_this->apply(num_elems);
  return;
// splicer end class.DataView.method.apply_nelems
}

void SIDRE_dataview_apply_nelems_offset(SIDRE_dataview * self,
                                        SIDRE_SidreLength num_elems,
                                        SIDRE_SidreLength offset)
{
  DataView * SH_this = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.apply_nelems_offset
  SH_this->apply(num_elems, offset);
  return;
// splicer end class.DataView.method.apply_nelems_offset
}

void SIDRE_dataview_apply_nelems_offset_stride(SIDRE_dataview * self,
                                               SIDRE_SidreLength num_elems,
                                               SIDRE_SidreLength offset,
                                               SIDRE_SidreLength stride)
{
  DataView * SH_this = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.apply_nelems_offset_stride
  SH_this->apply(num_elems, offset, stride);
  return;
// splicer end class.DataView.method.apply_nelems_offset_stride
}

void SIDRE_dataview_apply_type_nelems(SIDRE_dataview * self, int type,
                                      SIDRE_SidreLength num_elems)
{
  DataView * SH_this = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.apply_type_nelems
  SH_this->apply(getTypeID(type), num_elems);
  return;
// splicer end class.DataView.method.apply_type_nelems
}

void SIDRE_dataview_apply_type_nelems_offset(SIDRE_dataview * self, int type,
                                             SIDRE_SidreLength num_elems,
                                             SIDRE_SidreLength offset)
{
  DataView * SH_this = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.apply_type_nelems_offset
  SH_this->apply(getTypeID(type), num_elems, offset);
  return;
// splicer end class.DataView.method.apply_type_nelems_offset
}

void SIDRE_dataview_apply_type_nelems_offset_stride(SIDRE_dataview * self,
                                                    int type,
                                                    SIDRE_SidreLength num_elems,
                                                    SIDRE_SidreLength offset,
                                                    SIDRE_SidreLength stride)
{
  DataView * SH_this = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.apply_type_nelems_offset_stride
  SH_this->apply(getTypeID(type), num_elems, offset, stride);
  return;
// splicer end class.DataView.method.apply_type_nelems_offset_stride
}

void SIDRE_dataview_apply_type_shape(SIDRE_dataview * self, int type, int ndims,
                                     SIDRE_SidreLength * shape)
{
  DataView * SH_this = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.apply_type_shape
  SH_this->apply(getTypeID(type), ndims, shape);
  return;
// splicer end class.DataView.method.apply_type_shape
}

bool SIDRE_dataview_has_buffer(const SIDRE_dataview * self)
{
  const DataView * SH_this =
    static_cast<const DataView *>(static_cast<const void *>(self));
// splicer begin class.DataView.method.has_buffer
  bool SH_rv = SH_this->hasBuffer();
  return SH_rv;
// splicer end class.DataView.method.has_buffer
}

bool SIDRE_dataview_is_external(const SIDRE_dataview * self)
{
  const DataView * SH_this =
    static_cast<const DataView *>(static_cast<const void *>(self));
// splicer begin class.DataView.method.is_external
  bool SH_rv = SH_this->isExternal();
  return SH_rv;
// splicer end class.DataView.method.is_external
}

bool SIDRE_dataview_is_allocated(SIDRE_dataview * self)
{
  DataView * SH_this = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.is_allocated
  bool SH_rv = SH_this->isAllocated();
  return SH_rv;
// splicer end class.DataView.method.is_allocated
}

bool SIDRE_dataview_is_applied(const SIDRE_dataview * self)
{
  const DataView * SH_this =
    static_cast<const DataView *>(static_cast<const void *>(self));
// splicer begin class.DataView.method.is_applied
  bool SH_rv = SH_this->isApplied();
  return SH_rv;
// splicer end class.DataView.method.is_applied
}

bool SIDRE_dataview_is_described(const SIDRE_dataview * self)
{
  const DataView * SH_this =
    static_cast<const DataView *>(static_cast<const void *>(self));
// splicer begin class.DataView.method.is_described
  bool SH_rv = SH_this->isDescribed();
  return SH_rv;
// splicer end class.DataView.method.is_described
}

bool SIDRE_dataview_is_empty(const SIDRE_dataview * self)
{
  const DataView * SH_this =
    static_cast<const DataView *>(static_cast<const void *>(self));
// splicer begin class.DataView.method.is_empty
  bool SH_rv = SH_this->isEmpty();
  return SH_rv;
// splicer end class.DataView.method.is_empty
}

bool SIDRE_dataview_is_opaque(const SIDRE_dataview * self)
{
  const DataView * SH_this =
    static_cast<const DataView *>(static_cast<const void *>(self));
// splicer begin class.DataView.method.is_opaque
  bool SH_rv = SH_this->isOpaque();
  return SH_rv;
// splicer end class.DataView.method.is_opaque
}

bool SIDRE_dataview_is_scalar(const SIDRE_dataview * self)
{
  const DataView * SH_this =
    static_cast<const DataView *>(static_cast<const void *>(self));
// splicer begin class.DataView.method.is_scalar
  bool SH_rv = SH_this->isScalar();
  return SH_rv;
// splicer end class.DataView.method.is_scalar
}

bool SIDRE_dataview_is_string(const SIDRE_dataview * self)
{
  const DataView * SH_this =
    static_cast<const DataView *>(static_cast<const void *>(self));
// splicer begin class.DataView.method.is_string
  bool SH_rv = SH_this->isString();
  return SH_rv;
// splicer end class.DataView.method.is_string
}

const char * SIDRE_dataview_get_name(const SIDRE_dataview * self)
{
  const DataView * SH_this =
    static_cast<const DataView *>(static_cast<const void *>(self));
// splicer begin class.DataView.method.get_name
  const std::string & SH_rv = SH_this->getName();
  return SH_rv.c_str();
// splicer end class.DataView.method.get_name
}

void SIDRE_dataview_get_name_bufferify(const SIDRE_dataview * self,
                                       char * SH_F_rv, int LSH_F_rv)
{
  const DataView * SH_this =
    static_cast<const DataView *>(static_cast<const void *>(self));
// splicer begin class.DataView.method.get_name_bufferify
  const std::string & SH_rv = SH_this->getName();
  shroud_FccCopy(SH_F_rv, LSH_F_rv, SH_rv.c_str());
  return;
// splicer end class.DataView.method.get_name_bufferify
}

SIDRE_databuffer * SIDRE_dataview_get_buffer(SIDRE_dataview * self)
{
  DataView * SH_this = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.get_buffer
  DataBuffer * SH_rv = SH_this->getBuffer();
  return static_cast<SIDRE_databuffer *>(static_cast<void *>(SH_rv));
// splicer end class.DataView.method.get_buffer
}

void * SIDRE_dataview_get_void_ptr(const SIDRE_dataview * self)
{
  const DataView * SH_this =
    static_cast<const DataView *>(static_cast<const void *>(self));
// splicer begin class.DataView.method.get_void_ptr
  void * SH_rv = SH_this->getVoidPtr();
  return SH_rv;
// splicer end class.DataView.method.get_void_ptr
}

void SIDRE_dataview_set_scalar_int(SIDRE_dataview * self, int value)
{
  DataView * SH_this = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.set_scalar_int
  SH_this->setScalar<int>(value);
  return;
// splicer end class.DataView.method.set_scalar_int
}

void SIDRE_dataview_set_scalar_long(SIDRE_dataview * self, long value)
{
  DataView * SH_this = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.set_scalar_long
  SH_this->setScalar<long>(value);
  return;
// splicer end class.DataView.method.set_scalar_long
}

void SIDRE_dataview_set_scalar_float(SIDRE_dataview * self, float value)
{
  DataView * SH_this = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.set_scalar_float
  SH_this->setScalar<float>(value);
  return;
// splicer end class.DataView.method.set_scalar_float
}

void SIDRE_dataview_set_scalar_double(SIDRE_dataview * self, double value)
{
  DataView * SH_this = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.set_scalar_double
  SH_this->setScalar<double>(value);
  return;
// splicer end class.DataView.method.set_scalar_double
}

void SIDRE_dataview_set_external_data_ptr_only(SIDRE_dataview * self,
                                               void * external_ptr)
{
  DataView * SH_this = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.set_external_data_ptr_only
  SH_this->setExternalDataPtr(external_ptr);
  return;
// splicer end class.DataView.method.set_external_data_ptr_only
}

void SIDRE_dataview_set_external_data_ptr_type(SIDRE_dataview * self, int type,
                                               SIDRE_SidreLength num_elems,
                                               void * external_ptr)
{
  DataView * SH_this = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.set_external_data_ptr_type
  SH_this->setExternalDataPtr(getTypeID(type), num_elems, external_ptr);
  return;
// splicer end class.DataView.method.set_external_data_ptr_type
}

void SIDRE_dataview_set_string(SIDRE_dataview * self, const char * value)
{
  DataView * SH_this = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.set_string
  const std::string SH_value(value);
  SH_this->setString(SH_value);
  return;
// splicer end class.DataView.method.set_string
}

void SIDRE_dataview_set_string_bufferify(SIDRE_dataview * self,
                                         const char * value, int Lvalue)
{
  DataView * SH_this = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.set_string_bufferify
  const std::string SH_value(value, Lvalue);
  SH_this->setString(SH_value);
  return;
// splicer end class.DataView.method.set_string_bufferify
}

void SIDRE_dataview_set_external_data_ptr_shape(SIDRE_dataview * self, int type,
                                                int ndims,
                                                SIDRE_SidreLength * shape,
                                                void * external_ptr)
{
  DataView * SH_this = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.set_external_data_ptr_shape
  SH_this->setExternalDataPtr(getTypeID(type), ndims, shape, external_ptr);
  return;
// splicer end class.DataView.method.set_external_data_ptr_shape
}

const char * SIDRE_dataview_get_string(SIDRE_dataview * self)
{
  DataView * SH_this = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.get_string
  const char * SH_rv = SH_this->getString();
  return SH_rv;
// splicer end class.DataView.method.get_string
}

void SIDRE_dataview_get_string_bufferify(SIDRE_dataview * self, char * name,
                                         int Lname)
{
  DataView * SH_this = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.get_string_bufferify
  const char * SH_rv = SH_this->getString();
  shroud_FccCopy(name, Lname, SH_rv);
  return;
// splicer end class.DataView.method.get_string_bufferify
}

int SIDRE_dataview_get_data_int(SIDRE_dataview * self)
{
  DataView * SH_this = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.get_data_int
  int SH_rv = SH_this->getData<int>();
  return SH_rv;
// splicer end class.DataView.method.get_data_int
}

long SIDRE_dataview_get_data_long(SIDRE_dataview * self)
{
  DataView * SH_this = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.get_data_long
  long SH_rv = SH_this->getData<long>();
  return SH_rv;
// splicer end class.DataView.method.get_data_long
}

float SIDRE_dataview_get_data_float(SIDRE_dataview * self)
{
  DataView * SH_this = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.get_data_float
  float SH_rv = SH_this->getData<float>();
  return SH_rv;
// splicer end class.DataView.method.get_data_float
}

double SIDRE_dataview_get_data_double(SIDRE_dataview * self)
{
  DataView * SH_this = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.get_data_double
  double SH_rv = SH_this->getData<double>();
  return SH_rv;
// splicer end class.DataView.method.get_data_double
}

SIDRE_datagroup * SIDRE_dataview_get_owning_group(SIDRE_dataview * self)
{
  DataView * SH_this = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.get_owning_group
  DataGroup * SH_rv = SH_this->getOwningGroup();
  return static_cast<SIDRE_datagroup *>(static_cast<void *>(SH_rv));
// splicer end class.DataView.method.get_owning_group
}

int SIDRE_dataview_get_type_id(const SIDRE_dataview * self)
{
  const DataView * SH_this =
    static_cast<const DataView *>(static_cast<const void *>(self));
// splicer begin class.DataView.method.get_type_id
  TypeID SH_rv = SH_this->getTypeID();
  return static_cast<int>(SH_rv);
// splicer end class.DataView.method.get_type_id
}

size_t SIDRE_dataview_get_total_bytes(const SIDRE_dataview * self)
{
  const DataView * SH_this =
    static_cast<const DataView *>(static_cast<const void *>(self));
// splicer begin class.DataView.method.get_total_bytes
  size_t SH_rv = SH_this->getTotalBytes();
  return SH_rv;
// splicer end class.DataView.method.get_total_bytes
}

size_t SIDRE_dataview_get_bytes_per_element(const SIDRE_dataview * self)
{
  const DataView * SH_this =
    static_cast<const DataView *>(static_cast<const void *>(self));
// splicer begin class.DataView.method.get_bytes_per_element
  size_t SH_rv = SH_this->getBytesPerElement();
  return SH_rv;
// splicer end class.DataView.method.get_bytes_per_element
}

size_t SIDRE_dataview_get_num_elements(const SIDRE_dataview * self)
{
  const DataView * SH_this =
    static_cast<const DataView *>(static_cast<const void *>(self));
// splicer begin class.DataView.method.get_num_elements
  size_t SH_rv = SH_this->getNumElements();
  return SH_rv;
// splicer end class.DataView.method.get_num_elements
}

size_t SIDRE_dataview_get_offset(const SIDRE_dataview * self)
{
  const DataView * SH_this =
    static_cast<const DataView *>(static_cast<const void *>(self));
// splicer begin class.DataView.method.get_offset
  size_t SH_rv = SH_this->getOffset();
  return SH_rv;
// splicer end class.DataView.method.get_offset
}

size_t SIDRE_dataview_get_stride(const SIDRE_dataview * self)
{
  const DataView * SH_this =
    static_cast<const DataView *>(static_cast<const void *>(self));
// splicer begin class.DataView.method.get_stride
  size_t SH_rv = SH_this->getStride();
  return SH_rv;
// splicer end class.DataView.method.get_stride
}

int SIDRE_dataview_get_num_dimensions(const SIDRE_dataview * self)
{
  const DataView * SH_this =
    static_cast<const DataView *>(static_cast<const void *>(self));
// splicer begin class.DataView.method.get_num_dimensions
  int SH_rv = SH_this->getNumDimensions();
  return SH_rv;
// splicer end class.DataView.method.get_num_dimensions
}

int SIDRE_dataview_get_shape(const SIDRE_dataview * self, int ndims,
                             SIDRE_SidreLength * shape)
{
  const DataView * SH_this =
    static_cast<const DataView *>(static_cast<const void *>(self));
// splicer begin class.DataView.method.get_shape
  int SH_rv = SH_this->getShape(ndims, shape);
  return SH_rv;
// splicer end class.DataView.method.get_shape
}

bool SIDRE_dataview_rename(SIDRE_dataview * self, const char * new_name)
{
  DataView * SH_this = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.rename
  const std::string SH_new_name(new_name);
  bool SH_rv = SH_this->rename(SH_new_name);
  return SH_rv;
// splicer end class.DataView.method.rename
}

bool SIDRE_dataview_rename_bufferify(SIDRE_dataview * self,
                                     const char * new_name, int Lnew_name)
{
  DataView * SH_this = static_cast<DataView *>(static_cast<void *>(self));
// splicer begin class.DataView.method.rename_bufferify
  const std::string SH_new_name(new_name, Lnew_name);
  bool SH_rv = SH_this->rename(SH_new_name);
  return SH_rv;
// splicer end class.DataView.method.rename_bufferify
}

void SIDRE_dataview_print(const SIDRE_dataview * self)
{
  const DataView * SH_this =
    static_cast<const DataView *>(static_cast<const void *>(self));
// splicer begin class.DataView.method.print
  SH_this->print();
  return;
// splicer end class.DataView.method.print
}

// splicer begin class.DataView.additional_functions
// splicer end class.DataView.additional_functions

}  // namespace axom
}  // namespace sidre
}  // extern "C"

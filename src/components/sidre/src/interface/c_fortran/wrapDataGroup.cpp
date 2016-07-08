// wrapDataGroup.cpp
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
// wrapDataGroup.cpp
#include "wrapDataGroup.h"
#include <string>
#include "shroudrt.hpp"
#include "sidre/DataGroup.hpp"
#include "sidre/SidreTypes.hpp"

extern "C" {
namespace asctoolkit
{
namespace sidre
{

const char * SIDRE_datagroup_get_name(const SIDRE_datagroup * self)
{
  const DataGroup * selfobj =
    static_cast<const DataGroup *>(static_cast<const void *>(self));
// splicer begin class.DataGroup.method.get_name
  const std::string & rv = selfobj->getName();
  return rv.c_str();
// splicer end class.DataGroup.method.get_name
}

void SIDRE_datagroup_get_name_bufferify(const SIDRE_datagroup * self,
                                        char * SH_F_rv, int LSH_F_rv)
{
  const DataGroup * selfobj =
    static_cast<const DataGroup *>(static_cast<const void *>(self));
// splicer begin class.DataGroup.method.get_name_bufferify
  const std::string & rv = selfobj->getName();
  asctoolkit::shroud::FccCopy(SH_F_rv, LSH_F_rv, rv.c_str());
  return;
// splicer end class.DataGroup.method.get_name_bufferify
}

const SIDRE_datagroup * SIDRE_datagroup_get_parent(const SIDRE_datagroup * self)
{
  const DataGroup * selfobj =
    static_cast<const DataGroup *>(static_cast<const void *>(self));
// splicer begin class.DataGroup.method.get_parent
  const DataGroup * rv = selfobj->getParent();
  return static_cast<const SIDRE_datagroup *>(static_cast<const void *>(rv));
// splicer end class.DataGroup.method.get_parent
}

const SIDRE_datastore * SIDRE_datagroup_get_data_store(
  const SIDRE_datagroup * self)
{
  const DataGroup * selfobj =
    static_cast<const DataGroup *>(static_cast<const void *>(self));
// splicer begin class.DataGroup.method.get_data_store
  const DataStore * rv = selfobj->getDataStore();
  return static_cast<const SIDRE_datastore *>(static_cast<const void *>(rv));
// splicer end class.DataGroup.method.get_data_store
}

size_t SIDRE_datagroup_get_num_views(const SIDRE_datagroup * self)
{
  const DataGroup * selfobj =
    static_cast<const DataGroup *>(static_cast<const void *>(self));
// splicer begin class.DataGroup.method.get_num_views
  size_t rv = selfobj->getNumViews();
  return rv;
// splicer end class.DataGroup.method.get_num_views
}

size_t SIDRE_datagroup_get_num_groups(const SIDRE_datagroup * self)
{
  const DataGroup * selfobj =
    static_cast<const DataGroup *>(static_cast<const void *>(self));
// splicer begin class.DataGroup.method.get_num_groups
  size_t rv = selfobj->getNumGroups();
  return rv;
// splicer end class.DataGroup.method.get_num_groups
}

bool SIDRE_datagroup_has_view(SIDRE_datagroup * self, const char * name)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.has_view
  const std::string SH_name(name);
  bool rv = selfobj->hasView(SH_name);
  return rv;
// splicer end class.DataGroup.method.has_view
}

bool SIDRE_datagroup_has_view_bufferify(SIDRE_datagroup * self,
                                        const char * name, int Lname)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.has_view_bufferify
  const std::string SH_name(name, Lname);
  bool rv = selfobj->hasView(SH_name);
  return rv;
// splicer end class.DataGroup.method.has_view_bufferify
}

SIDRE_dataview * SIDRE_datagroup_get_view_from_name(SIDRE_datagroup * self,
                                                    const char * name)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.get_view_from_name
  const std::string SH_name(name);
  DataView * rv = selfobj->getView(SH_name);
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.get_view_from_name
}

SIDRE_dataview * SIDRE_datagroup_get_view_from_name_bufferify(
  SIDRE_datagroup * self, const char * name, int Lname)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.get_view_from_name_bufferify
  const std::string SH_name(name, Lname);
  DataView * rv = selfobj->getView(SH_name);
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.get_view_from_name_bufferify
}

SIDRE_dataview * SIDRE_datagroup_get_view_from_index(SIDRE_datagroup * self,
                                                     const SIDRE_IndexType idx)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.get_view_from_index
  DataView * rv = selfobj->getView(idx);
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.get_view_from_index
}

SIDRE_IndexType SIDRE_datagroup_get_view_index(const SIDRE_datagroup * self,
                                               const char * name)
{
  const DataGroup * selfobj =
    static_cast<const DataGroup *>(static_cast<const void *>(self));
// splicer begin class.DataGroup.method.get_view_index
  const std::string SH_name(name);
  IndexType rv = selfobj->getViewIndex(SH_name);
  return rv;
// splicer end class.DataGroup.method.get_view_index
}

SIDRE_IndexType SIDRE_datagroup_get_view_index_bufferify(
  const SIDRE_datagroup * self, const char * name, int Lname)
{
  const DataGroup * selfobj =
    static_cast<const DataGroup *>(static_cast<const void *>(self));
// splicer begin class.DataGroup.method.get_view_index_bufferify
  const std::string SH_name(name, Lname);
  IndexType rv = selfobj->getViewIndex(SH_name);
  return rv;
// splicer end class.DataGroup.method.get_view_index_bufferify
}

const char * SIDRE_datagroup_get_view_name(const SIDRE_datagroup * self,
                                           SIDRE_IndexType idx)
{
  const DataGroup * selfobj =
    static_cast<const DataGroup *>(static_cast<const void *>(self));
// splicer begin class.DataGroup.method.get_view_name
  const std::string & rv = selfobj->getViewName(idx);
// check for error
  if (!nameIsValid(rv))
  {
    return SIDRE_InvalidName;
  }

  return rv.c_str();
// splicer end class.DataGroup.method.get_view_name
}

void SIDRE_datagroup_get_view_name_bufferify(const SIDRE_datagroup * self,
                                             SIDRE_IndexType idx,
                                             char * SH_F_rv, int LSH_F_rv)
{
  const DataGroup * selfobj =
    static_cast<const DataGroup *>(static_cast<const void *>(self));
// splicer begin class.DataGroup.method.get_view_name_bufferify
  const std::string & rv = selfobj->getViewName(idx);
// check for error
  if (!nameIsValid(rv))
  {
    std::memset(SH_F_rv, ' ', LSH_F_rv);
    return;
  }

  asctoolkit::shroud::FccCopy(SH_F_rv, LSH_F_rv, rv.c_str());
  return;
// splicer end class.DataGroup.method.get_view_name_bufferify
}

SIDRE_IndexType SIDRE_datagroup_get_first_valid_view_index(
  const SIDRE_datagroup * self)
{
  const DataGroup * selfobj =
    static_cast<const DataGroup *>(static_cast<const void *>(self));
// splicer begin class.DataGroup.method.get_first_valid_view_index
  IndexType rv = selfobj->getFirstValidViewIndex();
  return rv;
// splicer end class.DataGroup.method.get_first_valid_view_index
}

SIDRE_IndexType SIDRE_datagroup_get_next_valid_view_index(
  const SIDRE_datagroup * self, SIDRE_IndexType idx)
{
  const DataGroup * selfobj =
    static_cast<const DataGroup *>(static_cast<const void *>(self));
// splicer begin class.DataGroup.method.get_next_valid_view_index
  IndexType rv = selfobj->getNextValidViewIndex(idx);
  return rv;
// splicer end class.DataGroup.method.get_next_valid_view_index
}

SIDRE_dataview * SIDRE_datagroup_create_view_and_allocate_nelems(
  SIDRE_datagroup * self, const char * name, int type,
  SIDRE_SidreLength num_elems)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_and_allocate_nelems
  const std::string SH_name(name);
  DataView * rv = selfobj->createViewAndAllocate(SH_name, getTypeID(
                                                   type), num_elems);
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_and_allocate_nelems
}

SIDRE_dataview * SIDRE_datagroup_create_view_and_allocate_nelems_bufferify(
  SIDRE_datagroup * self, const char * name, int Lname, int type,
  SIDRE_SidreLength num_elems)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_and_allocate_nelems_bufferify
  const std::string SH_name(name, Lname);
  DataView * rv = selfobj->createViewAndAllocate(SH_name, getTypeID(
                                                   type), num_elems);
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_and_allocate_nelems_bufferify
}

SIDRE_dataview * SIDRE_datagroup_create_view_and_allocate_shape(
  SIDRE_datagroup * self, const char * name, int type, int ndims,
  SIDRE_SidreLength * shape)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_and_allocate_shape
  const std::string SH_name(name);
  DataView * rv = selfobj->createViewAndAllocate(SH_name, getTypeID(
                                                   type), ndims, shape);
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_and_allocate_shape
}

SIDRE_dataview * SIDRE_datagroup_create_view_and_allocate_shape_bufferify(
  SIDRE_datagroup * self, const char * name, int Lname, int type,
  int ndims, SIDRE_SidreLength * shape)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_and_allocate_shape_bufferify
  const std::string SH_name(name, Lname);
  DataView * rv = selfobj->createViewAndAllocate(SH_name, getTypeID(
                                                   type), ndims, shape);
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_and_allocate_shape_bufferify
}

SIDRE_dataview * SIDRE_datagroup_create_view_scalar_int(SIDRE_datagroup * self,
                                                        const char * name,
                                                        int value)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_scalar_int
  const std::string SH_name(name);
  DataView * rv = selfobj->createViewScalar<int>(SH_name, value);
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_scalar_int
}

SIDRE_dataview * SIDRE_datagroup_create_view_scalar_int_bufferify(
  SIDRE_datagroup * self, const char * name, int Lname, int value)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_scalar_int_bufferify
  const std::string SH_name(name, Lname);
  DataView * rv = selfobj->createViewScalar<int>(SH_name, value);
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_scalar_int_bufferify
}

SIDRE_dataview * SIDRE_datagroup_create_view_scalar_long(SIDRE_datagroup * self,
                                                         const char * name,
                                                         long value)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_scalar_long
  const std::string SH_name(name);
  DataView * rv = selfobj->createViewScalar<long>(SH_name, value);
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_scalar_long
}

SIDRE_dataview * SIDRE_datagroup_create_view_scalar_long_bufferify(
  SIDRE_datagroup * self, const char * name, int Lname, long value)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_scalar_long_bufferify
  const std::string SH_name(name, Lname);
  DataView * rv = selfobj->createViewScalar<long>(SH_name, value);
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_scalar_long_bufferify
}

SIDRE_dataview * SIDRE_datagroup_create_view_scalar_float(
  SIDRE_datagroup * self, const char * name, float value)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_scalar_float
  const std::string SH_name(name);
  DataView * rv = selfobj->createViewScalar<float>(SH_name, value);
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_scalar_float
}

SIDRE_dataview * SIDRE_datagroup_create_view_scalar_float_bufferify(
  SIDRE_datagroup * self, const char * name, int Lname, float value)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_scalar_float_bufferify
  const std::string SH_name(name, Lname);
  DataView * rv = selfobj->createViewScalar<float>(SH_name, value);
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_scalar_float_bufferify
}

SIDRE_dataview * SIDRE_datagroup_create_view_scalar_double(
  SIDRE_datagroup * self, const char * name, double value)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_scalar_double
  const std::string SH_name(name);
  DataView * rv = selfobj->createViewScalar<double>(SH_name, value);
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_scalar_double
}

SIDRE_dataview * SIDRE_datagroup_create_view_scalar_double_bufferify(
  SIDRE_datagroup * self, const char * name, int Lname, double value)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_scalar_double_bufferify
  const std::string SH_name(name, Lname);
  DataView * rv = selfobj->createViewScalar<double>(SH_name, value);
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_scalar_double_bufferify
}

SIDRE_dataview * SIDRE_datagroup_create_view_string(SIDRE_datagroup * self,
                                                    const char * name,
                                                    const char * value)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_string
  const std::string SH_name(name);
  const std::string SH_value(value);
  DataView * rv = selfobj->createViewString(SH_name, SH_value);
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_string
}

SIDRE_dataview * SIDRE_datagroup_create_view_string_bufferify(
  SIDRE_datagroup * self, const char * name, int Lname,
  const char * value, int Lvalue)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_string_bufferify
  const std::string SH_name(name, Lname);
  const std::string SH_value(value, Lvalue);
  DataView * rv = selfobj->createViewString(SH_name, SH_value);
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_string_bufferify
}

SIDRE_dataview * SIDRE_datagroup_create_view_empty(SIDRE_datagroup * self,
                                                   const char * name)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_empty
  const std::string SH_name(name);
  DataView * rv = selfobj->createView(SH_name);
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_empty
}

SIDRE_dataview * SIDRE_datagroup_create_view_empty_bufferify(
  SIDRE_datagroup * self, const char * name, int Lname)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_empty_bufferify
  const std::string SH_name(name, Lname);
  DataView * rv = selfobj->createView(SH_name);
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_empty_bufferify
}

SIDRE_dataview * SIDRE_datagroup_create_view_from_type(SIDRE_datagroup * self,
                                                       const char * name,
                                                       int type,
                                                       SIDRE_SidreLength num_elems)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_from_type
  const std::string SH_name(name);
  DataView * rv = selfobj->createView(SH_name, getTypeID(type), num_elems);
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_from_type
}

SIDRE_dataview * SIDRE_datagroup_create_view_from_type_bufferify(
  SIDRE_datagroup * self, const char * name, int Lname, int type,
  SIDRE_SidreLength num_elems)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_from_type_bufferify
  const std::string SH_name(name, Lname);
  DataView * rv = selfobj->createView(SH_name, getTypeID(type), num_elems);
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_from_type_bufferify
}

SIDRE_dataview * SIDRE_datagroup_create_view_from_type_and_buffer(
  SIDRE_datagroup * self, const char * name, int type,
  SIDRE_SidreLength num_elems, SIDRE_databuffer * buff)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_from_type_and_buffer
  const std::string SH_name(name);
  DataView * rv = selfobj->createView(SH_name, getTypeID(
                                        type), num_elems,
                                      static_cast<DataBuffer *>(static_cast<void
                                                                            *>(
                                                                  buff)));
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_from_type_and_buffer
}

SIDRE_dataview * SIDRE_datagroup_create_view_from_type_and_buffer_bufferify(
  SIDRE_datagroup * self, const char * name, int Lname, int type,
  SIDRE_SidreLength num_elems, SIDRE_databuffer * buff)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_from_type_and_buffer_bufferify
  const std::string SH_name(name, Lname);
  DataView * rv = selfobj->createView(SH_name, getTypeID(
                                        type), num_elems,
                                      static_cast<DataBuffer *>(static_cast<void
                                                                            *>(
                                                                  buff)));
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_from_type_and_buffer_bufferify
}

SIDRE_dataview * SIDRE_datagroup_create_view_from_type_external(
  SIDRE_datagroup * self, const char * name, int type,
  SIDRE_SidreLength num_elems, void * external_ptr)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_from_type_external
  const std::string SH_name(name);
  DataView * rv = selfobj->createView(SH_name, getTypeID(
                                        type), num_elems, external_ptr);
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_from_type_external
}

SIDRE_dataview * SIDRE_datagroup_create_view_from_type_external_bufferify(
  SIDRE_datagroup * self, const char * name, int Lname, int type,
  SIDRE_SidreLength num_elems, void * external_ptr)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_from_type_external_bufferify
  const std::string SH_name(name, Lname);
  DataView * rv = selfobj->createView(SH_name, getTypeID(
                                        type), num_elems, external_ptr);
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_from_type_external_bufferify
}

SIDRE_dataview * SIDRE_datagroup_create_view_from_shape(SIDRE_datagroup * self,
                                                        const char * name,
                                                        int type, int ndims,
                                                        SIDRE_SidreLength * shape)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_from_shape
  const std::string SH_name(name);
  DataView * rv = selfobj->createView(SH_name, getTypeID(type), ndims, shape);
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_from_shape
}

SIDRE_dataview * SIDRE_datagroup_create_view_from_shape_bufferify(
  SIDRE_datagroup * self, const char * name, int Lname, int type,
  int ndims, SIDRE_SidreLength * shape)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_from_shape_bufferify
  const std::string SH_name(name, Lname);
  DataView * rv = selfobj->createView(SH_name, getTypeID(type), ndims, shape);
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_from_shape_bufferify
}

SIDRE_dataview * SIDRE_datagroup_create_view_from_shape_and_buffer(
  SIDRE_datagroup * self, const char * name, int type, int ndims,
  SIDRE_SidreLength * shape, SIDRE_databuffer * buff)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_from_shape_and_buffer
  const std::string SH_name(name);
  DataView * rv = selfobj->createView(SH_name, getTypeID(
                                        type), ndims, shape,
                                      static_cast<DataBuffer *>(static_cast<void
                                                                            *>(
                                                                  buff)));
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_from_shape_and_buffer
}

SIDRE_dataview * SIDRE_datagroup_create_view_from_shape_and_buffer_bufferify(
  SIDRE_datagroup * self, const char * name, int Lname, int type,
  int ndims, SIDRE_SidreLength * shape, SIDRE_databuffer * buff)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_from_shape_and_buffer_bufferify
  const std::string SH_name(name, Lname);
  DataView * rv = selfobj->createView(SH_name, getTypeID(
                                        type), ndims, shape,
                                      static_cast<DataBuffer *>(static_cast<void
                                                                            *>(
                                                                  buff)));
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_from_shape_and_buffer_bufferify
}

SIDRE_dataview * SIDRE_datagroup_create_view_from_shape_external(
  SIDRE_datagroup * self, const char * name, int type, int ndims,
  SIDRE_SidreLength * shape, void * external_ptr)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_from_shape_external
  const std::string SH_name(name);
  DataView * rv = selfobj->createView(SH_name, getTypeID(
                                        type), ndims, shape, external_ptr);
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_from_shape_external
}

SIDRE_dataview * SIDRE_datagroup_create_view_from_shape_external_bufferify(
  SIDRE_datagroup * self, const char * name, int Lname, int type,
  int ndims, SIDRE_SidreLength * shape, void * external_ptr)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_from_shape_external_bufferify
  const std::string SH_name(name, Lname);
  DataView * rv = selfobj->createView(SH_name, getTypeID(
                                        type), ndims, shape, external_ptr);
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_from_shape_external_bufferify
}

SIDRE_dataview * SIDRE_datagroup_create_view_into_buffer(SIDRE_datagroup * self,
                                                         const char * name,
                                                         SIDRE_databuffer * buff)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_into_buffer
  const std::string SH_name(name);
  DataView * rv =
    selfobj->createView(SH_name,
                        static_cast<DataBuffer *>(static_cast<void *>(buff)));
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_into_buffer
}

SIDRE_dataview * SIDRE_datagroup_create_view_into_buffer_bufferify(
  SIDRE_datagroup * self, const char * name, int Lname,
  SIDRE_databuffer * buff)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_into_buffer_bufferify
  const std::string SH_name(name, Lname);
  DataView * rv =
    selfobj->createView(SH_name,
                        static_cast<DataBuffer *>(static_cast<void *>(buff)));
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_into_buffer_bufferify
}

SIDRE_dataview * SIDRE_datagroup_create_view_external(SIDRE_datagroup * self,
                                                      const char * name,
                                                      void * external_ptr)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_external
  const std::string SH_name(name);
  DataView * rv = selfobj->createView(SH_name, external_ptr);
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_external
}

SIDRE_dataview * SIDRE_datagroup_create_view_external_bufferify(
  SIDRE_datagroup * self, const char * name, int Lname,
  void * external_ptr)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_external_bufferify
  const std::string SH_name(name, Lname);
  DataView * rv = selfobj->createView(SH_name, external_ptr);
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_external_bufferify
}

void SIDRE_datagroup_destroy_view(SIDRE_datagroup * self, const char * name)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.destroy_view
  const std::string SH_name(name);
  selfobj->destroyView(SH_name);
  return;
// splicer end class.DataGroup.method.destroy_view
}

void SIDRE_datagroup_destroy_view_bufferify(SIDRE_datagroup * self,
                                            const char * name, int Lname)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.destroy_view_bufferify
  const std::string SH_name(name, Lname);
  selfobj->destroyView(SH_name);
  return;
// splicer end class.DataGroup.method.destroy_view_bufferify
}

void SIDRE_datagroup_destroy_view_and_data_name(SIDRE_datagroup * self,
                                                const char * name)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.destroy_view_and_data_name
  const std::string SH_name(name);
  selfobj->destroyViewAndData(SH_name);
  return;
// splicer end class.DataGroup.method.destroy_view_and_data_name
}

void SIDRE_datagroup_destroy_view_and_data_name_bufferify(
  SIDRE_datagroup * self, const char * name, int Lname)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.destroy_view_and_data_name_bufferify
  const std::string SH_name(name, Lname);
  selfobj->destroyViewAndData(SH_name);
  return;
// splicer end class.DataGroup.method.destroy_view_and_data_name_bufferify
}

void SIDRE_datagroup_destroy_view_and_data_index(SIDRE_datagroup * self,
                                                 SIDRE_IndexType idx)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.destroy_view_and_data_index
  selfobj->destroyViewAndData(idx);
  return;
// splicer end class.DataGroup.method.destroy_view_and_data_index
}

SIDRE_dataview * SIDRE_datagroup_move_view(SIDRE_datagroup * self,
                                           SIDRE_dataview * view)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.move_view
  DataView * rv =
    selfobj->moveView(static_cast<DataView *>(static_cast<void *>(view)));
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.move_view
}

SIDRE_dataview * SIDRE_datagroup_copy_view(SIDRE_datagroup * self,
                                           SIDRE_dataview * view)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.copy_view
  DataView * rv =
    selfobj->copyView(static_cast<DataView *>(static_cast<void *>(view)));
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.copy_view
}

bool SIDRE_datagroup_has_group(SIDRE_datagroup * self, const char * name)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.has_group
  const std::string SH_name(name);
  bool rv = selfobj->hasGroup(SH_name);
  return rv;
// splicer end class.DataGroup.method.has_group
}

bool SIDRE_datagroup_has_group_bufferify(SIDRE_datagroup * self,
                                         const char * name, int Lname)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.has_group_bufferify
  const std::string SH_name(name, Lname);
  bool rv = selfobj->hasGroup(SH_name);
  return rv;
// splicer end class.DataGroup.method.has_group_bufferify
}

SIDRE_datagroup * SIDRE_datagroup_get_group(SIDRE_datagroup * self,
                                            const char * name)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.get_group
  const std::string SH_name(name);
  DataGroup * rv = selfobj->getGroup(SH_name);
  return static_cast<SIDRE_datagroup *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.get_group
}

SIDRE_datagroup * SIDRE_datagroup_get_group_bufferify(SIDRE_datagroup * self,
                                                      const char * name,
                                                      int Lname)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.get_group_bufferify
  const std::string SH_name(name, Lname);
  DataGroup * rv = selfobj->getGroup(SH_name);
  return static_cast<SIDRE_datagroup *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.get_group_bufferify
}

SIDRE_IndexType SIDRE_datagroup_get_group_index(const SIDRE_datagroup * self,
                                                const char * name)
{
  const DataGroup * selfobj =
    static_cast<const DataGroup *>(static_cast<const void *>(self));
// splicer begin class.DataGroup.method.get_group_index
  const std::string SH_name(name);
  IndexType rv = selfobj->getGroupIndex(SH_name);
  return rv;
// splicer end class.DataGroup.method.get_group_index
}

SIDRE_IndexType SIDRE_datagroup_get_group_index_bufferify(
  const SIDRE_datagroup * self, const char * name, int Lname)
{
  const DataGroup * selfobj =
    static_cast<const DataGroup *>(static_cast<const void *>(self));
// splicer begin class.DataGroup.method.get_group_index_bufferify
  const std::string SH_name(name, Lname);
  IndexType rv = selfobj->getGroupIndex(SH_name);
  return rv;
// splicer end class.DataGroup.method.get_group_index_bufferify
}

const char * SIDRE_datagroup_get_group_name(const SIDRE_datagroup * self,
                                            SIDRE_IndexType idx)
{
  const DataGroup * selfobj =
    static_cast<const DataGroup *>(static_cast<const void *>(self));
// splicer begin class.DataGroup.method.get_group_name
  const std::string & rv = selfobj->getGroupName(idx);
// check for error
  if (!nameIsValid(rv))
  {
    return SIDRE_InvalidName;
  }

  return rv.c_str();
// splicer end class.DataGroup.method.get_group_name
}

void SIDRE_datagroup_get_group_name_bufferify(const SIDRE_datagroup * self,
                                              SIDRE_IndexType idx,
                                              char * SH_F_rv, int LSH_F_rv)
{
  const DataGroup * selfobj =
    static_cast<const DataGroup *>(static_cast<const void *>(self));
// splicer begin class.DataGroup.method.get_group_name_bufferify
  const std::string & rv = selfobj->getGroupName(idx);
// check for error
  if (!nameIsValid(rv))
  {
    std::memset(SH_F_rv, ' ', LSH_F_rv);
    return;
  }

  asctoolkit::shroud::FccCopy(SH_F_rv, LSH_F_rv, rv.c_str());
  return;
// splicer end class.DataGroup.method.get_group_name_bufferify
}

SIDRE_IndexType SIDRE_datagroup_get_first_valid_group_index(
  const SIDRE_datagroup * self)
{
  const DataGroup * selfobj =
    static_cast<const DataGroup *>(static_cast<const void *>(self));
// splicer begin class.DataGroup.method.get_first_valid_group_index
  IndexType rv = selfobj->getFirstValidGroupIndex();
  return rv;
// splicer end class.DataGroup.method.get_first_valid_group_index
}

SIDRE_IndexType SIDRE_datagroup_get_next_valid_group_index(
  const SIDRE_datagroup * self, SIDRE_IndexType idx)
{
  const DataGroup * selfobj =
    static_cast<const DataGroup *>(static_cast<const void *>(self));
// splicer begin class.DataGroup.method.get_next_valid_group_index
  IndexType rv = selfobj->getNextValidGroupIndex(idx);
  return rv;
// splicer end class.DataGroup.method.get_next_valid_group_index
}

SIDRE_datagroup * SIDRE_datagroup_create_group(SIDRE_datagroup * self,
                                               const char * name)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_group
  const std::string SH_name(name);
  DataGroup * rv = selfobj->createGroup(SH_name);
  return static_cast<SIDRE_datagroup *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_group
}

SIDRE_datagroup * SIDRE_datagroup_create_group_bufferify(SIDRE_datagroup * self,
                                                         const char * name,
                                                         int Lname)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_group_bufferify
  const std::string SH_name(name, Lname);
  DataGroup * rv = selfobj->createGroup(SH_name);
  return static_cast<SIDRE_datagroup *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_group_bufferify
}

void SIDRE_datagroup_destroy_group_name(SIDRE_datagroup * self,
                                        const char * name)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.destroy_group_name
  const std::string SH_name(name);
  selfobj->destroyGroup(SH_name);
  return;
// splicer end class.DataGroup.method.destroy_group_name
}

void SIDRE_datagroup_destroy_group_name_bufferify(SIDRE_datagroup * self,
                                                  const char * name, int Lname)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.destroy_group_name_bufferify
  const std::string SH_name(name, Lname);
  selfobj->destroyGroup(SH_name);
  return;
// splicer end class.DataGroup.method.destroy_group_name_bufferify
}

void SIDRE_datagroup_destroy_group_index(SIDRE_datagroup * self,
                                         SIDRE_IndexType idx)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.destroy_group_index
  selfobj->destroyGroup(idx);
  return;
// splicer end class.DataGroup.method.destroy_group_index
}

SIDRE_datagroup * SIDRE_datagroup_move_group(SIDRE_datagroup * self,
                                             SIDRE_datagroup * grp)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.move_group
  DataGroup * rv =
    selfobj->moveGroup(static_cast<DataGroup *>(static_cast<void *>(grp)));
  return static_cast<SIDRE_datagroup *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.move_group
}

void SIDRE_datagroup_print(const SIDRE_datagroup * self)
{
  const DataGroup * selfobj =
    static_cast<const DataGroup *>(static_cast<const void *>(self));
// splicer begin class.DataGroup.method.print
  selfobj->print();
  return;
// splicer end class.DataGroup.method.print
}

bool SIDRE_datagroup_is_equivalent_to(const SIDRE_datagroup * self,
                                      const SIDRE_datagroup * other)
{
  const DataGroup * selfobj =
    static_cast<const DataGroup *>(static_cast<const void *>(self));
// splicer begin class.DataGroup.method.is_equivalent_to
  bool rv =
    selfobj->isEquivalentTo(static_cast<const DataGroup *>(static_cast<const
                                                                       void *>(
                                                             other)));
  return rv;
// splicer end class.DataGroup.method.is_equivalent_to
}

void SIDRE_datagroup_save(const SIDRE_datagroup * self, const char * file_path,
                          const char * protocol)
{
  const DataGroup * selfobj =
    static_cast<const DataGroup *>(static_cast<const void *>(self));
// splicer begin class.DataGroup.method.save
  const std::string SH_file_path(file_path);
  const std::string SH_protocol(protocol);
  selfobj->save(SH_file_path, SH_protocol);
  return;
// splicer end class.DataGroup.method.save
}

void SIDRE_datagroup_save_bufferify(const SIDRE_datagroup * self,
                                    const char * file_path, int Lfile_path,
                                    const char * protocol, int Lprotocol)
{
  const DataGroup * selfobj =
    static_cast<const DataGroup *>(static_cast<const void *>(self));
// splicer begin class.DataGroup.method.save_bufferify
  const std::string SH_file_path(file_path, Lfile_path);
  const std::string SH_protocol(protocol, Lprotocol);
  selfobj->save(SH_file_path, SH_protocol);
  return;
// splicer end class.DataGroup.method.save_bufferify
}

void SIDRE_datagroup_load(SIDRE_datagroup * self, const char * file_path,
                          const char * protocol)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.load
  const std::string SH_file_path(file_path);
  const std::string SH_protocol(protocol);
  selfobj->load(SH_file_path, SH_protocol);
  return;
// splicer end class.DataGroup.method.load
}

void SIDRE_datagroup_load_bufferify(SIDRE_datagroup * self,
                                    const char * file_path, int Lfile_path,
                                    const char * protocol, int Lprotocol)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.load_bufferify
  const std::string SH_file_path(file_path, Lfile_path);
  const std::string SH_protocol(protocol, Lprotocol);
  selfobj->load(SH_file_path, SH_protocol);
  return;
// splicer end class.DataGroup.method.load_bufferify
}

void SIDRE_datagroup_load_external_data(SIDRE_datagroup * self,
                                        const char * file_path,
                                        const char * protocol)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.load_external_data
  const std::string SH_file_path(file_path);
  const std::string SH_protocol(protocol);
  selfobj->loadExternalData(SH_file_path, SH_protocol);
  return;
// splicer end class.DataGroup.method.load_external_data
}

void SIDRE_datagroup_load_external_data_bufferify(SIDRE_datagroup * self,
                                                  const char * file_path,
                                                  int Lfile_path,
                                                  const char * protocol,
                                                  int Lprotocol)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.load_external_data_bufferify
  const std::string SH_file_path(file_path, Lfile_path);
  const std::string SH_protocol(protocol, Lprotocol);
  selfobj->loadExternalData(SH_file_path, SH_protocol);
  return;
// splicer end class.DataGroup.method.load_external_data_bufferify
}

// splicer begin class.DataGroup.additional_functions
// splicer end class.DataGroup.additional_functions

}  // namespace asctoolkit
}  // namespace sidre
}  // extern "C"

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

void SIDRE_datagroup_get_name_bufferify(SIDRE_datagroup * self, char * name,
                                        int Lname)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.get_name_bufferify
  const std::string & rv = selfobj->getName();
  asctoolkit::shroud::FccCopy(name, Lname, rv.c_str());
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

size_t SIDRE_datagroup_get_num_views(SIDRE_datagroup * self)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.get_num_views
  size_t rv = selfobj->getNumViews();
  return rv;
// splicer end class.DataGroup.method.get_num_views
}

size_t SIDRE_datagroup_get_num_groups(SIDRE_datagroup * self)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.get_num_groups
  size_t rv = selfobj->getNumGroups();
  return rv;
// splicer end class.DataGroup.method.get_num_groups
}

bool SIDRE_datagroup_has_view(SIDRE_datagroup * self, const char * name)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.has_view
  std::string SH_name(name);
  bool rv = selfobj->hasView(SH_name);
  return rv;
// splicer end class.DataGroup.method.has_view
}

bool SIDRE_datagroup_has_view_bufferify(SIDRE_datagroup * self,
                                        const char * name, int Lname)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.has_view_bufferify
  std::string SH_name(name, Lname);
  bool rv = selfobj->hasView(SH_name);
  return rv;
// splicer end class.DataGroup.method.has_view_bufferify
}

SIDRE_dataview * SIDRE_datagroup_get_view_from_name(SIDRE_datagroup * self,
                                                    const char * name)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.get_view_from_name
  std::string SH_name(name);
  DataView * rv = selfobj->getView(SH_name);
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.get_view_from_name
}

SIDRE_dataview * SIDRE_datagroup_get_view_from_name_bufferify(
  SIDRE_datagroup * self, const char * name, int Lname)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.get_view_from_name_bufferify
  std::string SH_name(name, Lname);
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

SIDRE_IndexType SIDRE_datagroup_get_view_index(SIDRE_datagroup * self,
                                               const char * name)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.get_view_index
  std::string SH_name(name);
  IndexType rv = selfobj->getViewIndex(SH_name);
  return rv;
// splicer end class.DataGroup.method.get_view_index
}

SIDRE_IndexType SIDRE_datagroup_get_view_index_bufferify(SIDRE_datagroup * self,
                                                         const char * name,
                                                         int Lname)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.get_view_index_bufferify
  std::string SH_name(name, Lname);
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

void SIDRE_datagroup_get_view_name_bufferify(SIDRE_datagroup * self,
                                             SIDRE_IndexType idx, char * name,
                                             int Lname)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.get_view_name_bufferify
  const std::string & rv = selfobj->getViewName(idx);
// check for error
  if (!nameIsValid(rv))
  {
    std::memset(name, ' ', Lname);
    return;
  }

  asctoolkit::shroud::FccCopy(name, Lname, rv.c_str());
  return;
// splicer end class.DataGroup.method.get_view_name_bufferify
}

SIDRE_IndexType SIDRE_datagroup_get_first_valid_view_index(
  SIDRE_datagroup * self)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.get_first_valid_view_index
  IndexType rv = selfobj->getFirstValidViewIndex();
  return rv;
// splicer end class.DataGroup.method.get_first_valid_view_index
}

SIDRE_IndexType SIDRE_datagroup_get_next_valid_view_index(
  SIDRE_datagroup * self, SIDRE_IndexType idx)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.get_next_valid_view_index
  IndexType rv = selfobj->getNextValidViewIndex(idx);
  return rv;
// splicer end class.DataGroup.method.get_next_valid_view_index
}

SIDRE_dataview * SIDRE_datagroup_create_view_and_allocate(
  SIDRE_datagroup * self, const char * name, int type,
  SIDRE_SidreLength num_elems)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_and_allocate
  std::string SH_name(name);
  DataView * rv = selfobj->createViewAndAllocate(SH_name, getTypeID(
                                                   type), num_elems);
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_and_allocate
}

SIDRE_dataview * SIDRE_datagroup_create_view_and_allocate_bufferify(
  SIDRE_datagroup * self, const char * name, int Lname, int type,
  SIDRE_SidreLength num_elems)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_and_allocate_bufferify
  std::string SH_name(name, Lname);
  DataView * rv = selfobj->createViewAndAllocate(SH_name, getTypeID(
                                                   type), num_elems);
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_and_allocate_bufferify
}

SIDRE_dataview * SIDRE_datagroup_create_view_empty(SIDRE_datagroup * self,
                                                   const char * name)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_empty
  std::string SH_name(name);
  DataView * rv = selfobj->createView(SH_name);
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_empty
}

SIDRE_dataview * SIDRE_datagroup_create_view_empty_bufferify(
  SIDRE_datagroup * self, const char * name, int Lname)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_empty_bufferify
  std::string SH_name(name, Lname);
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
  std::string SH_name(name);
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
  std::string SH_name(name, Lname);
  DataView * rv = selfobj->createView(SH_name, getTypeID(type), num_elems);
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_from_type_bufferify
}

SIDRE_dataview * SIDRE_datagroup_create_view_into_buffer(SIDRE_datagroup * self,
                                                         const char * name,
                                                         SIDRE_databuffer * buff)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_into_buffer
  std::string SH_name(name);
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
  std::string SH_name(name, Lname);
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
  std::string SH_name(name);
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
  std::string SH_name(name, Lname);
  DataView * rv = selfobj->createView(SH_name, external_ptr);
  return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_external_bufferify
}

void SIDRE_datagroup_destroy_view(SIDRE_datagroup * self, const char * name)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.destroy_view
  std::string SH_name(name);
  selfobj->destroyView(SH_name);
  return;
// splicer end class.DataGroup.method.destroy_view
}

void SIDRE_datagroup_destroy_view_bufferify(SIDRE_datagroup * self,
                                            const char * name, int Lname)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.destroy_view_bufferify
  std::string SH_name(name, Lname);
  selfobj->destroyView(SH_name);
  return;
// splicer end class.DataGroup.method.destroy_view_bufferify
}

void SIDRE_datagroup_destroy_view_and_data_name(SIDRE_datagroup * self,
                                                const char * name)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.destroy_view_and_data_name
  std::string SH_name(name);
  selfobj->destroyViewAndData(SH_name);
  return;
// splicer end class.DataGroup.method.destroy_view_and_data_name
}

void SIDRE_datagroup_destroy_view_and_data_name_bufferify(
  SIDRE_datagroup * self, const char * name, int Lname)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.destroy_view_and_data_name_bufferify
  std::string SH_name(name, Lname);
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
  std::string SH_name(name);
  bool rv = selfobj->hasGroup(SH_name);
  return rv;
// splicer end class.DataGroup.method.has_group
}

bool SIDRE_datagroup_has_group_bufferify(SIDRE_datagroup * self,
                                         const char * name, int Lname)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.has_group_bufferify
  std::string SH_name(name, Lname);
  bool rv = selfobj->hasGroup(SH_name);
  return rv;
// splicer end class.DataGroup.method.has_group_bufferify
}

SIDRE_datagroup * SIDRE_datagroup_get_group(SIDRE_datagroup * self,
                                            const char * name)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.get_group
  std::string SH_name(name);
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
  std::string SH_name(name, Lname);
  DataGroup * rv = selfobj->getGroup(SH_name);
  return static_cast<SIDRE_datagroup *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.get_group_bufferify
}

SIDRE_IndexType SIDRE_datagroup_get_group_index(SIDRE_datagroup * self,
                                                const char * name)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.get_group_index
  std::string SH_name(name);
  IndexType rv = selfobj->getGroupIndex(SH_name);
  return rv;
// splicer end class.DataGroup.method.get_group_index
}

SIDRE_IndexType SIDRE_datagroup_get_group_index_bufferify(
  SIDRE_datagroup * self, const char * name, int Lname)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.get_group_index_bufferify
  std::string SH_name(name, Lname);
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

void SIDRE_datagroup_get_group_name_bufferify(SIDRE_datagroup * self,
                                              SIDRE_IndexType idx, char * name,
                                              int Lname)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.get_group_name_bufferify
  const std::string & rv = selfobj->getGroupName(idx);
// check for error
  if (!nameIsValid(rv))
  {
    std::memset(name, ' ', Lname);
    return;
  }

  asctoolkit::shroud::FccCopy(name, Lname, rv.c_str());
  return;
// splicer end class.DataGroup.method.get_group_name_bufferify
}

SIDRE_IndexType SIDRE_datagroup_get_first_valid_group_index(
  SIDRE_datagroup * self)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.get_first_valid_group_index
  IndexType rv = selfobj->getFirstValidGroupIndex();
  return rv;
// splicer end class.DataGroup.method.get_first_valid_group_index
}

SIDRE_IndexType SIDRE_datagroup_get_next_valid_group_index(
  SIDRE_datagroup * self, SIDRE_IndexType idx)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
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
  std::string SH_name(name);
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
  std::string SH_name(name, Lname);
  DataGroup * rv = selfobj->createGroup(SH_name);
  return static_cast<SIDRE_datagroup *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_group_bufferify
}

void SIDRE_datagroup_destroy_group_name(SIDRE_datagroup * self,
                                        const char * name)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.destroy_group_name
  std::string SH_name(name);
  selfobj->destroyGroup(SH_name);
  return;
// splicer end class.DataGroup.method.destroy_group_name
}

void SIDRE_datagroup_destroy_group_name_bufferify(SIDRE_datagroup * self,
                                                  const char * name, int Lname)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.destroy_group_name_bufferify
  std::string SH_name(name, Lname);
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

void SIDRE_datagroup_print(SIDRE_datagroup * self)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.print
  selfobj->print();
  return;
// splicer end class.DataGroup.method.print
}

void SIDRE_datagroup_save(SIDRE_datagroup * self, const char * obase,
                          const char * protocol)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.save
  std::string SH_obase(obase);
  std::string SH_protocol(protocol);
  selfobj->save(SH_obase, SH_protocol);
  return;
// splicer end class.DataGroup.method.save
}

void SIDRE_datagroup_save_bufferify(SIDRE_datagroup * self, const char * obase,
                                    int Lobase, const char * protocol,
                                    int Lprotocol)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.save_bufferify
  std::string SH_obase(obase, Lobase);
  std::string SH_protocol(protocol, Lprotocol);
  selfobj->save(SH_obase, SH_protocol);
  return;
// splicer end class.DataGroup.method.save_bufferify
}

void SIDRE_datagroup_load(SIDRE_datagroup * self, const char * obase,
                          const char * protocol)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.load
  std::string SH_obase(obase);
  std::string SH_protocol(protocol);
  selfobj->load(SH_obase, SH_protocol);
  return;
// splicer end class.DataGroup.method.load
}

void SIDRE_datagroup_load_bufferify(SIDRE_datagroup * self, const char * obase,
                                    int Lobase, const char * protocol,
                                    int Lprotocol)
{
  DataGroup * selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.load_bufferify
  std::string SH_obase(obase, Lobase);
  std::string SH_protocol(protocol, Lprotocol);
  selfobj->load(SH_obase, SH_protocol);
  return;
// splicer end class.DataGroup.method.load_bufferify
}

// splicer begin class.DataGroup.additional_functions
// splicer end class.DataGroup.additional_functions

}  // namespace asctoolkit
}  // namespace sidre
}  // extern "C"

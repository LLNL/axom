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
#include "sidre/DataGroup.hpp"
#include "sidre/SidreTypes.hpp"

extern "C" {
namespace asctoolkit {
namespace sidre {

const char * SIDRE_datagroup_get_name(const SIDRE_datagroup * self)
{
const DataGroup *selfobj = static_cast<const DataGroup *>(static_cast<const void *>(self));
// splicer begin class.DataGroup.method.get_name
const std::string & rv = selfobj->getName();
return rv.c_str();
// splicer end class.DataGroup.method.get_name
}

const SIDRE_datagroup * SIDRE_datagroup_get_parent(const SIDRE_datagroup * self)
{
const DataGroup *selfobj = static_cast<const DataGroup *>(static_cast<const void *>(self));
// splicer begin class.DataGroup.method.get_parent
const DataGroup * rv = selfobj->getParent();
return static_cast<const SIDRE_datagroup *>(static_cast<const void *>(rv));
// splicer end class.DataGroup.method.get_parent
}

const SIDRE_datastore * SIDRE_datagroup_get_data_store(const SIDRE_datagroup * self)
{
const DataGroup *selfobj = static_cast<const DataGroup *>(static_cast<const void *>(self));
// splicer begin class.DataGroup.method.get_data_store
const DataStore * rv = selfobj->getDataStore();
return static_cast<const SIDRE_datastore *>(static_cast<const void *>(rv));
// splicer end class.DataGroup.method.get_data_store
}

size_t SIDRE_datagroup_get_num_views(SIDRE_datagroup * self)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.get_num_views
size_t rv = selfobj->getNumViews();
return rv;
// splicer end class.DataGroup.method.get_num_views
}

size_t SIDRE_datagroup_get_num_groups(SIDRE_datagroup * self)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.get_num_groups
size_t rv = selfobj->getNumGroups();
return rv;
// splicer end class.DataGroup.method.get_num_groups
}

bool SIDRE_datagroup_has_view(SIDRE_datagroup * self, const char * name)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.has_view
bool rv = selfobj->hasView(name);
return rv;
// splicer end class.DataGroup.method.has_view
}

bool SIDRE_datagroup_has_view_bufferify(SIDRE_datagroup * self, const char * name, int Lname)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.has_view_bufferify
bool rv = selfobj->hasView(std::string(name, Lname));
return rv;
// splicer end class.DataGroup.method.has_view_bufferify
}

SIDRE_dataview * SIDRE_datagroup_get_view_from_name(SIDRE_datagroup * self, const char * name)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.get_view_from_name
DataView * rv = selfobj->getView(name);
return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.get_view_from_name
}

SIDRE_dataview * SIDRE_datagroup_get_view_from_name_bufferify(SIDRE_datagroup * self, const char * name, int Lname)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.get_view_from_name_bufferify
DataView * rv = selfobj->getView(std::string(name, Lname));
return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.get_view_from_name_bufferify
}

SIDRE_dataview * SIDRE_datagroup_get_view_from_index(SIDRE_datagroup * self, const SIDRE_IndexType idx)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.get_view_from_index
DataView * rv = selfobj->getView(idx);
return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.get_view_from_index
}

SIDRE_IndexType SIDRE_datagroup_get_view_index(SIDRE_datagroup * self, const char * name)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.get_view_index
IndexType rv = selfobj->getViewIndex(name);
return rv;
// splicer end class.DataGroup.method.get_view_index
}

SIDRE_IndexType SIDRE_datagroup_get_view_index_bufferify(SIDRE_datagroup * self, const char * name, int Lname)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.get_view_index_bufferify
IndexType rv = selfobj->getViewIndex(std::string(name, Lname));
return rv;
// splicer end class.DataGroup.method.get_view_index_bufferify
}

const char * SIDRE_datagroup_get_view_name(const SIDRE_datagroup * self, SIDRE_IndexType idx)
{
const DataGroup *selfobj = static_cast<const DataGroup *>(static_cast<const void *>(self));
// splicer begin class.DataGroup.method.get_view_name
const std::string & rv = selfobj->getViewName(idx);
if (! nameIsValid(rv)) {
    return SIDRE_InvalidName;
}

return rv.c_str();
// splicer end class.DataGroup.method.get_view_name
}

SIDRE_dataview * SIDRE_datagroup_create_view_and_allocate_from_type(SIDRE_datagroup * self, const char * name, int type, SIDRE_SidreLength num_elems)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_and_allocate_from_type
DataView * rv = selfobj->createViewAndAllocate(name, getTypeID(type), num_elems);
return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_and_allocate_from_type
}

SIDRE_dataview * SIDRE_datagroup_create_view_and_allocate_from_type_bufferify(SIDRE_datagroup * self, const char * name, int Lname, int type, SIDRE_SidreLength num_elems)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_and_allocate_from_type_bufferify
DataView * rv = selfobj->createViewAndAllocate(std::string(name, Lname), getTypeID(type), num_elems);
return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_and_allocate_from_type_bufferify
}

SIDRE_dataview * SIDRE_datagroup_create_view_empty(SIDRE_datagroup * self, const char * name)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_empty
DataView * rv = selfobj->createView(name);
return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_empty
}

SIDRE_dataview * SIDRE_datagroup_create_view_empty_bufferify(SIDRE_datagroup * self, const char * name, int Lname)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_empty_bufferify
DataView * rv = selfobj->createView(std::string(name, Lname));
return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_empty_bufferify
}

SIDRE_dataview * SIDRE_datagroup_create_view_from_type(SIDRE_datagroup * self, const char * name, int type, SIDRE_SidreLength num_elems)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_from_type
DataView * rv = selfobj->createView(name, getTypeID(type), num_elems);
return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_from_type
}

SIDRE_dataview * SIDRE_datagroup_create_view_from_type_bufferify(SIDRE_datagroup * self, const char * name, int Lname, int type, SIDRE_SidreLength num_elems)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_from_type_bufferify
DataView * rv = selfobj->createView(std::string(name, Lname), getTypeID(type), num_elems);
return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_from_type_bufferify
}

SIDRE_dataview * SIDRE_datagroup_create_view_into_buffer(SIDRE_datagroup * self, const char * name, SIDRE_databuffer * buff)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_into_buffer
DataView * rv = selfobj->createView(name, static_cast<DataBuffer *>(static_cast<void *>(buff)));
return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_into_buffer
}

SIDRE_dataview * SIDRE_datagroup_create_view_into_buffer_bufferify(SIDRE_datagroup * self, const char * name, int Lname, SIDRE_databuffer * buff)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_into_buffer_bufferify
DataView * rv = selfobj->createView(std::string(name, Lname), static_cast<DataBuffer *>(static_cast<void *>(buff)));
return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_into_buffer_bufferify
}

SIDRE_dataview * SIDRE_datagroup_create_view_external(SIDRE_datagroup * self, const char * name, void * external_ptr)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_external
DataView * rv = selfobj->createView(name, external_ptr);
return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_external
}

SIDRE_dataview * SIDRE_datagroup_create_view_external_bufferify(SIDRE_datagroup * self, const char * name, int Lname, void * external_ptr)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_view_external_bufferify
DataView * rv = selfobj->createView(std::string(name, Lname), external_ptr);
return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_view_external_bufferify
}

SIDRE_dataview * SIDRE_datagroup_create_external_view(SIDRE_datagroup * self, const char * name, void * external_data, int type, SIDRE_SidreLength num_elems)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_external_view
DataView * rv = selfobj->createExternalView(name, external_data, getTypeID(type), num_elems);
return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_external_view
}

SIDRE_dataview * SIDRE_datagroup_create_external_view_bufferify(SIDRE_datagroup * self, const char * name, int Lname, void * external_data, int type, SIDRE_SidreLength num_elems)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_external_view_bufferify
DataView * rv = selfobj->createExternalView(std::string(name, Lname), external_data, getTypeID(type), num_elems);
return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_external_view_bufferify
}

SIDRE_dataview * SIDRE_datagroup_create_external_view_with_shape(SIDRE_datagroup * self, const char * name, void * external_data, int type, int ndims, SIDRE_SidreLength * shape)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_external_view_with_shape
DataView * rv = selfobj->createExternalView(name, external_data, getTypeID(type), ndims, shape);
return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_external_view_with_shape
}

SIDRE_dataview * SIDRE_datagroup_create_external_view_with_shape_bufferify(SIDRE_datagroup * self, const char * name, int Lname, void * external_data, int type, int ndims, SIDRE_SidreLength * shape)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_external_view_with_shape_bufferify
DataView * rv = selfobj->createExternalView(std::string(name, Lname), external_data, getTypeID(type), ndims, shape);
return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_external_view_with_shape_bufferify
}

void SIDRE_datagroup_destroy_view(SIDRE_datagroup * self, const char * name)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.destroy_view
selfobj->destroyView(name);
return;
// splicer end class.DataGroup.method.destroy_view
}

void SIDRE_datagroup_destroy_view_bufferify(SIDRE_datagroup * self, const char * name, int Lname)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.destroy_view_bufferify
selfobj->destroyView(std::string(name, Lname));
return;
// splicer end class.DataGroup.method.destroy_view_bufferify
}

void SIDRE_datagroup_destroy_view_and_data(SIDRE_datagroup * self, const char * name)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.destroy_view_and_data
selfobj->destroyViewAndData(name);
return;
// splicer end class.DataGroup.method.destroy_view_and_data
}

void SIDRE_datagroup_destroy_view_and_data_bufferify(SIDRE_datagroup * self, const char * name, int Lname)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.destroy_view_and_data_bufferify
selfobj->destroyViewAndData(std::string(name, Lname));
return;
// splicer end class.DataGroup.method.destroy_view_and_data_bufferify
}

SIDRE_dataview * SIDRE_datagroup_move_view(SIDRE_datagroup * self, SIDRE_dataview * view)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.move_view
DataView * rv = selfobj->moveView(static_cast<DataView *>(static_cast<void *>(view)));
return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.move_view
}

SIDRE_dataview * SIDRE_datagroup_copy_view(SIDRE_datagroup * self, SIDRE_dataview * view)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.copy_view
DataView * rv = selfobj->copyView(static_cast<DataView *>(static_cast<void *>(view)));
return static_cast<SIDRE_dataview *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.copy_view
}

bool SIDRE_datagroup_has_group(SIDRE_datagroup * self, const char * name)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.has_group
bool rv = selfobj->hasGroup(name);
return rv;
// splicer end class.DataGroup.method.has_group
}

bool SIDRE_datagroup_has_group_bufferify(SIDRE_datagroup * self, const char * name, int Lname)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.has_group_bufferify
bool rv = selfobj->hasGroup(std::string(name, Lname));
return rv;
// splicer end class.DataGroup.method.has_group_bufferify
}

SIDRE_datagroup * SIDRE_datagroup_get_group(SIDRE_datagroup * self, const char * name)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.get_group
DataGroup * rv = selfobj->getGroup(name);
return static_cast<SIDRE_datagroup *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.get_group
}

SIDRE_datagroup * SIDRE_datagroup_get_group_bufferify(SIDRE_datagroup * self, const char * name, int Lname)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.get_group_bufferify
DataGroup * rv = selfobj->getGroup(std::string(name, Lname));
return static_cast<SIDRE_datagroup *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.get_group_bufferify
}

SIDRE_IndexType SIDRE_datagroup_get_group_index(SIDRE_datagroup * self, const char * name)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.get_group_index
IndexType rv = selfobj->getGroupIndex(name);
return rv;
// splicer end class.DataGroup.method.get_group_index
}

SIDRE_IndexType SIDRE_datagroup_get_group_index_bufferify(SIDRE_datagroup * self, const char * name, int Lname)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.get_group_index_bufferify
IndexType rv = selfobj->getGroupIndex(std::string(name, Lname));
return rv;
// splicer end class.DataGroup.method.get_group_index_bufferify
}

const char * SIDRE_datagroup_get_group_name(const SIDRE_datagroup * self, SIDRE_IndexType idx)
{
const DataGroup *selfobj = static_cast<const DataGroup *>(static_cast<const void *>(self));
// splicer begin class.DataGroup.method.get_group_name
const std::string & rv = selfobj->getGroupName(idx);
if (! nameIsValid(rv)) {
    return SIDRE_InvalidName;
}

return rv.c_str();
// splicer end class.DataGroup.method.get_group_name
}

SIDRE_datagroup * SIDRE_datagroup_create_group(SIDRE_datagroup * self, const char * name)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_group
DataGroup * rv = selfobj->createGroup(name);
return static_cast<SIDRE_datagroup *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_group
}

SIDRE_datagroup * SIDRE_datagroup_create_group_bufferify(SIDRE_datagroup * self, const char * name, int Lname)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.create_group_bufferify
DataGroup * rv = selfobj->createGroup(std::string(name, Lname));
return static_cast<SIDRE_datagroup *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.create_group_bufferify
}

void SIDRE_datagroup_destroy_group(SIDRE_datagroup * self, const char * name)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.destroy_group
selfobj->destroyGroup(name);
return;
// splicer end class.DataGroup.method.destroy_group
}

void SIDRE_datagroup_destroy_group_bufferify(SIDRE_datagroup * self, const char * name, int Lname)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.destroy_group_bufferify
selfobj->destroyGroup(std::string(name, Lname));
return;
// splicer end class.DataGroup.method.destroy_group_bufferify
}

SIDRE_datagroup * SIDRE_datagroup_move_group(SIDRE_datagroup * self, SIDRE_datagroup * grp)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.move_group
DataGroup * rv = selfobj->moveGroup(static_cast<DataGroup *>(static_cast<void *>(grp)));
return static_cast<SIDRE_datagroup *>(static_cast<void *>(rv));
// splicer end class.DataGroup.method.move_group
}

void SIDRE_datagroup_print(SIDRE_datagroup * self)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.print
selfobj->print();
return;
// splicer end class.DataGroup.method.print
}

void SIDRE_datagroup_save(SIDRE_datagroup * self, const char * obase, const char * protocol)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.save
selfobj->save(obase, protocol);
return;
// splicer end class.DataGroup.method.save
}

void SIDRE_datagroup_save_bufferify(SIDRE_datagroup * self, const char * obase, int Lobase, const char * protocol, int Lprotocol)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.save_bufferify
selfobj->save(std::string(obase, Lobase), std::string(protocol, Lprotocol));
return;
// splicer end class.DataGroup.method.save_bufferify
}

void SIDRE_datagroup_load(SIDRE_datagroup * self, const char * obase, const char * protocol)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.load
selfobj->load(obase, protocol);
return;
// splicer end class.DataGroup.method.load
}

void SIDRE_datagroup_load_bufferify(SIDRE_datagroup * self, const char * obase, int Lobase, const char * protocol, int Lprotocol)
{
DataGroup *selfobj = static_cast<DataGroup *>(static_cast<void *>(self));
// splicer begin class.DataGroup.method.load_bufferify
selfobj->load(std::string(obase, Lobase), std::string(protocol, Lprotocol));
return;
// splicer end class.DataGroup.method.load_bufferify
}

// splicer begin class.DataGroup.additional_functions
// splicer end class.DataGroup.additional_functions

}  // namespace asctoolkit
}  // namespace sidre
}  // extern "C"

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
#define EXAMPLE_WRAPPER_IMPL
#include "wrapDataGroup.h"
#include "sidre/DataGroup.hpp"
#include "sidre/SidreWrapperHelpers.hpp"

extern "C" {
namespace asctoolkit {
namespace sidre {

const char * ATK_datagroup_get_name(const ATK_datagroup * self)
{
const DataGroup *selfobj = static_cast<const DataGroup *>(self);
// splicer begin class.DataGroup.method.getName
const std::string & rv = selfobj->getName();
return isNameValid(rv) ? rv.c_str() : ATK_InvalidName;
// splicer end class.DataGroup.method.getName
}

const ATK_datagroup * ATK_datagroup_get_parent(const ATK_datagroup * self)
{
const DataGroup *selfobj = static_cast<const DataGroup *>(self);
// splicer begin class.DataGroup.method.getParent
const DataGroup * rv = selfobj->getParent();
return rv;
// splicer end class.DataGroup.method.getParent
}

const ATK_datastore * ATK_datagroup_get_data_store(const ATK_datagroup * self)
{
const DataGroup *selfobj = static_cast<const DataGroup *>(self);
// splicer begin class.DataGroup.method.getDataStore
const DataStore * rv = selfobj->getDataStore();
return rv;
// splicer end class.DataGroup.method.getDataStore
}

size_t ATK_datagroup_get_num_views(ATK_datagroup * self)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin class.DataGroup.method.getNumViews
size_t rv = selfobj->getNumViews();
return rv;
// splicer end class.DataGroup.method.getNumViews
}

size_t ATK_datagroup_get_num_groups(ATK_datagroup * self)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin class.DataGroup.method.getNumGroups
size_t rv = selfobj->getNumGroups();
return rv;
// splicer end class.DataGroup.method.getNumGroups
}

bool ATK_datagroup_has_view(ATK_datagroup * self, const char * name)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin class.DataGroup.method.hasView
bool rv = selfobj->hasView(name);
return rv;
// splicer end class.DataGroup.method.hasView
}

ATK_dataview * ATK_datagroup_create_view_and_buffer_simple(ATK_datagroup * self, const char * name)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin class.DataGroup.method.createViewAndBuffer
DataView * rv = selfobj->createViewAndBuffer(name);
return rv;
// splicer end class.DataGroup.method.createViewAndBuffer
}

ATK_dataview * ATK_datagroup_create_view_and_buffer_from_type(ATK_datagroup * self, const char * name, int type, ATK_SidreLength len)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin class.DataGroup.method.createViewAndBuffer
DataView * rv = selfobj->createViewAndBuffer(name, getTypeID(type), len);
return rv;
// splicer end class.DataGroup.method.createViewAndBuffer
}

ATK_dataview * ATK_datagroup_create_opaque_view(ATK_datagroup * self, const char * name, void * opaque_ptr)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin class.DataGroup.method.createOpaqueView
DataView * rv = selfobj->createOpaqueView(name, opaque_ptr);
return rv;
// splicer end class.DataGroup.method.createOpaqueView
}

ATK_dataview * ATK_datagroup_create_view(ATK_datagroup * self, const char * name, ATK_databuffer * buff)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin class.DataGroup.method.createView
DataView * rv = selfobj->createView(name, static_cast<DataBuffer *>(buff));
return rv;
// splicer end class.DataGroup.method.createView
}

ATK_dataview * ATK_datagroup_create_external_view(ATK_datagroup * self, const char * name, void * external_data, int type, ATK_SidreLength len)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin class.DataGroup.method.createExternalView
DataView * rv = selfobj->createExternalView(name, external_data, getTypeID(type), len);
return rv;
// splicer end class.DataGroup.method.createExternalView
}

ATK_dataview * ATK_datagroup_move_view(ATK_datagroup * self, ATK_dataview * view)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin class.DataGroup.method.moveView
DataView * rv = selfobj->moveView(static_cast<DataView *>(view));
return rv;
// splicer end class.DataGroup.method.moveView
}

ATK_dataview * ATK_datagroup_copy_view(ATK_datagroup * self, ATK_dataview * view)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin class.DataGroup.method.copyView
DataView * rv = selfobj->copyView(static_cast<DataView *>(view));
return rv;
// splicer end class.DataGroup.method.copyView
}

void ATK_datagroup_destroy_view_and_buffer(ATK_datagroup * self, const char * name)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin class.DataGroup.method.destroyViewAndBuffer
selfobj->destroyViewAndBuffer(name);
return;
// splicer end class.DataGroup.method.destroyViewAndBuffer
}

ATK_dataview * ATK_datagroup_get_view(ATK_datagroup * self, const char * name)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin class.DataGroup.method.getView
DataView * rv = selfobj->getView(name);
return rv;
// splicer end class.DataGroup.method.getView
}

ATK_IndexType ATK_datagroup_get_view_index(ATK_datagroup * self, const char * name)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin class.DataGroup.method.getViewIndex
IndexType rv = selfobj->getViewIndex(name);
return rv;
// splicer end class.DataGroup.method.getViewIndex
}

const char * ATK_datagroup_get_view_name(const ATK_datagroup * self, ATK_IndexType idx)
{
const DataGroup *selfobj = static_cast<const DataGroup *>(self);
// splicer begin class.DataGroup.method.getViewName
const std::string & rv = selfobj->getViewName(idx);
return isNameValid(rv) ? rv.c_str() : ATK_InvalidName;
// splicer end class.DataGroup.method.getViewName
}

bool ATK_datagroup_has_group(ATK_datagroup * self, const char * name)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin class.DataGroup.method.hasGroup
bool rv = selfobj->hasGroup(name);
return rv;
// splicer end class.DataGroup.method.hasGroup
}

ATK_datagroup * ATK_datagroup_create_group(ATK_datagroup * self, const char * name)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin class.DataGroup.method.createGroup
DataGroup * rv = selfobj->createGroup(name);
return rv;
// splicer end class.DataGroup.method.createGroup
}

ATK_datagroup * ATK_datagroup_move_group(ATK_datagroup * self, ATK_datagroup * grp)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin class.DataGroup.method.moveGroup
DataGroup * rv = selfobj->moveGroup(static_cast<DataGroup *>(grp));
return rv;
// splicer end class.DataGroup.method.moveGroup
}

void ATK_datagroup_destroy_group(ATK_datagroup * self, const char * name)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin class.DataGroup.method.destroyGroup
selfobj->destroyGroup(name);
return;
// splicer end class.DataGroup.method.destroyGroup
}

ATK_datagroup * ATK_datagroup_get_group(ATK_datagroup * self, const char * name)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin class.DataGroup.method.getGroup
DataGroup * rv = selfobj->getGroup(name);
return rv;
// splicer end class.DataGroup.method.getGroup
}

ATK_IndexType ATK_datagroup_get_group_index(ATK_datagroup * self, const char * name)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin class.DataGroup.method.getGroupIndex
IndexType rv = selfobj->getGroupIndex(name);
return rv;
// splicer end class.DataGroup.method.getGroupIndex
}

const char * ATK_datagroup_get_group_name(const ATK_datagroup * self, ATK_IndexType idx)
{
const DataGroup *selfobj = static_cast<const DataGroup *>(self);
// splicer begin class.DataGroup.method.getGroupName
const std::string & rv = selfobj->getGroupName(idx);
return isNameValid(rv) ? rv.c_str() : ATK_InvalidName;
// splicer end class.DataGroup.method.getGroupName
}

int ATK_datagroup_get_group_name_length(ATK_datagroup * self, ATK_IndexType idx)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin class.DataGroup.method.GetGroupNameLength
return selfobj->getGroupName(idx).length();
// splicer end class.DataGroup.method.GetGroupNameLength
}

void ATK_datagroup_print(ATK_datagroup * self)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin class.DataGroup.method.print
selfobj->print();
return;
// splicer end class.DataGroup.method.print
}

void ATK_datagroup_save(ATK_datagroup * self, const char * obase, const char * protocol)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin class.DataGroup.method.save
selfobj->save(obase, protocol);
return;
// splicer end class.DataGroup.method.save
}

void ATK_datagroup_load(ATK_datagroup * self, const char * obase, const char * protocol)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin class.DataGroup.method.load
selfobj->load(obase, protocol);
return;
// splicer end class.DataGroup.method.load
}

// splicer begin class.DataGroup.additional_functions
// identical to ATK_datagroup_get_group_name except return a pointer to an empty string.
static const char * empty_string  = "";

const char * ATK_datagroup_get_group_name_with_error_check(const ATK_datagroup * self, ATK_IndexType idx)
{
    const DataGroup *selfobj = static_cast<const DataGroup *>(self);
    const std::string & rv = selfobj->getGroupName(idx);
    return isNameValid(rv) ? rv.c_str() : empty_string;
}
// splicer end class.DataGroup.additional_functions

}  // namespace asctoolkit
}  // namespace sidre
}  // extern "C"

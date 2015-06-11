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

extern "C" {
namespace asctoolkit {
namespace sidre {

const char * ATK_datagroup_get_name(const ATK_datagroup * self)
{
const DataGroup *selfobj = static_cast<const DataGroup *>(self);
// splicer begin
const std::string & rv = selfobj->getName();
return rv.c_str();
// splicer end
}

const ATK_datagroup * ATK_datagroup_get_parent(const ATK_datagroup * self)
{
const DataGroup *selfobj = static_cast<const DataGroup *>(self);
// splicer begin
const DataGroup * rv = selfobj->getParent();
return rv;
// splicer end
}

const ATK_datastore * ATK_datagroup_get_data_store(const ATK_datagroup * self)
{
const DataGroup *selfobj = static_cast<const DataGroup *>(self);
// splicer begin
const DataStore * rv = selfobj->getDataStore();
return rv;
// splicer end
}

bool ATK_datagroup_has_view(ATK_datagroup * self, const char * name)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin
bool rv = selfobj->hasView(name);
return rv;
// splicer end
}

ATK_dataview * ATK_datagroup_create_view_and_buffer(ATK_datagroup * self, const char * name)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin
DataView * rv = selfobj->createViewAndBuffer(name);
return rv;
// splicer end
}

ATK_dataview * ATK_datagroup_create_view_and_buffer_from_type(ATK_datagroup * self, const char * name, int type, long len)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin
DataView * rv = selfobj->createViewAndBuffer(name, getTypeID(type), len);
return rv;
// splicer end
}

ATK_dataview * ATK_datagroup_create_view(ATK_datagroup * self, const char * name, ATK_databuffer * buff)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin
DataView * rv = selfobj->createView(name, static_cast<DataBuffer *>(buff));
return rv;
// splicer end
}

ATK_dataview * ATK_datagroup_move_view(ATK_datagroup * self, ATK_dataview * view)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin
DataView * rv = selfobj->moveView(static_cast<DataView *>(view));
return rv;
// splicer end
}

ATK_dataview * ATK_datagroup_copy_view(ATK_datagroup * self, ATK_dataview * view)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin
DataView * rv = selfobj->copyView(static_cast<DataView *>(view));
return rv;
// splicer end
}

void ATK_datagroup_destroy_view_and_buffer(ATK_datagroup * self, const char * name)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin
selfobj->destroyViewAndBuffer(name);
return;
// splicer end
}

ATK_dataview * ATK_datagroup_get_view(ATK_datagroup * self, const char * name)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin
DataView * rv = selfobj->getView(name);
return rv;
// splicer end
}

ATK_IndexType ATK_datagroup_get_view_index(ATK_datagroup * self, const char * name)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin
IndexType rv = selfobj->getViewIndex(name);
return rv;
// splicer end
}

const char * ATK_datagroup_get_view_name(const ATK_datagroup * self, ATK_IndexType idx)
{
const DataGroup *selfobj = static_cast<const DataGroup *>(self);
// splicer begin
const std::string & rv = selfobj->getViewName(idx);
return rv.c_str();
// splicer end
}

size_t ATK_datagroup_get_num_views(ATK_datagroup * self)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin
size_t rv = selfobj->getNumViews();
return rv;
// splicer end
}

bool ATK_datagroup_has_group(ATK_datagroup * self, const char * name)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin
bool rv = selfobj->hasGroup(name);
return rv;
// splicer end
}

ATK_datagroup * ATK_datagroup_create_group(ATK_datagroup * self, const char * name)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin
DataGroup * rv = selfobj->createGroup(name);
return rv;
// splicer end
}

ATK_datagroup * ATK_datagroup_move_group(ATK_datagroup * self, ATK_datagroup * grp)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin
DataGroup * rv = selfobj->moveGroup(static_cast<DataGroup *>(grp));
return rv;
// splicer end
}

void ATK_datagroup_destroy_group(ATK_datagroup * self, const char * name)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin
selfobj->destroyGroup(name);
return;
// splicer end
}

ATK_datagroup * ATK_datagroup_get_group(ATK_datagroup * self, const char * name)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin
DataGroup * rv = selfobj->getGroup(name);
return rv;
// splicer end
}

ATK_IndexType ATK_datagroup_get_group_index(ATK_datagroup * self, const char * name)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin
IndexType rv = selfobj->getGroupIndex(name);
return rv;
// splicer end
}

const char * ATK_datagroup_get_group_name(const ATK_datagroup * self, ATK_IndexType idx)
{
const DataGroup *selfobj = static_cast<const DataGroup *>(self);
// splicer begin
const std::string & rv = selfobj->getGroupName(idx);
return rv.c_str();
// splicer end
}

size_t ATK_datagroup_get_num_groups(ATK_datagroup * self)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin
size_t rv = selfobj->getNumGroups();
return rv;
// splicer end
}

void ATK_datagroup_print(ATK_datagroup * self)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin
selfobj->print();
return;
// splicer end
}

void ATK_datagroup_save(ATK_datagroup * self, const char * obase, const char * protocol)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin
selfobj->save(obase, protocol);
return;
// splicer end
}

void ATK_datagroup_load(ATK_datagroup * self, const char * obase, const char * protocol)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin
selfobj->load(obase, protocol);
return;
// splicer end
}

}  // namespace asctoolkit
}  // namespace sidre
}  // extern "C"

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

ATK_datagroup * ATK_datagroup_get_parent(ATK_datagroup * self)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin
DataGroup * rv = selfobj->getParent();
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

ATK_datagroup * ATK_datagroup_create_group(ATK_datagroup * self, const char * name)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin
DataGroup * rv = selfobj->createGroup(name);
return rv;
// splicer end
}

}  // namespace asctoolkit
}  // namespace sidre
}  // extern "C"

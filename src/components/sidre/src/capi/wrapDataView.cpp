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
#define EXAMPLE_WRAPPER_IMPL
#include "wrapDataView.h"
#include "sidre/DataView.hpp"

extern "C" {
namespace asctoolkit {
namespace sidre {

ATK_dataview * ATK_dataview_declare(ATK_dataview * self, int type, long len)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin
DataView * rv = selfobj->declare(getTypeID(type), len);
return rv;
// splicer end
}

ATK_dataview * ATK_dataview_allocate(ATK_dataview * self, int type, long len)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin
DataView * rv = selfobj->allocate(getTypeID(type), len);
return rv;
// splicer end
}

bool ATK_dataview_has_buffer(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin
bool rv = selfobj->hasBuffer();
return rv;
// splicer end
}

const char * ATK_dataview_get_name(const ATK_dataview * self)
{
const DataView *selfobj = static_cast<const DataView *>(self);
// splicer begin
const std::string & rv = selfobj->getName();
return rv.c_str();
// splicer end
}

ATK_databuffer * ATK_dataview_get_buffer(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin
DataBuffer * rv = selfobj->getBuffer();
return rv;
// splicer end
}

ATK_datagroup * ATK_dataview_get_owning_group(ATK_dataview * self)
{
DataView *selfobj = static_cast<DataView *>(self);
// splicer begin
DataGroup * rv = selfobj->getOwningGroup();
return rv;
// splicer end
}

}  // namespace asctoolkit
}  // namespace sidre
}  // extern "C"

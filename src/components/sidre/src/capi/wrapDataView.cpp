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

const char * ATK_dataview_get_name(const ATK_dataview * self)
{
const DataView *selfobj = static_cast<const DataView *>(self);
// splicer begin
const std::string & rv = selfobj->getName();
return rv.c_str();
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

//
// Copyright (c) 2015, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory.
//
// All rights reserved.
//
// This source code cannot be distributed without permission and
// further review from Lawrence Livermore National Laboratory.
//
// wrapDataBuffer.cpp
#define EXAMPLE_WRAPPER_IMPL
#include "wrapDataBuffer.h"
#include "sidre/DataBuffer.hpp"

extern "C" {
namespace asctoolkit {
namespace sidre {

ATK_IndexType ATK_databuffer_get_index(ATK_databuffer * self)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(self);
// splicer begin
IndexType rv = selfobj->getIndex();
return rv;
// splicer end
}

ATK_databuffer * ATK_databuffer_declare(ATK_databuffer * self, int type, ATK_SidreLength len)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(self);
// splicer begin
DataBuffer * rv = selfobj->declare(getTypeID(type), len);
return rv;
// splicer end
}

ATK_databuffer * ATK_databuffer_allocate_noargs(ATK_databuffer * self)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(self);
// splicer begin
DataBuffer * rv = selfobj->allocate();
return rv;
// splicer end
}

ATK_databuffer * ATK_databuffer_allocate_from_type(ATK_databuffer * self, int type, ATK_SidreLength len)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(self);
// splicer begin
DataBuffer * rv = selfobj->allocate(getTypeID(type), len);
return rv;
// splicer end
}

ATK_databuffer * ATK_databuffer_reallocate(ATK_databuffer * self, int type, ATK_SidreLength len)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(self);
// splicer begin
DataBuffer * rv = selfobj->reallocate(getTypeID(type), len);
return rv;
// splicer end
}

void * ATK_databuffer_get_data(ATK_databuffer * self)
{
DataBuffer *selfobj = static_cast<DataBuffer *>(self);
// splicer begin
void * rv = selfobj->getData();
return rv;
// splicer end
}

}  // namespace asctoolkit
}  // namespace sidre
}  // extern "C"

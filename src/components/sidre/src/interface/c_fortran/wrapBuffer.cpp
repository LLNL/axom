// wrapBuffer.cpp
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
// wrapBuffer.cpp
#include "wrapBuffer.h"
#include "sidre/Buffer.hpp"
#include "sidre/SidreTypes.hpp"

extern "C" {
namespace axom {
namespace sidre {

SIDRE_IndexType SIDRE_buffer_get_index(const SIDRE_buffer * self)
{
const Buffer *selfobj = static_cast<const Buffer *>(static_cast<const void *>(self));
// splicer begin class.Buffer.method.get_index
    IndexType rv = selfobj->getIndex();
    return rv;
// splicer end class.Buffer.method.get_index
}

size_t SIDRE_buffer_get_num_views(const SIDRE_buffer * self)
{
const Buffer *selfobj = static_cast<const Buffer *>(static_cast<const void *>(self));
// splicer begin class.Buffer.method.get_num_views
    size_t rv = selfobj->getNumViews();
    return rv;
// splicer end class.Buffer.method.get_num_views
}

void SIDRE_buffer_describe(SIDRE_buffer * self, int type, SIDRE_SidreLength num_elems)
{
Buffer *selfobj = static_cast<Buffer *>(static_cast<void *>(self));
// splicer begin class.Buffer.method.describe
    selfobj->describe(getTypeID(type), num_elems);
    return;
// splicer end class.Buffer.method.describe
}

void SIDRE_buffer_allocate_existing(SIDRE_buffer * self)
{
Buffer *selfobj = static_cast<Buffer *>(static_cast<void *>(self));
// splicer begin class.Buffer.method.allocate_existing
    selfobj->allocate();
    return;
// splicer end class.Buffer.method.allocate_existing
}

void SIDRE_buffer_allocate_from_type(SIDRE_buffer * self, int type, SIDRE_SidreLength num_elems)
{
Buffer *selfobj = static_cast<Buffer *>(static_cast<void *>(self));
// splicer begin class.Buffer.method.allocate_from_type
    selfobj->allocate(getTypeID(type), num_elems);
    return;
// splicer end class.Buffer.method.allocate_from_type
}

void SIDRE_buffer_reallocate(SIDRE_buffer * self, SIDRE_SidreLength num_elems)
{
Buffer *selfobj = static_cast<Buffer *>(static_cast<void *>(self));
// splicer begin class.Buffer.method.reallocate
    selfobj->reallocate(num_elems);
    return;
// splicer end class.Buffer.method.reallocate
}

void * SIDRE_buffer_get_void_ptr(SIDRE_buffer * self)
{
Buffer *selfobj = static_cast<Buffer *>(static_cast<void *>(self));
// splicer begin class.Buffer.method.get_void_ptr
    void * rv = selfobj->getVoidPtr();
    return rv;
// splicer end class.Buffer.method.get_void_ptr
}

int SIDRE_buffer_get_type_id(const SIDRE_buffer * self)
{
const Buffer *selfobj = static_cast<const Buffer *>(static_cast<const void *>(self));
// splicer begin class.Buffer.method.get_type_id
    TypeID rv = selfobj->getTypeID();
    return static_cast<int>(rv);
// splicer end class.Buffer.method.get_type_id
}

size_t SIDRE_buffer_get_num_elements(const SIDRE_buffer * self)
{
const Buffer *selfobj = static_cast<const Buffer *>(static_cast<const void *>(self));
// splicer begin class.Buffer.method.get_num_elements
    size_t rv = selfobj->getNumElements();
    return rv;
// splicer end class.Buffer.method.get_num_elements
}

size_t SIDRE_buffer_get_total_bytes(const SIDRE_buffer * self)
{
const Buffer *selfobj = static_cast<const Buffer *>(static_cast<const void *>(self));
// splicer begin class.Buffer.method.get_total_bytes
    size_t rv = selfobj->getTotalBytes();
    return rv;
// splicer end class.Buffer.method.get_total_bytes
}

size_t SIDRE_buffer_get_bytes_per_element(const SIDRE_buffer * self)
{
const Buffer *selfobj = static_cast<const Buffer *>(static_cast<const void *>(self));
// splicer begin class.Buffer.method.get_bytes_per_element
    size_t rv = selfobj->getBytesPerElement();
    return rv;
// splicer end class.Buffer.method.get_bytes_per_element
}

void SIDRE_buffer_print(const SIDRE_buffer * self)
{
const Buffer *selfobj = static_cast<const Buffer *>(static_cast<const void *>(self));
// splicer begin class.Buffer.method.print
    selfobj->print();
    return;
// splicer end class.Buffer.method.print
}

// splicer begin class.Buffer.additional_functions
// splicer end class.Buffer.additional_functions

}  // namespace axom
}  // namespace sidre
}  // extern "C"

/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

/*!
 ******************************************************************************
 *
 * \file
 *
 * \brief   Implementation file for DataBuffer class.
 *
 ******************************************************************************
 */


// Associated header file
#include "DataBuffer.hpp"

// Standard C++ headers
#include <algorithm>

// Other CS Toolkit headers
#include "common/CommonTypes.hpp"
#include "common/Utilities.hpp"

// SiDRe project headers
#include "DataGroup.hpp"
#include "DataView.hpp"
#include "SidreTypes.hpp"


namespace asctoolkit
{
namespace sidre
{


/*
*************************************************************************
*
* Declare and allocate data described using a Conduit schema.
*
*************************************************************************
*/
DataBuffer* DataBuffer::declare(ATK_TypeID type, long len)
{
#if 0
    conduit::index_t dtype_id = type;
    conduit::index_t num_elements = len;
    const DataType* dtype = DataType(dtype_id, num_elements, 0, 0, 0, 0);
    //	const DataType& dtype = DataType(type, len, 0, 0, 0, 0, 0);
    m_schema.set(dtype);
#endif
    ATK_ASSERT_MSG(type == 0, "Bad Type");
    ATK_ASSERT_MSG(len < 0, "Bad Length");
    return this;
}

DataBuffer* DataBuffer::allocate(const Schema& schema)
{
    declare(schema);
    allocate();
    return this;
}

/*
*************************************************************************
*
* Declare and allocate data described using a Conduit pre-defined data type.
*
*************************************************************************
*/
DataBuffer* DataBuffer::allocate(const DataType& dtype)
{
    declare(dtype);
    allocate();
    return this;
}

/*
*************************************************************************
*
* Reallocate data using a Conduit schema.
*
*************************************************************************
*/
DataBuffer* DataBuffer::reallocate(const Schema& schema)
{
    //  make sure realloc actually makes sense
    ATK_ASSERT_MSG(m_data != ATK_NULLPTR,
                   "Attempting to reallocate an unallocated buffer");

    std::size_t realloc_size = schema.total_bytes();
    ATK_ASSERT_MSG(realloc_size > 0,
                   "Attempting to reallocate buffer to size 0");

    void* realloc_data = allocateBytes(realloc_size);

    // use conduit to get data from old to new schema.
    Node n;
    n.set_external(schema,realloc_data);
    // use conduit update, need more error checking.
    n.update(m_node);

    // cleanup old data
    cleanup();

    // set the buffer to use the new schema
    m_schema = schema;
    
    // let the buffer hold the new data
    m_data = realloc_data;
    
    // update the buffer's Conduit Node
    m_node.set_external(m_schema,m_data);
    return this;
}

/*
*************************************************************************
*
* Reallocate data using a basic Conduit data type.
*
*************************************************************************
*/
DataBuffer* DataBuffer::reallocate(const DataType& dtype)
{
    Schema s(dtype);
    reallocate(s);
    return this;
}

/*
*************************************************************************
*
* Allocate data previously declared.
*
*************************************************************************
*/
DataBuffer* DataBuffer::allocate()
{
    // cleanup old data
    cleanup();
    std::size_t alloc_size = m_schema.total_bytes();
    ATK_ASSERT_MSG(alloc_size > 0,
                   "Attempting to allocate buffer of size 0");
    m_data = allocateBytes(alloc_size);
    m_node.set_external(m_schema,m_data);
    return this;
}


/*
*************************************************************************
*
* Copy data buffer description to given Conduit node.
*
*************************************************************************
*/
void DataBuffer::info(Node &n) const
{
    n["uid"].set(m_uid);
    n["schema"].set(m_schema.to_json());
    n["node"].set(m_node.to_json());
}

/*
*************************************************************************
*   
* Print JSON description of data buffer to stdout.
*   
*************************************************************************
*/
void DataBuffer::print() const
{
    Node n;
    info(n);
    n.print();
}


/*
*************************************************************************
*   
* PRIVATE ctor taking unique id.
*   
*************************************************************************
*/
DataBuffer::DataBuffer( IDType uid ) :
    m_uid(uid),
    m_views(),
    m_data(ATK_NULLPTR),
    m_node(),
    m_schema()
{
    
}


/*
*************************************************************************
*   
* PRIVATE copy ctor.
*   
*************************************************************************
*/
DataBuffer::DataBuffer(const DataBuffer& source ) :
    m_uid(source.m_uid),
    m_views(source.m_views),
    m_data(source.m_data),
    m_node(source.m_node),
    m_schema(source.m_schema)
{
// disallow?
}


/*
*************************************************************************
*   
* PRIVATE dtor.
*   
*************************************************************************
*/
DataBuffer::~DataBuffer()
{
    cleanup();
}


/*
*************************************************************************
*   
* PRIVATE method to attach data view.
*   
*************************************************************************
*/
void DataBuffer::attachView( DataView* view )
{
    m_views.push_back( view );
}


/*
*************************************************************************
*   
* PRIVATE method to detach data view.
*   
*************************************************************************
*/
void DataBuffer::detachView( DataView* view )
{
    //Find new end iterator
    std::vector<DataView*>::iterator pos = std::remove(m_views.begin(),
                                                       m_views.end(),
                                                       view);
    // check if pos is ok?
    //Erase the "removed" elements.
    m_views.erase(pos, m_views.end());
}

/*
*************************************************************************
*   
* PRIVATE cleanup
*   
*************************************************************************
*/
void DataBuffer::cleanup()
{
    // cleanup alloced data
    if(m_data != ATK_NULLPTR)
    {
        releaseBytes(m_data);
        m_data = ATK_NULLPTR;
    }
}

/*
*************************************************************************
*   
* PRIVATE allocateBytes
*   
*************************************************************************
*/
void *DataBuffer::allocateBytes(std::size_t numBytes)
{
    ATK_ASSERT_MSG(numBytes > 0,
                   "Attempting to allocate 0 bytes");

    char *data = new char[numBytes];
    return ((void *)data);
}

/*
*************************************************************************
*   
* PRIVATE releaseBytes
*   
*************************************************************************
*/
void DataBuffer::releaseBytes(void *ptr)
{
    delete [] ((char*)ptr);
}



} /* end namespace sidre */
} /* end namespace asctoolkit */


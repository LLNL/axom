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
#include "common/Utilities.hpp"

// SiDRe project headers
#include "DataGroup.hpp"
#include "DataView.hpp"


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
* Allocate data previosly declared.
*
*************************************************************************
*/
DataBuffer* DataBuffer::allocate()
{
    std::size_t alloc_size = m_schema.total_bytes();

    ATK_ASSERT_MSG(alloc_size > 0, "Attempting to allocate buffer of size 0");   
    m_memblob.resize(alloc_size);
    m_data = m_memblob.data();

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
    n["descriptor"].set(m_schema.to_json());
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
DataBuffer::DataBuffer( common::IDType uid ) :
    m_uid(uid),
    m_views(),
    m_data(ATK_NULLPTR),
    m_memblob(),
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
    m_memblob(source.m_memblob),
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


} /* end namespace sidre */
} /* end namespace asctoolkit */


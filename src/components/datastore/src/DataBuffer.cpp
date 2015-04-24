
#include "DataBuffer.hpp"
#include "DataGroup.hpp"
#include "DataView.hpp"
#include <algorithm>

#include "Utilities.hpp"

namespace DataStoreNS
{


DataBuffer::DataBuffer( const IDType uid ) :
    m_uid(uid),
    m_views(),
    m_data(nullptr),
    m_memblob(),
    m_node(),
    m_schema()
{
    
}


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


DataBuffer::~DataBuffer()
{

}


/// init calls set declare, allocate
DataBuffer* DataBuffer::Allocate(const Schema &schema)
{
    Declare(schema);
    Allocate();
    return this;
}

/// init calls set declare, allocate
DataBuffer* DataBuffer::Allocate(const DataType &dtype)
{
    Declare(dtype);
    Allocate();
    return this;
}


DataBuffer* DataBuffer::Allocate()
{
    std::size_t alloc_size = m_schema.total_bytes();

    ATK_ASSERT_MSG(alloc_size > 0, "Attempting to allocate buffer of size 0");   
    m_memblob.resize(alloc_size);
    m_data = m_memblob.data();
    
    m_node.set_external(m_schema,m_data);
    return this;
}

void DataBuffer::Info(Node &n) const
{
    n["uid"].set(m_uid);
    n["descriptor"].set(m_schema.to_json());
    n["node"].set(m_node.to_json());
}


void DataBuffer::Print() const
{
    Node n;
    Info(n);
    n.print();
}



void DataBuffer::AttachView( DataView* view )
{
    m_views.push_back( view );
}


void DataBuffer::DetachView( DataView* view )
{
    //Find new end iterator
    std::vector<DataView*>::iterator pos = std::remove(m_views.begin(),
                                                       m_views.end(),
                                                       view);
    // check if pos is ok?
    //Erase the "removed" elements.
    m_views.erase(pos, m_views.end());
}


} /* namespace Datastore */



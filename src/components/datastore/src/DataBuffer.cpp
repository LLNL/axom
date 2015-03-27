
#include "DataBuffer.hpp"
#include "DataGroup.hpp"
#include "DataView.hpp"

namespace DataStoreNS
{


DataBuffer::DataBuffer( const IDType uid ) :
    m_uid(uid),
    m_ViewContainer(),
    m_data(nullptr),
    m_memblob(),
    m_node(),
    m_schema()
{
    
}


DataBuffer::DataBuffer(const DataBuffer& source ) :
    m_uid(source.m_uid),
    m_ViewContainer(source.m_ViewContainer),
    m_data(source.m_data),
    m_memblob(source.m_memblob),
    m_node(source.m_node),
    m_schema(source.m_schema)
        
{

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

    if(alloc_size == 0)
    {
        // ?
        throw std::exception();
    }
            
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



void DataBuffer::AttachView( DataView* dataView )
{
  m_ViewContainer.insert( dataView );
}


void DataBuffer::DetachView( DataView* dataView )
{
  m_ViewContainer.erase( dataView );
}


} /* namespace Datastore */



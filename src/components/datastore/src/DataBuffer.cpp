
#include "DataBuffer.hpp"
#include "DataGroup.hpp"
#include "DataView.hpp"

namespace DataStoreNS
{


DataBuffer::DataBuffer( const IDType uid ) :
    m_uid(uid),
    m_stringDescriptor(),
    m_ViewContainer(),
    m_data(nullptr),
    m_memblob(),
    m_node(),
    m_schema()
{
    
}

DataBuffer::DataBuffer( const IDType uid,
                        const std::string& stringDescriptor ) :
    m_uid(uid),
    m_stringDescriptor(stringDescriptor),
    m_ViewContainer(),
    m_data(nullptr),
    m_memblob(),
    m_node(),
    m_schema()
{
      
}

DataBuffer::DataBuffer(const DataBuffer& source ) :
    m_uid(source.m_uid),
    m_stringDescriptor(source.m_stringDescriptor),
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



/// init calls set descriptor, allocate, and apply descriptor  
DataBuffer* DataBuffer::Init(const Schema &schema)
{
    SetDescriptor(schema);
    Allocate();
    ApplyDescriptor();
    return this;
}

DataBuffer* DataBuffer::Init(const DataType &dtype)
{
    SetDescriptor(dtype);
    Allocate();
    ApplyDescriptor();
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
    
    ReconcileDataViews();
    return this;
}

void DataBuffer::Print(Node &n) const
{
    n["DataBuffer/descriptor"] = m_schema.to_json();
    n["DataBuffer/node"] = m_node.to_json();
}

void DataBuffer::Print() const
{
    Node n;
    Print(n);
    n.print();
}



DataBuffer* DataBuffer::ApplyDescriptor()
{
    m_node.set_external(m_schema,m_data);
    return this;
}


void DataBuffer::AttachView( DataView* dataView )
{
  m_ViewContainer.insert( dataView );
}


void DataBuffer::DetachView( DataView* dataView )
{
  m_ViewContainer.erase( dataView );
}

void DataBuffer::ReconcileDataViews()
{
  for( ViewContainerType::iterator iterView=m_ViewContainer.begin() ;
       iterView != m_ViewContainer.end() ; ++iterView )
  {
    (*iterView)->ReconcileWithBuffer();
  }
}

} /* namespace Datastore */

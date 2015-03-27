/*
 * DataSet.cpp
 *
 *  Created on: Dec 1, 2014
 *      Author: settgast
 */

#include "DataView.hpp"
#include "DataGroup.hpp"
#include "DataStore.hpp"
#include "DataBuffer.hpp"

namespace DataStoreNS
{


DataView::DataView( const std::string& name,
                    DataGroup* const parentGroup,
                    DataBuffer* const dataBuffer ):
    m_name(name),
    m_parentGroup(parentGroup),
    m_dataBuffer(dataBuffer),
    m_schema(),
    m_node()
{

}

DataView::DataView( const std::string& name,
                    DataGroup* const parentGroup) :
  m_name(name),
  m_parentGroup(parentGroup),
  m_dataBuffer(nullptr),
  m_schema(),
  m_node()
{
    m_dataBuffer = parentGroup->GetDataStore()->CreateBuffer();
    m_dataBuffer->AttachView(this);
}

DataView::DataView(const DataView& source ) :
    m_name(source.m_name),
    m_parentGroup(source.m_parentGroup),
    m_dataBuffer(source.m_dataBuffer),
    m_schema(source.m_schema),
    m_node(source.m_node)
{
    
}


DataView::~DataView()
{
    
}


DataView* DataView::Apply()
{
    m_node.set_external(m_schema,m_dataBuffer->GetData());
    return this;
}

DataView* DataView::Apply(const Schema &schema)
{
    m_schema.set(schema);
    Apply();
    return this;
}

DataView* DataView::Apply(const DataType &dtype)
{
    m_schema.set(dtype);
    Apply();
    return this;
}


DataView* DataView::Declare(const Schema &schema)
{
    m_schema.set(schema);
    return this;
}

DataView* DataView::Declare(const DataType &dtype)
{
    m_schema.set(dtype);
    return this;
}

DataView* DataView::Allocate()
{
    // we only force alloc if there is a 1-1 between the view and buffer
    if(m_dataBuffer->CountViews() != 1)
    {
        throw std::exception();
    }
    
    m_dataBuffer->Allocate(m_schema);
    Apply();
    return this;
}

DataView* DataView::Allocate(const Schema &schema)
{
    Declare(schema);
    Allocate();
    Apply();
    return this;
}

DataView* DataView::Allocate(const DataType &dtype)
{
    
    Declare(dtype);
    Allocate();
    Apply();
    return this;
}



void DataView::Info(Node &n) const
{
    n["name"] = m_name;
    n["descriptor"] = m_schema.to_json();
    n["node"] = m_node.to_json();
}

void DataView::Print() const
{
    Node n;
    Info(n);
    n.print();
}


} /* namespace Datastore */

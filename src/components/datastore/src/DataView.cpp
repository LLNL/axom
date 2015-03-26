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
    m_node(),
    m_schema()
{

}

DataView::DataView( const std::string& name,
                    DataGroup* const parentGroup,
                    DataStore* const dataStore ) :
  m_name(name),
  m_parentGroup(parentGroup),
  m_dataBuffer(nullptr),
  m_node(),
  m_schema()
{
    m_dataBuffer = dataStore->CreateBuffer();
    m_dataBuffer->AttachView(this);
}

DataView::DataView(const DataView& source ) :
    m_name(source.m_name),
    m_parentGroup(source.m_parentGroup),
    m_dataBuffer(source.m_dataBuffer),
    m_node(source.m_node),
    m_schema(source.m_schema)
{
}


DataView::~DataView()
{
    
}


DataView* DataView::Init(const Schema &schema)
{
    SetDescriptor(schema);
    m_dataBuffer->Init(m_schema);
    ApplyDescriptor();
    return this;
}

DataView* DataView::Init(const DataType &dtype)
{
    SetDescriptor(dtype);
    m_dataBuffer->Init(m_schema);
    ApplyDescriptor();
    return this;
}


DataView* DataView::ApplyDescriptor()
{
    m_node.set_external(m_schema,m_dataBuffer->GetData());
    return this;
}


DataView* DataView::SetDescriptor(const Schema &schema)
{
    m_schema.set(schema);
    return this;
}

DataView* DataView::SetDescriptor(const DataType &dtype)
{
    m_schema.set(dtype);
    return this;
}



void DataView::ReconcileWithBuffer()
{
    /// TODO:
}


void DataView::Print(Node &n) const
{
    n["DataView/name"] = m_name;
    n["DataView/descriptor"] = m_schema.to_json();
}


void DataView::Print() const
{
    Node n;
    Print(n);
    n.print();
}


} /* namespace Datastore */

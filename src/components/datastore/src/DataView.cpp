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


DataView* DataView::Allocate()
{
    // this implies a push alloc?
    m_dataBuffer->SetDescriptor(m_schema);
    m_dataBuffer->Allocate();
    m_dataBuffer->ApplyDescriptor();
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
  m_schema.set(m_dataBuffer->GetDescriptor());
  ApplyDescriptor();
}

} /* namespace Datastore */

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
                    DataGroup* const parent,
                    DataBuffer* const buffer ):
    m_name(name),
    m_group(parent),
    m_buffer(buffer),
    m_schema(),
    m_node(),
    m_applied(false)
{

}

DataView::DataView( const std::string& name,
                    DataGroup* const parent) :
  m_name(name),
  m_group(parent),
  m_buffer(nullptr),
  m_schema(),
  m_node(),
  m_applied(false)
{
    m_buffer = parent->GetDataStore()->CreateBuffer();
    m_buffer->AttachView(this);
}

/// Note: we can't simply set the group pointer here
DataView::DataView( const DataView& source ) :
    m_name(source.m_name),
    m_group(nullptr),
    m_buffer(source.m_buffer),
    m_schema(source.m_schema),
    m_node(source.m_node),
    m_applied(source.m_applied)
{
    throw std::exception();
}


DataView::~DataView()
{
    
}

DataView* DataView::Apply()
{
    m_node.set_external(m_schema,m_buffer->GetData());
    m_applied = true;
    return this;
}

DataView* DataView::Apply(const Schema &schema)
{
    Declare(schema);
    Apply();
    return this;
}

DataView* DataView::Apply(const DataType &dtype)
{
    Declare(dtype);
    Apply();
    return this;
}


DataView* DataView::Declare(const Schema &schema)
{
    m_schema.set(schema);
    m_applied = false;
    return this;
}

DataView* DataView::Declare(const DataType &dtype)
{
    m_schema.set(dtype);
    m_applied = false;
    return this;
}

DataView* DataView::Allocate()
{
    // we only force alloc if there is a 1-1 between the view and buffer
    if(m_buffer->CountViews() != 1)
    {
        throw std::exception();
    }
    
    m_buffer->Allocate(m_schema);
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
    n["applied"] = m_applied;
}

void DataView::Print() const
{
    Node n;
    Info(n);
    n.print();
}


} /* namespace Datastore */

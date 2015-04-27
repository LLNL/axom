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

#include "Utilities.hpp"

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
    m_applied(false),
    m_opaque(false)
{

}

DataView::DataView( const std::string& name,
                    DataGroup* const parent,
                    void* opaque) :
  m_name(name),
  m_group(parent),
  m_buffer(nullptr),
  m_schema(),
  m_node(),
  m_applied(false),
  m_opaque(true)
{
    // todo, conduit should provide a check for if uint64 is a
    // good enough type to rep void *
    GetNode().set((conduit::uint64)opaque);
}


DataView::~DataView()
{
    if(m_buffer != nullptr)
    {
        m_buffer->DetachView(this);
    }
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
    ATK_ASSERT_MSG( m_buffer->CountViews() == 1, \
                      "Can only allocate from a view if it's the only view associated with its buffer");
    
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

void *DataView::GetOpaque() const
{
    // if(!m_opaque) error?
    return (void*)(GetNode().as_uint64());
}

void DataView::Info(Node &n) const
{
    n["name"] = m_name;
    n["descriptor"] = m_schema.to_json();
    n["node"] = m_node.to_json();
    n["applied"] = m_applied;
    n["opaque"] = m_opaque;
}

void DataView::Print() const
{
    Node n;
    Info(n);
    n.print();
}


} /* namespace Datastore */

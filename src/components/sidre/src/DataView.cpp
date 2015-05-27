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
 * \brief   Implementation file for DataView class.
 *
 ******************************************************************************
 */


// Associated header file
#include "DataView.hpp"

// Other toolkit project headers
#include "common/CommonTypes.hpp"

// SiDRe project headers
#include "DataBuffer.hpp"
#include "DataGroup.hpp"
#include "DataStore.hpp"


namespace asctoolkit
{
namespace sidre
{


/*
*************************************************************************
*
* Declare a data view as a Conduit schema.
*
*************************************************************************
*/
DataView* DataView::declare(const Schema& schema)
{
    m_schema.set(schema);
    m_is_applied = false;
    return this;
}

/*
*************************************************************************
*
* Declare a data view as a pre-defined Conduit data type.
*
*************************************************************************
*/
DataView* DataView::declare(const DataType& dtype)
{
    m_schema.set(dtype);
    m_is_applied = false;
    return this;
}

/*
*************************************************************************
*
* Allocate data for a DataView, previously declared.
*
*************************************************************************
*/
DataView* DataView::allocate()
{
    // we only force alloc if there is a 1-1 between the view and buffer
    ATK_ASSERT_MSG( m_data_buffer->getNumViews() == 1, \
                    "Data allocation on a view allowed only if it's the only view associated with its buffer");
    
    m_data_buffer->allocate(m_schema);
    apply();
    return this;
}

/*
*************************************************************************
*
* Declare a data view as a Conduit schema and allocate the data.
*
*************************************************************************
*/
DataView* DataView::allocate(const Schema& schema)
{
    declare(schema);
    allocate();
    apply();
    return this;
}

/*
*************************************************************************
*
* Declare a data view as a Conduit pre-defined data type and allocate the data.
*
*************************************************************************
*/
DataView* DataView::allocate(const DataType& dtype)
{
    
    declare(dtype);
    allocate();
    apply();
    return this;
}


/*
*************************************************************************
*
* Reallocate the data view's buffer according to a Conduit schema.
*
*************************************************************************
*/
DataView* DataView::reallocate(const Schema& schema)
{
    // in these cases the view does not have a buffer
    ATK_ASSERT_MSG( !isExternal() && !isOpaque(),
                   "Attempting to reallocate an external or opaque view");
    
    // following alloc,  we only force realloc if there is a 1-1 between 
    // the view and its buffer.
    ATK_ASSERT_MSG( m_data_buffer->getNumViews() == 1, \
                    "Data reallocation on a view allowed only if it's the only view associated with its buffer");
    
    declare(schema);
    m_data_buffer->reallocate(m_schema);
    apply();
    return this;
}

/*
*************************************************************************
*
* Reallocate the data view's buffer according to a Conduit data type.
*
*************************************************************************
*/
DataView* DataView::reallocate(const DataType& dtype)
{
    
    // in these cases the view does not have a buffer
    ATK_ASSERT_MSG( !isExternal() && !isOpaque(),
                   "Attempting to reallocate an external or opaque view");

    // following alloc,  we only force realloc if there is a 1-1 between 
    // the view and its buffer.
    ATK_ASSERT_MSG( m_data_buffer->getNumViews() == 1, \
                    "Data reallocation on a view allowed only if it's the only view associated with its buffer");
    
    declare(dtype);
    m_data_buffer->reallocate(m_schema);
    apply();
    return this;
}


/*
*************************************************************************
*
* Apply a previously declared data view to data held in the buffer.
*
*************************************************************************
*/
DataView* DataView::apply()
{
    m_node.set_external(m_schema, m_data_buffer->getData());
    m_is_applied = true;
    return this;
}

/*
*************************************************************************
*
* Apply a Conduit Schema to data held in the buffer.
*
*************************************************************************
*/
DataView* DataView::apply(const Schema& schema)
{
    declare(schema);
    apply();
    return this;
}

/*
*************************************************************************
*
* Apply a Conduit pre-defined data type description to data held in the buffer.
*
*************************************************************************
*/
DataView* DataView::apply(const DataType &dtype)
{
    declare(dtype);
    apply();
    return this;
}

/*
*************************************************************************
*
* Return void* pointer to opaque data.
*
*************************************************************************
*/
void *DataView::getOpaque() const
{
    // if (!m_is_opaque) error?
    return (void*)(getNode().as_uint64());
}


/*
*************************************************************************
*
* Copy data view description to given Conduit node.
*
*************************************************************************
*/
void DataView::info(Node &n) const
{
    n["name"] = m_name;
    n["schema"] = m_schema.to_json();
    n["node"] = m_node.to_json();
    n["is_opaque"] = m_is_opaque;
    n["is_external"] = m_is_external;
    n["is_applied"] = m_is_applied;
}


/*
*************************************************************************
*
* Print JSON description of data view to stdout.
*
*************************************************************************
*/
void DataView::print() const
{
    Node n;
    info(n);
    n.print();
}


/*
*************************************************************************
*
* PRIVATE ctor for DataView associated with DataBuffer.
*
*************************************************************************
*/
DataView::DataView( const std::string& name,
                    DataGroup* const owning_group,
                    DataBuffer* const data_buffer) 
:   m_name(name),
    m_owning_group(owning_group),
    m_data_buffer(data_buffer),
    m_schema(),
    m_node(),
    m_is_opaque(false),
    m_is_external(false),
    m_is_applied(false)
{

}

/*
*************************************************************************
*
* PRIVATE ctor for DataView associated with opaque data.
*
*************************************************************************
*/
DataView::DataView( const std::string& name,
                    DataGroup* const owning_group,
                    void* opaque_ptr) 
: m_name(name),
  m_owning_group(owning_group),
  m_data_buffer(ATK_NULLPTR),
  m_schema(),
  m_node(),
  m_is_opaque(true),
  m_is_external(false),
  m_is_applied(false)
{
    // todo, conduit should provide a check for if uint64 is a
    // good enough type to rep void *
    m_node.set((conduit::uint64)opaque_ptr);
}

/*
*************************************************************************
*
* PRIVATE ctor for DataView associated with external data described by
*         given data type.
*
*************************************************************************
*/
DataView::DataView( const std::string& name,
                    DataGroup* const owning_group,
                    void* external_data,
                    const DataType& dtype) 
: m_name(name),
  m_owning_group(owning_group),
  m_data_buffer(ATK_NULLPTR),
  m_schema(),
  m_node(),
  m_is_opaque(false),
  m_is_external(true),
  m_is_applied(false)
{
    m_schema.set(dtype);
    m_node.set_external(m_schema, external_data);
    m_is_applied = true;
}

/*
*************************************************************************
*
* PRIVATE ctor for DataView associated with external data described by
*         given schema.
*
*************************************************************************
*/
DataView::DataView( const std::string& name,
                    DataGroup* const owning_group,
                    void* external_data,
                    const Schema& schema)
: m_name(name),
  m_owning_group(owning_group),
  m_data_buffer(ATK_NULLPTR),
  m_schema(),
  m_node(),
  m_is_opaque(false),
  m_is_external(true),
  m_is_applied(false)
{
    m_schema.set(schema);
    m_node.set_external(m_schema, external_data);
    m_is_applied = true;
}


/*
*************************************************************************
*
* PRIVATE dtor.
*
*************************************************************************
*/
DataView::~DataView()
{
    if (m_data_buffer != ATK_NULLPTR)
    {
        m_data_buffer->detachView(this);
    }
}


} /* end namespace sidre */
} /* end namespace asctoolkit */


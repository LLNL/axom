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
#include "slic/slic.hpp"

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
 * Declare a data view from sidre type and length.
 *
 *************************************************************************
 */
DataView * DataView::declare(TypeID type, SidreLength len)
{
  SLIC_ASSERT_MSG( !isOpaque(),
                  "Cannot call declare on an opaque view");
  SLIC_ASSERT_MSG(len >= 0, "Must declare number of elements in view >=0");

  if ( !isOpaque() && len >= 0 ) 
  {
    DataType dtype = conduit::DataType::default_dtype(type);
    dtype.set_number_of_elements(len);

    m_schema.set(dtype);
    m_is_applied = false;
  }
  return this;
}

/*
 *************************************************************************
 *
 * Declare a data view as a Conduit schema.
 *
 *************************************************************************
 */
DataView * DataView::declare(const Schema& schema)
{
  SLIC_ASSERT_MSG( !isOpaque(),
                  "Cannot call declare on an opaque view");

  if ( !isOpaque() ) 
  {
    m_schema.set(schema);
    m_is_applied = false;
  }
  return this;
}

/*
 *************************************************************************
 *
 * Declare a data view as a pre-defined Conduit data type.
 *
 *************************************************************************
 */
DataView * DataView::declare(const DataType& dtype)
{
  SLIC_ASSERT_MSG( !isOpaque(),
                  "Cannot call declare on an opaque view");

  if ( !isOpaque() ) 
  {
    m_schema.set(dtype);
    m_is_applied = false;
  }
  return this;
}

/*
 *************************************************************************
 *
 * Allocate data for a DataView, previously declared.
 *
 *************************************************************************
 */
DataView * DataView::allocate()
{
  SLIC_ASSERT_MSG( !isOpaque(),
                  "Cannot call allocate an external or opaque view");
  SLIC_ASSERT_MSG( m_data_buffer->getNumViews() == 1, \
                  "Data allocation on a view allowed only if it's the only view associated with its buffer");

  if ( !isOpaque() && m_data_buffer->getNumViews() == 1 )
  {
    m_data_buffer->allocate(m_schema);
    apply();
  }
  return this;
}

/*
 *************************************************************************
 *
 * Declare a data view using sidre type and length.
 *
 *************************************************************************
 */
DataView * DataView::allocate( TypeID type, SidreLength len)
{
  SLIC_ASSERT_MSG( !isOpaque(),
                  "Attempting to allocate an external or opaque view");
  SLIC_ASSERT_MSG(len >= 0, "Must allocate number of elements in view >=0");
  SLIC_ASSERT_MSG( m_data_buffer->getNumViews() == 1, \
                  "Data allocation on a view allowed only if it's the only view associated with its buffer");

  if ( !isOpaque() && len >= 0 && m_data_buffer->getNumViews() == 1 ) 
  {
    declare(type, len);
    allocate();
    apply();
  }
  return this;
}

/*
 *************************************************************************
 *
 * Declare a data view as a Conduit schema and allocate the data.
 *
 *************************************************************************
 */
DataView * DataView::allocate(const Schema& schema)
{
  SLIC_ASSERT_MSG( !isOpaque(),
                  "Attempting to allocate an external or opaque view");
  SLIC_ASSERT_MSG( m_data_buffer->getNumViews() == 1, \
                  "Data allocation on a view allowed only if it's the only view associated with its buffer");

  if ( !isOpaque() && m_data_buffer->getNumViews() == 1 )
  {
    declare(schema);
    allocate();
    apply();
  }
  return this;
}

/*
 *************************************************************************
 *
 * Declare a data view as a Conduit pre-defined data type and allocate the data.
 *
 *************************************************************************
 */
DataView * DataView::allocate(const DataType& dtype)
{
  SLIC_ASSERT_MSG( !isOpaque(),
                  "Attempting to allocate an external or opaque view");
  SLIC_ASSERT_MSG( m_data_buffer->getNumViews() == 1, \
                  "Data allocation on a view allowed only if it's the only view associated with its buffer");

  if ( !isOpaque() && m_data_buffer->getNumViews() == 1 )
  {
    declare(dtype);
    allocate();
    apply();
  }
  return this;
}


/*
 *************************************************************************
 *
 * Reallocate the data view's buffer according to a Conduit schema.
 *
 *************************************************************************
 */
DataView * DataView::reallocate(TypeID type, SidreLength len)
{
  SLIC_ASSERT_MSG( !isOpaque(),
                  "Attempting to reallocate an external or opaque view");
  SLIC_ASSERT_MSG(len >= 0, "Must re-allocate number of elements in view >=0");
  SLIC_ASSERT_MSG( m_data_buffer->getNumViews() == 1, \
                  "Data reallocation on a view allowed only if it's the only view associated with its buffer");

  if ( !isOpaque() && len >= 0 && m_data_buffer->getNumViews() == 1 )
  {
    declare(type, len);
    m_data_buffer->reallocate(type, len);
    apply();
  }
  return this;
}

/*
 *************************************************************************
 *
 * Reallocate the data view's buffer according to a Conduit schema.
 *
 *************************************************************************
 */
DataView * DataView::reallocate(const Schema& schema)
{
  SLIC_ASSERT_MSG( !isOpaque(),
                  "Attempting to reallocate an external or opaque view");
  SLIC_ASSERT_MSG( m_data_buffer->getNumViews() == 1, \
                  "Data reallocation on a view allowed only if it's the only view associated with its buffer");

  if ( !isOpaque() && m_data_buffer->getNumViews() == 1 )
  {
    declare(schema);
    m_data_buffer->reallocate(m_schema);
    apply();
  }
  return this;
}

/*
 *************************************************************************
 *
 * Reallocate the data view's buffer according to a Conduit data type.
 *
 *************************************************************************
 */
DataView * DataView::reallocate(const DataType& dtype)
{
  SLIC_ASSERT_MSG( !isOpaque(),
                  "Attempting to reallocate an external or opaque view");
  SLIC_ASSERT_MSG( m_data_buffer->getNumViews() == 1, \
                  "Data reallocation on a view allowed only if it's the only view associated with its buffer");

  if ( !isOpaque() && m_data_buffer->getNumViews() == 1 )
  {
    declare(dtype);
    m_data_buffer->reallocate(m_schema);
    apply();
  }
  return this;
}


/*
 *************************************************************************
 *
 * Apply a previously declared data view to data held in the buffer.
 *
 *************************************************************************
 */
DataView * DataView::apply()
{
  SLIC_ASSERT_MSG( !isOpaque(),
                  "Cannot call apply() on an opaque view");
 
  if ( !isOpaque() )
  { 
    m_node.set_external(m_schema, m_data_buffer->getData());
    m_is_applied = true;
  }
  return this;
}

/*
 *************************************************************************
 *
 * Apply a Conduit Schema to data held in the buffer.
 *
 *************************************************************************
 */
DataView * DataView::apply(const Schema& schema)
{
  SLIC_ASSERT_MSG( !isOpaque(),
                  "Cannot call apply() on an opaque view");
 
  if ( !isOpaque() )
  { 
    declare(schema);
    apply();
  }
  return this;
}

/*
 *************************************************************************
 *
 * Apply a Conduit pre-defined data type description to data held in the buffer.
 *
 *************************************************************************
 */
DataView * DataView::apply(const DataType &dtype)
{
  SLIC_ASSERT_MSG( !isOpaque(),
                  "Cannot call apply() on an opaque view");
 
  if ( !isOpaque() )
  { 
    declare(dtype);
    apply();
  }
  return this;
}

/*
 *************************************************************************
 *
 * Return void* pointer to buffer data.
 *
 *************************************************************************
 */
void * DataView::getDataPointer() const
{
  if ( isOpaque() ) {
      return (void *)(getNode().as_uint64());
  } else {
      return m_data_buffer->getData();
  }
}

/*
 *************************************************************************
 *
 * Return void* pointer to opaque data.
 *
 *************************************************************************
 */
void * DataView::getOpaque() const
{
  if ( isOpaque() ) 
  {
    return (void *)(getNode().as_uint64());
  } 
  else 
  {
    return ATK_NULLPTR; 
  }
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
  print(std::cout);
}

/*
 *************************************************************************
 *
 * Print JSON description of data view to an  ostream.
 *
 *************************************************************************
 */
void DataView::print(std::ostream& os) const
{
  Node n;
  info(n);
  /// TODO: after conduit update, use new ostream variant of to_json.
  std::ostringstream oss;
  n.json_to_stream(oss);
  os << oss.str();
}

/*
 *************************************************************************
 *
 * PRIVATE ctor for DataView associated with DataBuffer.
 *
 *************************************************************************
 */
DataView::DataView( const std::string& name,
                    DataGroup * const owning_group,
                    DataBuffer * const data_buffer)
  :   m_name(name),
  m_owning_group(owning_group),
  m_data_buffer(data_buffer),
  m_schema(),
  m_node(),
  m_is_opaque(false),
  m_is_applied(false)
{}

/*
 *************************************************************************
 *
 * PRIVATE ctor for DataView associated with opaque data.
 *
 *************************************************************************
 */
DataView::DataView( const std::string& name,
                    DataGroup * const owning_group,
                    void * opaque_ptr)
  : m_name(name),
  m_owning_group(owning_group),
  m_data_buffer(ATK_NULLPTR),
  m_schema(),
  m_node(),
  m_is_opaque(true),
  m_is_applied(false)
{
  // todo, conduit should provide a check for if uint64 is a
  // good enough type to rep void *
  m_node.set((conduit::uint64)opaque_ptr);
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

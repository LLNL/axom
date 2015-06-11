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
 * \brief   Implementation file for DataBuffer class.
 *
 ******************************************************************************
 */


// Associated header file
#include "DataBuffer.hpp"

// Standard C++ headers
#include <algorithm>

// Other CS Toolkit headers
#include "common/CommonTypes.hpp"
#include "common/Utilities.hpp"

// SiDRe project headers
#include "DataGroup.hpp"
#include "DataView.hpp"
#include "SidreTypes.hpp"


namespace asctoolkit
{
namespace sidre
{

/*
 *************************************************************************
 *
 * Set buffer to externally-owned data.
 *
 *************************************************************************
 */
DataBuffer * DataBuffer::declareExternal(void * external_data,
                                         const Schema& schema)
{
  ATK_ASSERT_MSG( m_data == ATK_NULLPTR,
                  "Attempting to declare buffer external, but buffer has already been allocated" );
  ATK_ASSERT_MSG( external_data != ATK_NULLPTR,
                  "Attempting to set buffer to null external data" );
  m_schema.set(schema);
  m_data = external_data;
  m_node.set_external(m_schema, m_data);
  m_is_data_external = true;
  return this;
}

/*
 *************************************************************************
 *
 * Set buffer to externally-owned data.
 *
 *************************************************************************
 */
DataBuffer * DataBuffer::declareExternal(void * external_data,
                                         const DataType& dtype)
{
  ATK_ASSERT_MSG( m_data == ATK_NULLPTR,
                  "Attempting to declare buffer external, but buffer has already been allocated" );
  ATK_ASSERT_MSG( external_data != ATK_NULLPTR,
                  "Attempting to set buffer to null external data" );
  m_schema.set(dtype);
  m_data = external_data;
  m_node.set_external(m_schema, m_data);
  m_is_data_external = true;
  return this;
}

/*
 *************************************************************************
 *
 * Declare and allocate data described using type enum and length.
 *
 *************************************************************************
 */
DataBuffer * DataBuffer::declare(const TypeID type, const SidreLength len)
{
  ATK_ASSERT_MSG(len >= 0, "Bad Length");
  DataType dtype = conduit::DataType::default_dtype(type);
  dtype.set_number_of_elements(len);
  m_schema.set(dtype);
  return this;
}

/*
 *************************************************************************
 *
 * Allocate data previously declared.
 *
 *************************************************************************
 */
DataBuffer * DataBuffer::allocate()
{
  ATK_ASSERT_MSG( !m_is_data_external,
                  "Attempting to allocate buffer holding external data");

  if ( !m_is_data_external )
  {
    // cleanup old data
    cleanup();
    std::size_t alloc_size = m_schema.total_bytes();
    ATK_ASSERT_MSG(alloc_size > 0,
                   "Attempting to allocate buffer of size 0");
    m_data = allocateBytes(alloc_size);
    m_node.set_external(m_schema,m_data);
  }

  return this;
}

/*
 *************************************************************************
 *
 * Declare and allocate data described using type and length.
 *
 *************************************************************************
 */
DataBuffer * DataBuffer::allocate(const TypeID type, const SidreLength len)
{
  ATK_ASSERT_MSG( !m_is_data_external,
                  "Attempting to allocate buffer holding external data");

  if ( !m_is_data_external )
  {
    declare(type, len);
    allocate();
  }

  return this;
}

/*
 *************************************************************************
 *
 * Declare and allocate data described using a Conduit schema.
 *
 *************************************************************************
 */
DataBuffer * DataBuffer::allocate(const Schema& schema)
{
  ATK_ASSERT_MSG( !m_is_data_external,
                  "Attempting to allocate buffer holding external data");

  if ( !m_is_data_external )
  {
    declare(schema);
    allocate();
  }

  return this;
}

/*
 *************************************************************************
 *
 * Declare and allocate data described using a Conduit pre-defined data type.
 *
 *************************************************************************
 */
DataBuffer * DataBuffer::allocate(const DataType& dtype)
{
  ATK_ASSERT_MSG( !m_is_data_external,
                  "Attempting to allocate buffer holding external data");

  if ( !m_is_data_external )
  {
    declare(dtype);
    allocate();
  }

  return this;
}

/*
 *************************************************************************
 *
 * Reallocate data using a Sidre type and length.
 *
 *************************************************************************
 */
DataBuffer * DataBuffer::reallocate( const TypeID type, const SidreLength len)
{
  ATK_ASSERT_MSG( !m_is_data_external,
                  "Attempting to re-allocate buffer holding external data");

  if ( !m_is_data_external )
  {
    DataType dtype = conduit::DataType::default_dtype(type);
    dtype.set_number_of_elements(len);

    Schema s(dtype);
    reallocate(s);
  }

  return this;
}

/*
 *************************************************************************
 *
 * Reallocate data using a Conduit schema.
 *
 *************************************************************************
 */
DataBuffer * DataBuffer::reallocate(const Schema& schema)
{
  ATK_ASSERT_MSG( !m_is_data_external,
                  "Attempting to re-allocate buffer holding external data");

  if ( !m_is_data_external )
  {

    //  make sure realloc actually makes sense
    ATK_ASSERT_MSG(m_data != ATK_NULLPTR,
                   "Attempting to reallocate an unallocated buffer");

    std::size_t realloc_size = schema.total_bytes();
    ATK_ASSERT_MSG(realloc_size > 0,
                   "Attempting to reallocate buffer to size 0");

    void * realloc_data = allocateBytes(realloc_size);

    // use conduit to get data from old to new schema.
    Node n;
    n.set_external(schema, realloc_data);
    // use conduit update, need more error checking.
    n.update(m_node);

    // cleanup old data
    cleanup();

    // set the buffer to use the new schema
    m_schema = schema;

    // let the buffer hold the new data
    m_data = realloc_data;

    // update the buffer's Conduit Node
    m_node.set_external(m_schema, m_data);
  }

  return this;
}

/*
 *************************************************************************
 *
 * Reallocate data using a basic Conduit data type.
 *
 *************************************************************************
 */
DataBuffer * DataBuffer::reallocate(const DataType& dtype)
{
  ATK_ASSERT_MSG( !m_is_data_external,
                  "Attempting to re-allocate buffer holding external data");

  if ( !m_is_data_external )
  {
    Schema s(dtype);
    reallocate(s);
  }

  return this;
}


/*
 *************************************************************************
 *
 * Copy data buffer description to given Conduit node.
 *
 *************************************************************************
 */
void DataBuffer::info(Node &n) const
{
  n["index"].set(m_index);
  n["is_data_external"].set(m_is_data_external);
  n["schema"].set(m_schema.to_json());
  n["node"].set(m_node.to_json());
}

/*
 *************************************************************************
 *
 * Print JSON description of data buffer to stdout.
 *
 *************************************************************************
 */
void DataBuffer::print() const
{
  print(std::cout);
}

/*
 *************************************************************************
 *
 * Print JSON description of data buffer to an ostream.
 *
 *************************************************************************
 */
void DataBuffer::print(std::ostream& os) const
{
  Node n;
  info(n);
  /// TODO: after conduit update, use new ostream variant of to_json.
  std::ostringstream oss;
  n.to_pure_json(oss);
  os << oss.str();
}



/*
 *************************************************************************
 *
 * PRIVATE ctor taking unique index.
 *
 *************************************************************************
 */
DataBuffer::DataBuffer( IndexType index )
  : m_index(index),
  m_views(),
  m_data(ATK_NULLPTR),
  m_node(),
  m_schema(),
  m_is_data_external(false)
{}


/*
 *************************************************************************
 *
 * PRIVATE copy ctor.
 *
 *************************************************************************
 */
DataBuffer::DataBuffer(const DataBuffer& source )
  : m_index(source.m_index),
  m_views(source.m_views),
  m_data(source.m_data),
  m_node(source.m_node),
  m_schema(source.m_schema),
  m_is_data_external(source.m_is_data_external)
{
// disallow?
}


/*
 *************************************************************************
 *
 * PRIVATE dtor.
 *
 *************************************************************************
 */
DataBuffer::~DataBuffer()
{
  cleanup();
}


/*
 *************************************************************************
 *
 * PRIVATE method to attach data view.
 *
 *************************************************************************
 */
void DataBuffer::attachView( DataView * view )
{
  m_views.push_back( view );
}


/*
 *************************************************************************
 *
 * PRIVATE method to detach data view.
 *
 *************************************************************************
 */
void DataBuffer::detachView( DataView * view )
{
  //Find new end iterator
  std::vector<DataView *>::iterator pos = std::remove(m_views.begin(),
                                                      m_views.end(),
                                                      view);
  // check if pos is ok?
  //Erase the "removed" elements.
  m_views.erase(pos, m_views.end());
}

/*
 *************************************************************************
 *
 * PRIVATE cleanup
 *
 *************************************************************************
 */
void DataBuffer::cleanup()
{
  // cleanup alloced data
  if ( m_data != ATK_NULLPTR && !m_is_data_external )
  {
    releaseBytes(m_data);
    m_data = ATK_NULLPTR;
  }
}

/*
 *************************************************************************
 *
 * PRIVATE allocateBytes
 *
 *************************************************************************
 */
void * DataBuffer::allocateBytes(std::size_t num_bytes)
{
  ATK_ASSERT_MSG(num_bytes > 0,
                 "Attempting to allocate 0 bytes");

  char * data = new char[num_bytes];
  return ((void *)data);
}

/*
 *************************************************************************
 *
 * PRIVATE releaseBytes
 *
 *************************************************************************
 */
void DataBuffer::releaseBytes(void * ptr)
{
  if ( !m_is_data_external )
  {
    delete [] ((char *)ptr);
  }
}



} /* end namespace sidre */
} /* end namespace asctoolkit */

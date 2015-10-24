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
#include "slic/slic.hpp"

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
 * Return total number of bytes associated with this DataBuffer object.
 *
 *************************************************************************
 */
size_t DataBuffer::getTotalBytes() const
{
  static size_t bytes_per_item[] = {
    0, // CONDUIT_EMPTY_T
    0, // CONDUIT_OBJECT_T
    0, // CONDUIT_LIST_T
    1, // CONDUIT_INT8_T
    2, // CONDUIT_INT16_T
    4, // CONDUIT_INT32_T
    8, // CONDUIT_INT64_T
    1, // CONDUIT_UINT8_T
    2, // CONDUIT_UINT16_T
    4, // CONDUIT_UINT32_T
    8, // CONDUIT_UINT64_T
    4, // CONDUIT_FLOAT32_T
    8, // CONDUIT_FLOAT64_T
    1, // CONDUIT_CHAR8_STR_T
  };

  return bytes_per_item[m_type] * m_nitems;
}

/*
 *************************************************************************
 *
 * Return non-cost pointer to view with given index or null ptr.
 *
 *************************************************************************
 */
DataView * DataBuffer::getView( IndexType idx )
{
  SLIC_CHECK_MSG(hasView(idx), "no view exists with index == " << idx);

  if ( hasView(idx) )
  {
    return m_views[idx];
  }
  else
  {
    return ATK_NULLPTR;
  }
}


/*
 *************************************************************************
 *
 * Declare buffer to OWN data of given type and number of elements.
 *
 *************************************************************************
 */
DataBuffer * DataBuffer::declare(TypeID type, SidreLength len)
{
  SLIC_ASSERT_MSG(len >= 0, "Must declare number of elements >=0");

  if ( len >= 0 )
  {
    m_type = type;
    m_nitems = len;

    DataType dtype = conduit::DataType::default_dtype(type);
    dtype.set_number_of_elements(len);
    m_schema.set(dtype);
  }
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
  SLIC_ASSERT_MSG( !m_is_data_external,
                  "Attempting to allocate buffer holding external data");

  if ( !m_is_data_external )
  {
    // cleanup old data
    cleanup();
    std::size_t alloc_size = m_schema.total_bytes();
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
DataBuffer * DataBuffer::allocate(TypeID type, SidreLength len)
{
  SLIC_ASSERT_MSG(len >= 0, "Must allocate number of elements >=0");
  SLIC_ASSERT_MSG( !m_is_data_external,
                  "Attempting to allocate buffer holding external data");

  if ( len >= 0 && !m_is_data_external )
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
  SLIC_ASSERT_MSG( !m_is_data_external,
                  "Attempting to allocate buffer holding external data");

  if ( !m_is_data_external )
  {
    TypeID type = static_cast<TypeID>(schema.dtype().id());
    SidreLength nitems = schema.dtype().number_of_elements();
    declare(type, nitems);
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
  SLIC_ASSERT_MSG( !m_is_data_external,
                  "Attempting to allocate buffer holding external data");

  if ( !m_is_data_external )
  {
    TypeID type = static_cast<TypeID>(dtype.id());
    SidreLength nitems = dtype.number_of_elements();
    declare(type, nitems);
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
DataBuffer * DataBuffer::reallocate( TypeID type, SidreLength len)
{
  SLIC_ASSERT_MSG(len >= 0, "Must re-allocate number of elements >=0");
  SLIC_ASSERT_MSG( !m_is_data_external,
                  "Attempting to re-allocate buffer holding external data");

  if ( len >= 0 && !m_is_data_external )
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
  SLIC_ASSERT_MSG( !m_is_data_external,
                  "Attempting to re-allocate buffer holding external data");

  if ( !m_is_data_external )
  {
    TypeID type = static_cast<TypeID>(schema.dtype().id());
    SidreLength nitems = schema.dtype().number_of_elements();

    //  make sure realloc actually makes sense
    SLIC_ASSERT_MSG(m_data != ATK_NULLPTR,
                   "Attempting to reallocate an unallocated buffer");

    std::size_t realloc_size = schema.total_bytes();
    void * realloc_data = allocateBytes(realloc_size);

    // use conduit to get data from old to new schema.
    Node n;
    n.set_external(schema, realloc_data);
    // use conduit update, need more error checking.
    n.update(m_node);

    // XXX update m_type and m_nitems too
    declare(type, nitems);

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
  SLIC_ASSERT_MSG( !m_is_data_external,
                  "Attempting to re-allocate buffer holding external data");

  if ( !m_is_data_external )
  {
    TypeID type = static_cast<TypeID>(dtype.id());
    SidreLength nitems = dtype.number_of_elements();
    reallocate(type, nitems);
  }

  return this;
}


/*
 *************************************************************************
 *
 * Set buffer to externally-owned data.
 *
 *************************************************************************
 */
DataBuffer * DataBuffer::setExternalData(void * external_data)
{
  SLIC_ASSERT_MSG( external_data != ATK_NULLPTR, 
                  "Attempting to set buffer to external data given null pointer" );

  if ( external_data != ATK_NULLPTR )
  {
    m_data = external_data;
    m_node.set_external(m_schema, m_data);
    m_is_data_external = true;
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
  n.json_to_stream(oss);
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
  m_type(CONDUIT_EMPTY_T),
  m_nitems(0),
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
  m_type(CONDUIT_EMPTY_T),
  m_nitems(0),
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
  SLIC_ASSERT_MSG(num_bytes > 0,
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
    m_data = ATK_NULLPTR;
  }
}



} /* end namespace sidre */
} /* end namespace asctoolkit */

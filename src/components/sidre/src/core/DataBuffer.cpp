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
#include <cstring> //for memcpy

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
  return m_schema.total_bytes();
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
  if ( !hasView(idx) )
  {
    SLIC_CHECK_MSG(hasView(idx), "no view exists with index == " << idx);
    return ATK_NULLPTR;
  }

  return m_views[idx];
}


/*
 *************************************************************************
 *
 * Declare buffer to OWN data of given type and number of elements.
 *
 *************************************************************************
 */
DataBuffer * DataBuffer::declare(TypeID type, SidreLength num_elems)
{
  if ( num_elems < 0 )
  {
    SLIC_CHECK_MSG(num_elems >= 0, "Must declare number of elements >=0");
    return this;
  }

  m_type = type;

  DataType dtype = conduit::DataType::default_dtype(type);
  dtype.set_number_of_elements(num_elems);
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
  // cleanup old data
  cleanup();
  std::size_t alloc_size = getTotalBytes();
  m_data = allocateBytes(alloc_size);
  m_node.set_external(m_schema, m_data);

  return this;
}

/*
 *************************************************************************
 *
 * Declare and allocate data described by type and num elements.
 *
 *************************************************************************
 */
DataBuffer * DataBuffer::allocate(TypeID type, SidreLength num_elems)
{
  if ( num_elems < 0 )
  {
    SLIC_CHECK_MSG(num_elems >= 0, "Must allocate number of elements >=0");
    return this;
  }

  declare(type, num_elems);
  allocate();

  return this;
}

/*
 *************************************************************************
 *
 * Reallocate data to given number of elements.
 *
 *************************************************************************
 */
DataBuffer * DataBuffer::reallocate( SidreLength num_elems)
{
  if ( num_elems < 0 || m_data == ATK_NULLPTR )
  {
    SLIC_CHECK_MSG(num_elems >= 0, "Must re-allocate number of elements >=0");
    SLIC_CHECK_MSG( m_data != ATK_NULLPTR,
                     "Attempting to reallocate an unallocated buffer");
    return this;
  }

  std::size_t old_size = getTotalBytes();
  // update the buffer's Conduit Node
  DataType dtype = conduit::DataType::default_dtype(m_type);
  dtype.set_number_of_elements( num_elems );
  m_schema.set(dtype);

  std::size_t new_size = getTotalBytes();

  void * realloc_data = allocateBytes(new_size);

  std::memcpy(realloc_data, m_data, std::min(old_size, new_size));

  // cleanup old data
  cleanup();

  // let the buffer hold the new data
  m_data = realloc_data;

  // update the conduit node data pointer
  m_node.set_external(m_schema, m_data);

  return this;
}

/*
 *************************************************************************
 *
 * Update contents of buffer from src and which is nbytes long.
 *
 *************************************************************************
 */
DataBuffer * DataBuffer::update(const void * src, size_t nbytes)
{
  if ( nbytes > getTotalBytes() )
  {
    SLIC_CHECK_MSG(nbytes <= getTotalBytes(), "Unable to copy data into buffer, size exceeds available # bytes in buffer.");
    return this;
  }

  std::memcpy(m_data, src, nbytes);

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
  n.to_json_stream(os);
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
  m_type(EMPTY_ID),
  m_data(ATK_NULLPTR),
  m_node(),
  m_schema()
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
  m_type(EMPTY_ID),
  m_data(source.m_data),
  m_node(source.m_node),
  m_schema(source.m_schema)
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
  // cleanup allocated data
  if ( m_data != ATK_NULLPTR )
  {
    releaseBytes(m_data);
  }
}

/*
 *************************************************************************
 *
 * PRIVATE allocateBytes
 * Note: We allow a zero bytes allocation ( since it's legal for new() ).
 *************************************************************************
 */
void * DataBuffer::allocateBytes(std::size_t num_bytes)
{
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
  delete [] ((char *)ptr);
  m_data = ATK_NULLPTR;
}


} /* end namespace sidre */
} /* end namespace asctoolkit */

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
#include <cstring> // for std::memcpy

// Sidre project headers
#include "DataGroup.hpp"
#include "DataView.hpp"

namespace asctoolkit
{
namespace sidre
{

/*
 *************************************************************************
 *
 * Describe Buffer with given data type and number of elements.
 *
 *************************************************************************
 */
DataBuffer * DataBuffer::describe(TypeID type, SidreLength num_elems)
{
  if ( isAllocated() || num_elems < 0 )
  {
    SLIC_CHECK_MSG(!isAllocated(), "Cannot describe an allocated Buffer");
    SLIC_CHECK_MSG(num_elems >= 0, "Must describe Buffer with num elems >= 0");
    return this;
  }

  DataType& dtype = const_cast<DataType&>(m_node.dtype());
  dtype.set( dtype.default_dtype(type) );
  dtype.set_number_of_elements(num_elems);

  return this;
}

/*
 *************************************************************************
 *
 * Allocate data previously described.
 *
 *************************************************************************
 */
DataBuffer * DataBuffer::allocate()
{
  if ( !isDescribed() || isAllocated() )
  {
    SLIC_CHECK_MSG(isDescribed(), "Buffer is not described, cannot allocate.");
    SLIC_CHECK_MSG(!isAllocated(), "Buffer is already allocated.");

    return this;
  }

  void * data = allocateBytes( getTotalBytes() );

  SLIC_CHECK_MSG( data != ATK_NULLPTR,
                  "Buffer failed to allocate memory of size " <<
                  getTotalBytes() );

  if (data != ATK_NULLPTR)
  {
    m_node.set_external( DataType( m_node.dtype() ), data );
  }
  return this;
}

/*
 *************************************************************************
 *
 * Describe and allocate data described by type and num elements.
 *
 *************************************************************************
 */
DataBuffer * DataBuffer::allocate(TypeID type, SidreLength num_elems)
{
  if (isAllocated())
  {
    SLIC_CHECK_MSG(!isAllocated(), "Buffer is already allocated.");

    return this;
  }
 
  describe(type, num_elems);
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
  if (!isAllocated())
  {
    if (isDescribed())
    {
      allocate();
    } 
    else
    {
       SLIC_CHECK_MSG(isDescribed(),
                      "Can't re-allocate Buffer with no type description.");
    }

    return this;
  }

  if ( num_elems < 0 )
  {
    SLIC_CHECK_MSG(num_elems >= 0,
                   "Cannot re-allocate with number of elements < 0");
    return this;
  }

  SidreLength old_size = getTotalBytes();
  void * old_data_ptr = getVoidPtr();

  DataType dtype( m_node.dtype() );
  dtype.set_number_of_elements( num_elems );
  SidreLength new_size = dtype.total_bytes();
  void * new_data_ptr = allocateBytes(new_size);

  if ( new_data_ptr != ATK_NULLPTR )
  {
    m_node.reset();
    m_node.set_external(dtype, new_data_ptr);
    copyBytesIntoBuffer(old_data_ptr, std::min(old_size, new_size) );
    releaseBytes( old_data_ptr);
  } 
  else 
  {
     SLIC_CHECK_MSG(new_data_ptr != ATK_NULLPTR,
                    "Buffer re-allocate failed with " << new_size << " bytes.");
  }

  return this;
}

/*
 *************************************************************************
 *
 * Deallocate data in a Buffer.
 *
 *************************************************************************
 */
DataBuffer * DataBuffer::deallocate()
{
  if (!isAllocated())
  {
    return this;
  }

  releaseBytes(getVoidPtr());
  m_node.set_external( DataType( m_node.dtype() ), ATK_NULLPTR );

  for (size_t i = 0 ; i < m_views.size() ; ++i)
  {
    m_views[i]->unapply();
  }

  return this;
}

/*
 *************************************************************************
 *
 * Update contents of Buffer by copying nbytes of data into the Buffer 
 * from src.
 *
 *************************************************************************
 */
DataBuffer * DataBuffer::copyBytesIntoBuffer(const void * src, 
                                             SidreLength nbytes)
{
  if ( src == ATK_NULLPTR || nbytes < 0 || nbytes > getTotalBytes() )
  {
    SLIC_CHECK_MSG(src != ATK_NULLPTR, 
      "Cannot copy data into Buffer from null pointer.");
    SLIC_CHECK_MSG(nbytes >= 0, "Cannot copy < 0 bytes of data into Buffer.");
    SLIC_CHECK_MSG(nbytes <= getTotalBytes(),
      "Unable to copy " << nbytes << " bytes of data into Buffer with " << 
      getTotalBytes() << " bytes allocated.");

    return this;
  }

  copyBytes(src, getVoidPtr(), nbytes);

  return this;
}

/*
 *************************************************************************
 *
 * Copy data Buffer description to given Conduit node.
 *
 *************************************************************************
 */
void DataBuffer::copyToConduitNode(Node &n) const
{
  n["index"].set(m_index);
  n["node"].set(m_node.to_json());
}

/*
 *************************************************************************
 *
 * Print JSON description of data Buffer to stdout.
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
 * Print JSON description of data Buffer to an ostream.
 *
 *************************************************************************
 */
void DataBuffer::print(std::ostream& os) const
{
  Node n;
  copyToConduitNode(n);
  n.to_json_stream(os);
}



/*
 *************************************************************************
 *
 * PRIVATE ctor taking unique index.
 *
 *************************************************************************
 */
DataBuffer::DataBuffer( IndexType uid )
  : m_index(uid),
  m_views(),
  m_node()
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
  m_node(source.m_node)
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
  releaseBytes(getVoidPtr());
}

/*
 *************************************************************************
 *
 * PRIVATE method to attach Buffer to View (Buffer bookkeeping).
 *
 *************************************************************************
 */
void DataBuffer::attachToView( DataView * view )
{
  SLIC_ASSERT(view->m_data_buffer == this);

  if (view->m_data_buffer == this) 
  {
    m_views.push_back( view );
  }
}

/*
 *************************************************************************
 *
 * PRIVATE method to Buffer from View (Buffer bookkeeping).
 *
 *************************************************************************
 */
void DataBuffer::detachFromView( DataView * view )
{
  SLIC_ASSERT(view->m_data_buffer == this);

  if (view->m_data_buffer == this)
  {
    std::vector<DataView *>::iterator pos = std::find(m_views.begin(),
                                                      m_views.end(),
                                                      view);
    if ( pos != m_views.end() )
    {
      SLIC_ASSERT(pos != m_views.end());
      m_views.erase(pos);
      view->setBufferViewToEmpty();
    } 
  }
}

/*
 *************************************************************************
 *
 * PRIVATE method to detach Buffer from all Views (Buffer bookkeeping).
 *
 *************************************************************************
 */
void DataBuffer::detachFromAllViews()
{
  for (size_t i = 0 ; i < m_views.size() ; ++i)
  {
    m_views[i]->setBufferViewToEmpty();
  }

  m_views.clear();
}

/*
 *************************************************************************
 *
 * PRIVATE allocateBytes
 * Note: We allow a zero bytes allocation ( since it's legal for new() ).
 *
 *************************************************************************
 */
void * DataBuffer::allocateBytes(std::size_t num_bytes)
{
  return new(std::nothrow) detail::sidre_int8[num_bytes];
}

/*
 *************************************************************************
 *
 * PRIVATE copyBytes
 *
 *************************************************************************
 */
void DataBuffer::copyBytes( const void * src, void * dst, size_t num_bytes )
{
  std::memcpy( dst, src, num_bytes );
}

/*
 *************************************************************************
 *
 * PRIVATE releaseBytes
 *
 *************************************************************************
 */
void DataBuffer::releaseBytes( void * ptr)
{
  // Pointer type here should always match new call in allocateBytes.
  delete[] static_cast<detail::sidre_int8 *>(ptr);
}

/*
 *************************************************************************
 *
 * PRIVATE exportTo
 * Serializes Buffer state to a conduit node
 *
 *************************************************************************
 */
void DataBuffer::exportTo( conduit::Node& data_holder)
{
  data_holder["id"] = m_index;
  data_holder["schema_json"] = m_node.schema().to_json();

  // If Buffer is allocated, export it's node's data
  if ( isAllocated() )
  {
    // Do this instead of using the node copy constructor ( keep it zero-copy ).
    data_holder["data"].set_external( m_node.schema(), getVoidPtr() );

    // TODO - Ask Cyrus why he had following way previously.  Are we creating a default dtype
    // to remove any striding, offset?  Our Buffer does not allow those, only type and length.
    // DataType& dtype = conduit::DataType::default_dtype(ds_buff->getTypeID());
    // dtype.set_number_of_elements(ds_buff->getNumElements());
  }
}

/*
 *************************************************************************
 *
 * PRIVATE importFrom
 * Import Buffer state from a conduit node
 *
 *************************************************************************
 */
void DataBuffer::importFrom( conduit::Node& buffer_holder)
{
  Schema schema( buffer_holder["schema_json"].as_string() );
  TypeID type = static_cast<TypeID>( schema.dtype().id() );
  SidreLength num_elems = schema.dtype().number_of_elements();

  describe(type, num_elems);

  // If Buffer was allocated, the conduit node will have the entry "data".
  // Allocate and copy in that data.
  if (buffer_holder.has_path("data"))
  {
    allocate();
    conduit::Node& buffer_data_holder = buffer_holder["data"];
    copyBytesIntoBuffer(buffer_data_holder.element_ptr(0),
           buffer_data_holder.total_bytes() );
  }
  else
  {
    std::cerr << "NO PATH??" << std::endl;
  }
}


} /* end namespace sidre */
} /* end namespace asctoolkit */

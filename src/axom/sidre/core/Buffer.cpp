/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

// Associated header file
#include "Buffer.hpp"

// Standard C++ headers
#include <algorithm>
#include <cstring> // for std::memcpy

// Sidre project headers
#include "Group.hpp"
#include "View.hpp"

#include "axom/core/Types.hpp" // for common::int8

namespace axom
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
Buffer* Buffer::describe(TypeID type, IndexType num_elems)
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
Buffer* Buffer::allocate()
{
  if ( !isDescribed() || isAllocated() )
  {
    SLIC_CHECK_MSG(isDescribed(), "Buffer is not described, cannot allocate.");
    SLIC_CHECK_MSG(!isAllocated(), "Buffer is already allocated.");

    return this;
  }

  void* data = allocateBytes( getTotalBytes() );

  SLIC_CHECK_MSG( data != nullptr,
                  "Buffer failed to allocate memory of size " <<
                  getTotalBytes() );

  if (data != nullptr)
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
Buffer* Buffer::allocate(TypeID type, IndexType num_elems)
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
Buffer* Buffer::reallocate( IndexType num_elems)
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

  IndexType old_size = getTotalBytes();
  void* old_data_ptr = getVoidPtr();

  DataType dtype( m_node.dtype() );
  dtype.set_number_of_elements( num_elems );
  IndexType new_size = dtype.strided_bytes();
  void* new_data_ptr = allocateBytes(new_size);

  if ( new_data_ptr != nullptr )
  {
    m_node.reset();
    m_node.set_external(dtype, new_data_ptr);
    copyBytesIntoBuffer(old_data_ptr, std::min(old_size, new_size) );
    releaseBytes( old_data_ptr);
  }
  else
  {
    SLIC_CHECK_MSG(new_data_ptr != nullptr,
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
Buffer* Buffer::deallocate()
{
  if (!isAllocated())
  {
    return this;
  }

  releaseBytes(getVoidPtr());
  m_node.set_external( DataType( m_node.dtype() ), nullptr );

  std::set<View*>::iterator vit = m_views.begin();
  for ( ; vit != m_views.end() ; ++vit)
  {
    (*vit)->unapply();
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
Buffer* Buffer::copyBytesIntoBuffer(const void* src,
                                    IndexType nbytes)
{
  if ( src == nullptr || nbytes < 0 || nbytes > getTotalBytes() )
  {
    SLIC_CHECK_MSG(src != nullptr,
                   "Cannot copy data into Buffer from null pointer.");
    SLIC_CHECK_MSG(nbytes >= 0, "Cannot copy < 0 bytes of data into Buffer.");
    SLIC_CHECK_MSG(
      nbytes <= getTotalBytes(),
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
void Buffer::copyToConduitNode(Node &n) const
{
  n["index"].set(m_index);
  n["value"].set(m_node.to_json());
}

/*
 *************************************************************************
 *
 * Print JSON description of data Buffer to stdout.
 *
 *************************************************************************
 */
void Buffer::print() const
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
void Buffer::print(std::ostream& os) const
{
  Node n;
  copyToConduitNode(n);
  n.to_json_stream(os);
}


/*
 *************************************************************************
 *
 * Serializes Buffer state to a conduit node
 *
 *************************************************************************
 */
void Buffer::exportTo( conduit::Node& data_holder)
{
  data_holder["id"] = m_index;

  if ( isDescribed() )
  {
    data_holder["schema"] = m_node.schema().to_json();
  }

  // If Buffer is allocated, export it's node's data
  if ( isAllocated() )
  {
    // Do this instead of using the node copy constructor ( keep it zero-copy ).
    data_holder["data"].set_external( m_node.schema(), getVoidPtr() );

    // TODO - Ask Cyrus why he had following way previously.  Are we creating a
    // default dtype
    // to remove any striding, offset?  Our Buffer does not allow those, only
    // type and length.
    // DataType& dtype = conduit::DataType::default_dtype(ds_buff->getTypeID());
    // dtype.set_number_of_elements(ds_buff->getNumElements());
  }
}


/*
 *************************************************************************
 *
 * Import Buffer state from a conduit node
 *
 *************************************************************************
 */
void Buffer::importFrom( conduit::Node& buffer_holder)
{
  if (buffer_holder.has_path("schema"))
  {
    Schema schema( buffer_holder["schema"].as_string() );
    TypeID type = static_cast<TypeID>( schema.dtype().id() );
    IndexType num_elems = schema.dtype().number_of_elements();
    describe(type, num_elems);
  }

  // If Buffer was allocated, the conduit node will have the entry "data".
  // Allocate and copy in that data.
  if (buffer_holder.has_path("data"))
  {
    allocate();
    conduit::Node& buffer_data_holder = buffer_holder["data"];
    copyBytesIntoBuffer(buffer_data_holder.element_ptr(0),
                        buffer_data_holder.total_strided_bytes() );
  }
}


/*
 *************************************************************************
 *
 * PRIVATE ctor taking unique index.
 *
 *************************************************************************
 */
Buffer::Buffer( IndexType uid )
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
Buffer::Buffer(const Buffer& source )
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
Buffer::~Buffer()
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
void Buffer::attachToView( View* view )
{
  SLIC_ASSERT(view->m_data_buffer == this);

  if (view->m_data_buffer == this)
  {
    m_views.insert( view );
  }
}

/*
 *************************************************************************
 *
 * PRIVATE method to Buffer from View (Buffer bookkeeping).
 *
 *************************************************************************
 */
void Buffer::detachFromView( View* view )
{
  SLIC_ASSERT(view->m_data_buffer == this);
  SLIC_ASSERT(m_views.count(view) > 0);

  if (view->m_data_buffer == this && m_views.count(view) > 0)
  {
    m_views.erase(view);
    view->setBufferViewToEmpty();
  }
}

/*
 *************************************************************************
 *
 * PRIVATE method to detach Buffer from all Views (Buffer bookkeeping).
 *
 *************************************************************************
 */
void Buffer::detachFromAllViews()
{
  std::set<View*>::iterator vit = m_views.begin();
  for ( ; vit != m_views.end() ; ++vit)
  {
    (*vit)->setBufferViewToEmpty();
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
void* Buffer::allocateBytes(IndexType num_bytes)
{
  return new(std::nothrow) common::int8[num_bytes];
}

/*
 *************************************************************************
 *
 * PRIVATE copyBytes
 *
 *************************************************************************
 */
void Buffer::copyBytes( const void* src, void* dst, IndexType num_bytes )
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
void Buffer::releaseBytes( void* ptr)
{
  // Pointer type here should always match new call in allocateBytes.
  delete[] static_cast<common::int8*>(ptr);
}


} /* end namespace sidre */
} /* end namespace axom */

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
    DataType dtype = conduit::DataType::default_dtype(type);
    dtype.set_number_of_elements(len);
    m_schema.set(dtype);
  }
  return this;
}

/*
 *************************************************************************
 *
 * Declare buffer to OWN data described as a Conduit schema.
 *
 *************************************************************************
 */
DataBuffer * DataBuffer::declare(const Schema& schema)
{
  m_schema.set(schema);
  return this;
}

/*
 *************************************************************************
 *
 * Declare buffer to OWN data described as a Conduit data type.
 *
 *************************************************************************
 */
DataBuffer * DataBuffer::declare(const DataType& dtype)
{
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
  SLIC_ASSERT_MSG( !m_is_data_external,
                  "Attempting to allocate buffer holding external data");

  if(m_fortran_allocatable != ATK_NULLPTR)
  {
    // cleanup old data
    cleanup();
    TypeID type = static_cast<TypeID>(m_schema.dtype().id());
    SidreLength nitems = m_schema.dtype().number_of_elements();
    m_data = AllocateAllocatable(m_fortran_allocatable, type, m_rank, nitems);
    m_node.set_external(m_schema,m_data);
  }
  else if ( !m_is_data_external )
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
  SLIC_ASSERT_MSG( !m_is_data_external,
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
    Schema s(dtype);
    reallocate(s);
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
 * Set as Fortran allocatable buffer.
 *
 *************************************************************************
 */
DataBuffer * DataBuffer::setFortranAllocatable(void * array, TypeID type, int rank)
{
  SLIC_ASSERT_MSG( array != ATK_NULLPTR, 
		   "Attempting to set buffer to Fortran allocatable given null pointer" );
  // XXX check rank too

  if ( array != ATK_NULLPTR )
  {
#ifdef ATK_ENABLE_FORTRAN
    m_fortran_allocatable = array;
    m_rank = rank;
    m_data = AddressAllocatable(array, type, rank);
    m_node.set_external(m_schema, m_data);
#else
    SLIC_ASSERT_MSG( true ,
		     "Fortran support is not compiled into this version of Sidre" );
#endif
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
  m_rank(0),
  m_data(ATK_NULLPTR),
  m_node(),
  m_schema(),
  m_is_data_external(false),
  m_fortran_allocatable(ATK_NULLPTR)
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
  m_rank(0),
  m_data(source.m_data),
  m_node(source.m_node),
  m_schema(source.m_schema),
  m_is_data_external(source.m_is_data_external),
  m_fortran_allocatable(ATK_NULLPTR)
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
    if (m_fortran_allocatable != ATK_NULLPTR)
    {
#ifdef ATK_ENABLE_FORTRAN
      TypeID type = static_cast<TypeID>(m_node.dtype().id());
      DeallocateAllocatable(m_fortran_allocatable, type, m_rank);
#else
    SLIC_ASSERT_MSG( true ,
		     "Fortran support is not compiled into this version of Sidre" );
#endif
    }
    else if (!m_is_data_external )
    {
      releaseBytes(m_data);
    }
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

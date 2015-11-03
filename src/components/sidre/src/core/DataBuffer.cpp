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
 * Return number of bytes associated with a single item of DataBuffer's type.
 *
 *************************************************************************
 */
size_t DataBuffer::getBytesPerItem() const
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

  return bytes_per_item[m_type];
}

/*
 *************************************************************************
 *
 * Return total number of bytes associated with this DataBuffer object.
 *
 *************************************************************************
 */
size_t DataBuffer::getTotalBytes() const
{
  return getBytesPerItem() * m_nitems;
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

  if(m_fortran_allocatable != ATK_NULLPTR)
  {
#ifdef ATK_ENABLE_FORTRAN
    // cleanup old data
    cleanup();
    m_data = AllocateAllocatable(m_fortran_allocatable, m_type, m_fortran_rank, m_nitems);
#else
    SLIC_ERROR("Fortran support is not compiled into this version of Sidre");
#endif
  }
  else if ( !m_is_data_external )
  {
    // cleanup old data
    cleanup();
    std::size_t alloc_size = getTotalBytes();
    m_data = allocateBytes(alloc_size);
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
 * Reallocate data using a length.
 *
 *************************************************************************
 */
DataBuffer * DataBuffer::reallocate( SidreLength len)
{
  SLIC_ASSERT_MSG(len >= 0, "Must re-allocate number of elements >=0");
  SLIC_ASSERT_MSG( !m_is_data_external,
                  "Attempting to re-allocate buffer holding external data");

  if ( !m_is_data_external )
  {
    //  make sure realloc actually makes sense
    SLIC_ASSERT_MSG(m_data != ATK_NULLPTR,
                   "Attempting to reallocate an unallocated buffer");

    std::size_t realloc_size = len * getBytesPerItem();
    void * realloc_data = allocateBytes(realloc_size);

    memcpy(realloc_data, m_data, std::min(getTotalBytes(), realloc_size));

    // cleanup old data
    cleanup();

    // let the buffer hold the new data
    m_data = realloc_data;
    m_nitems = len;
  }

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
  size_t buff_nbytes = getTotalBytes();
  SLIC_ASSERT_MSG(nbytes <= buff_nbytes, "Must allocate number of elements >=0");

  if ( src != ATK_NULLPTR && nbytes <= buff_nbytes)
  {
    memcpy(m_data, src, nbytes);
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
    m_fortran_rank = rank;
    m_data = AddressAllocatable(array, type, rank);
#else
    SLIC_ERROR("Fortran support is not compiled into this version of Sidre");
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
  // Create a conduit node
  DataType dtype = conduit::DataType::default_dtype(m_type);
  dtype.set_number_of_elements(m_nitems);
  Schema schema(dtype);
  Node node;
  node.set_external(schema, m_data);

  n["index"].set(m_index);
  n["is_data_external"].set(m_is_data_external);
  n["schema"].set(schema.to_json());
  n["node"].set(node.to_json());
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
  m_is_data_external(false),
  m_fortran_rank(0),
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
  m_type(CONDUIT_EMPTY_T),
  m_nitems(0),
  m_data(source.m_data),
  m_is_data_external(source.m_is_data_external),
  m_fortran_rank(0),
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
      DeallocateAllocatable(m_fortran_allocatable, m_type, m_fortran_rank);
#else
      SLIC_ERROR("Fortran support is not compiled into this version of Sidre");
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

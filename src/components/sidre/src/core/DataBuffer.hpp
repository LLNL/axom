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
 * \brief   Header file containing definition of DataBuffer class.
 *
 ******************************************************************************
 */

#ifndef DATABUFFER_HPP_
#define DATABUFFER_HPP_

// Standard C++ headers
#include <vector>

// Other toolkit component headers
#include "common/CommonTypes.hpp"
#include "slic/slic.hpp"

// SiDRe project headers
#include "SidreTypes.hpp"

namespace asctoolkit
{
namespace sidre
{
// using directives to make Conduit usage easier and less visible
using conduit::Node;
using conduit::Schema;

class DataStore;
class DataView;

/*!
 * \class DataBuffer
 *
 * \brief DataBuffer holds a data object, which it owns (and allocates!)
 *
 * The DataBuffer class has the following properties:
 *
 *    - DataBuffer objects can only be created via the DataStore interface,
 *      not directly.
 *    - A DataBuffer object has a unique identifier within a DataStore,
 *      which is assigned by the DataStore when the buffer is created.
 *    - The data object owned by a DataBuffer is unique to that DataBuffer
 *      object; i.e.,  DataBuffers that own data do not share their data.
 *    - Typical usage is to describe the data a DataBuffer will hold and then
 *      allocate it by calling one of the DataBuffer allocate or
 *      reallocate methods.
 *    - A DataBuffer object maintains a collection of DataViews that
 *      refer to its data.
 *
 */
class DataBuffer
{

public:

  /*!
   * Friend declarations to constrain usage via controlled access to
   * private members.
   */
  friend class DataStore;
  friend class DataGroup;
  friend class DataView;

//@{
//!  @name Accessor methods

  /*!
   * \brief Return the unique index of this buffer object.
   */
  IndexType getIndex() const
  {
    return m_index;
  }

  /*!
   * \brief Return number of views attached to this buffer.
   */
  IndexType getNumViews() const
  {
    return m_views.size();
  }

  /*!
   * \brief Return void-pointer to data held by DataBuffer.
   */
  void * getVoidPtr()
  {
    return m_node.data_ptr();
  }

  /*!
   * \brief Returns data held by node (or pointer to data if array).
   */
  Node::Value getData()
  {
    if ( !isAllocated() )
    {
      SLIC_CHECK_MSG( isAllocated(), "Buffer data is not allocated.");
      return Node().value();
    }

    return m_node.value();
  }

  /*!
   * \brief Return type of data for this DataBuffer object.
   */
  TypeID getTypeID() const
  {
    return static_cast<TypeID>(m_node.dtype().id());
  }

  /*!
   * \brief Return total number of elements allocated by this DataBuffer object.
   */
  SidreLength getNumElements() const
  {
    return m_node.dtype().number_of_elements();
  }

  /*!
   * \brief Return total number of bytes associated with this DataBuffer object.
   */
  SidreLength getTotalBytes() const
  {
    return m_node.dtype().total_bytes();
  }
  //@}

  /*!
   * \brief Return true if buffer contains allocated data of > 0 bytes.
   */
  bool isAllocated() const
  {
    return (m_node.data_ptr() != ATK_NULLPTR) && (getTotalBytes() > 0);
  }

  /*!
   * \brief Return true if data description exists.
   */
  bool isDescribed() const
  {
    return !m_node.dtype().is_empty();
  }

  /*
   * \brief Return true if DataBuffer has an associated DataView with given
   *        index; else false.
   */
  bool hasView( IndexType idx ) const
  {
    return ( 0 <= idx && static_cast<unsigned>(idx) < m_views.size() &&
             m_views[idx] != ATK_NULLPTR );
  }

  /*!
   * \brief Return (non-const) pointer to data view object with given index
   *        associated with buffer, or ATK_NULLPTR if none exists.
   */
  DataView * getView( IndexType idx);

//@}

//@{
//!  @name Data declaration and allocation methods

  /*!
   * \brief Describe a buffer with data given type and number of elements.
   *
   * To use the buffer, the data must be allocated by calling allocate().
   *
   * If given number of elements is < 0, method does nothing.
   *
   * \return pointer to this DataBuffer object.
   */
  DataBuffer * describe(TypeID type, SidreLength num_elems);

  /*!
   * \brief Allocate data previously described using a describe() method.
   *
   * It is the responsibility of the caller to make sure that the buffer
   * object was previously described.  If the the buffer is already
   * holding data that it owns, that data will be deallocated and new data
   * will be allocated according to the current described state.
   *
   * \return pointer to this DataBuffer object.
   */
  DataBuffer * allocate();

  /*!
   * \brief Describe and allocate data described by type and number of elements.
   *
   * This is equivalent to calling describe(type, num_elems), then allocate().
   * on this DataBuffer object.
   *
   * \return pointer to this DataBuffer object.
   */
  DataBuffer * allocate(TypeID type, SidreLength num_elems);

  /*!
   * \brief Reallocate data to given number of elements.
   *
   *        Equivalent to calling describe(type), then allocate().
   *
   * If given number of elements < 0, this method does nothing.
   *
   * \return pointer to this DataBuffer object.
   */
  DataBuffer * reallocate(SidreLength num_elems);

  /*!
   * \brief Deallocate data in a buffer.
   *
   * If the buffer has no data, the routine does nothing.
   * All attached views will continue to be attached;
   * however, their data address will be set to ATK_NULLPTR
   * with a data length of 0.
   *
   * \return pointer to this DataBuffer object.
   */
  DataBuffer * deallocate();

  /*!
   * \brief Update contents of buffer memory.
   *
   * This will copy nbytes of data into the buffer.  nbytes must be greater
   * than 0 and less than getTotalBytes().
   *
   * If given pointer is null, this method does nothing.
   *
   * \return pointer to this DataBuffer object.
   */
  DataBuffer * update(const void * src, SidreLength nbytes);

  //@}

  /*!
   * \brief Copy data buffer description to given Conduit node.
   */
  void info(Node& n) const;

  /*!
   * \brief Print JSON description of data buffer to stdout.
   */
  void print() const;

  /*!
   * \brief Print JSON description of data buffer to an ostream.
   */
  void print(std::ostream& os) const;


private:

  /*!
   *  \brief Private ctor that assigns unique id.
   */
  DataBuffer( IndexType uid );

  /*!
   * \brief Private copy ctor.
   */
  DataBuffer(const DataBuffer& source );

  /*!
   * \brief Private dtor.
   */
  ~DataBuffer();

  /*!
   * \brief Private methods to attach/detach data view to buffer.
   */
  void attachView( DataView * dataView );
  ///
  void detachView( DataView * dataView );

  /*!
   * \brief Private method to detach all views to the buffer.
   */
  void detachAllViews( );

  /*!
   * \brief Allocate bytes for data in data buffer
   */
  void * allocateBytes(std::size_t num_bytes);

  /*!
   * \brief Copy bytes from one memory location to another.
   */
  void copyBytes( const void * src, void * dst, size_t num_bytes );

  /*!
   * \brief Release any allocated bytes pointed to by ptr.
   */
  void  releaseBytes(void * ptr);

  /*!
   * \brief Exports buffer's state to a conduit node.
   */
  void exportTo( conduit::Node& data_holder );

  /*!
   * \brief Import buffer's state from a conduit node.
   */
  void importFrom( conduit::Node& data_holder );

  /// Index Identifier - unique within a dataStore.
  IndexType m_index;

  /// Container of DataViews attached to this buffer.
  std::vector<DataView *> m_views;

  /// Conduit Node that holds buffer data.
  Node m_node;

  /*!
   *  Unimplemented ctors and copy-assignment operators.
   */
#ifdef USE_CXX11
  DataBuffer() = delete;
  DataBuffer( DataBuffer&& ) = delete;

  DataBuffer& operator=( const DataBuffer& ) = delete;
  DataBuffer& operator=( DataBuffer&& ) = delete;
#else
  DataBuffer();
  DataBuffer& operator=( const DataBuffer& );
#endif

};

} /* end namespace sidre */
} /* end namespace asctoolkit */

#endif /* DATABUFFER_HPP_ */

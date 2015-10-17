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
#include <map>
#include <set>
#include <vector>

// Other toolkit component headers
#include "common/CommonTypes.hpp"

// SiDRe project headers
#include "SidreTypes.hpp"
#include "MetaBuffer.hpp"
#include "sidre/SidreAllocatable.hpp"


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
 *    - The data object owned by a DataBuffer is unique to that DataDuffer
 *      object; i.e.,  DataBuffers that own data not share their data.
 *    - A DataBuffer may hold a pointer to externally-owned data. When this
 *      is the case, the buffer cannot be used to (re)allocate or deallocate
 *      the data. However, the external data can be desribed and accessed
 *      via the buffer object in a similar to data that is owned by a buffer.
 *    - Typical usage is to declare the data a DataBuffer will hold and then
 *      either allocate it by calling one of the DataBuffer allocate or
 *      reallocate methods, or set the buffer to reference externally-owned
 *      data by calling setExternalData().
 *    - A DataBuffer object maintains a collection of DataViews that
 *      refer to its data.
 *
 */
class DataBuffer
{
public:

  //
  // Friend declarations to constrain usage via controlled access to
  // private members.
  //
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
  size_t getNumViews() const
  {
    return m_views.size();
  }

  /*!
   * \brief Return true if buffer holds externally-owned data, or
   * false if buffer owns the data it holds (default case).
   */
  bool isExternal() const
  {
    return m_is_data_external;
  }

  /*!
   * \brief Return void-pointer to data held by DataBuffer.
   */
  void * getData()
  {
    return m_data;
  }

  /*!
   * \brief Return total number of bytes associated with this DataBuffer object.
   */
  size_t getTotalBytes() const
  {
    return m_schema.total_bytes();
  }

  /*!
   * \brief Return non-const reference to Conduit node holding data.
   */
  Node& getNode()
  {
    return m_node;
  }

  /*!
   * \brief Return const reference to Conduit node holding data.
   */
  const Node& getNode() const
  {
    return m_node;
  }

  /*!
   * \brief Returns Value class instance that supports casting to the appropriate data return type.  This function
   * version does require enough type information for the compiler to know what to cast the Value class to.
   * Example:
   * int* myptr = getValue();
   * int myint = getValue();
   */
  Node::Value getValue()
  {
    return m_node.value();
  }

  /*!
   * \brief Set value in conduit node.
   */
  template<typename ValueType>
  void setValue(ValueType value)
  {
    m_node.set(value);
  }


  /*!
   * \brief Lightweight templated wrapper around getValue that returns a Value class.  This function can be used in cases
   * were not enough information is provided to the compiler to cast the Value class based on the caller code line.  The
   * function template type must be explicitly provided on call.
   *
   * Example:
   * // will not work, compiler does not know what type to cast to for above getValue function.
   * assert( getValue() == 10 );
   * // use the templated version instead
   * assert (getValue<int>() == 10);
   */
  template<typename ValueType>
  ValueType getValue()
  {
    ValueType valueptr = m_node.value();
    return valueptr;
  }



  /*!
   * \brief Return const reference to Conduit schema describing data.
   */
  const Schema& getSchema() const
  {
    return m_schema;
  }

  /*!
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
   * \brief Declare a buffer with data given type and number of elements.
   *
   * To use the buffer, the data must be allocated by calling allocate()
   * or set to external data by calling setExternalData().
   *
   * If given length is < 0, method does nothing.
   *
   * \return pointer to this DataBuffer object.
   */
  DataBuffer * declare(TypeID type, SidreLength len);

  /*!
   * \brief Declare a buffer with data described as a Conduit schema.
   *
   * To use the buffer, the data must be allocated by calling allocate()
   * or set to external data by calling setExternalData().
   *
   * \return pointer to this DataBuffer object.
   */
  DataBuffer * declare(const Schema& schema);

  /*!
   * \brief Declare a buffer with data described as a pre-defined
   *        Conduit data type.
   *
   * To use the buffer, the data must be allocated by calling allocate()
   * or set to external data by calling setExternalData().
   *
   * \return pointer to this DataBuffer object.
   */
  DataBuffer * declare(const DataType& dtype);

  /*!
   * \brief Allocate data previously declared using a declare() method.
   *
   * It is the responsibility of the caller to make sure that the buffer
   * object was previously declared.  If the the buffer is already
   * holding data that it owns, that data will be deallocated and new data
   * will be allocated according to the current declared state.
   *
   * If buffer is already set to externally-owned data, this method
   * does nothing.
   *
   * \return pointer to this DataBuffer object.
   */
  DataBuffer * allocate();

  /*!
   * \brief Declare and allocate data described with type and length.
   *
   * This is equivalent to calling declare(type, len), then allocate().
   * on this DataBuffer object.
   *
   * If buffer is already set to externally-owned data, this method
   * does nothing.
   *
   * \return pointer to this DataBuffer object.
   */
  DataBuffer * allocate(TypeID type, SidreLength len);

  /*!
   * \brief Declare and allocate data described as a Conduit schema.
   *
   * This is equivalent to calling declare(schema), then allocate().
   * on this DataBuffer object.
   *
   * If buffer is already set to externally-owned data, this method
   * does nothing.
   *
   * \return pointer to this DataBuffer object.
   */
  DataBuffer * allocate(const Schema &schema);

  /*!
   * \brief Declare and allocate data described as a pre-defined
   *        Conduit data type.
   *
   * This is equivalent to calling declare(dtype), then allocate().
   * on this DataBuffer object.
   *
   * If buffer is already set to externally-owned data, this method
   * does nothing.
   *
   * \return pointer to this DataBuffer object.
   */
  DataBuffer * allocate(const DataType& dtype);

  /*!
   * \brief Reallocate data described with Sidre type and length.
   *
   *        Equivalent to calling declare(type), then allocate().
   *
   * If buffer is already set to externally-owned data or given length < 0,
   * this method does nothing.
   *
   * \return pointer to this DataBuffer object.
   */
  DataBuffer * reallocate(TypeID type, SidreLength len);

  /*!
   * \brief Reallocate data described as a Conduit schema.
   *
   * If buffer is already set to externally-owned data, this method
   * does nothing.
   *
   * \return pointer to this DataBuffer object.
   */
  DataBuffer * reallocate(const Schema& schema);

  /*!
   * \brief Reallocate data described as a pre-defined
   *        Conduit data type.
   *
   * If buffer is already set to externally-owned data, this method
   * does nothing.
   *
   * \return pointer to this DataBuffer object.
   */
  DataBuffer * reallocate(const DataType& dtype);

  /*!
   * \brief Set buffer to external data.
   *
   * It is the responsibility of the caller to make sure that the buffer
   * object was previously declared, that the data pointer is consistent
   * with how the buffer was declared, and that the buffer is not already
   * holding data that it owns.
   *
   * If given pointer is null, this method does nothing.
   *
   * \return pointer to this DataBuffer object.
   */
  DataBuffer * setExternalData(void * external_data);

  /*!
   * \brief Set as Fortran allocatable.
   *
   * If given pointer is null, this method does nothing.
   *
   * \return pointer to this DataBuffer object.
   */
  DataBuffer * setFortranAllocatable(void * array, TypeID type, int rank);

  /*!
   * \brief Set buffer to external data.
   *
   * \return pointer to this DataBuffer object.
   */
  DataBuffer * setMetaBuffer(MetaBuffer * meta_buffer)
  {
    m_meta_buffer = meta_buffer;
    return this;
  };


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
   * \brief Private methods allocate and deallocate data owned by buffer.
   */
  void cleanup();
  ///
  void * allocateBytes(std::size_t num_bytes);
  ///
  void  releaseBytes(void * );


  /// Index Identifier - unique within a dataStore.
  IndexType m_index;

  /// Container of DataViews attached to this buffer.
  std::vector<DataView *> m_views;

  /// Number of dimensions
  int m_rank;

  /// Pointer to the data owned by DataBuffer.
  void * m_data;

  /// Conduit Node that holds buffer data.
  Node m_node;

  /// Conduit Schema that describes buffer data.
  Schema m_schema;

  /// Is buffer holding externally-owned data?
  bool m_is_data_external;

  /// Is buffer holding Fortran allocatable?
  void *m_fortran_allocatable;

  MetaBuffer * m_meta_buffer;
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

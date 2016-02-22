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
 * \brief   Header file containing definition of DataStore class.
 *
 ******************************************************************************
 */

#ifndef DATASTORE_HPP_
#define DATASTORE_HPP_

// Standard C++ headers
#include <vector>
#include <stack>

// Other CS Toolkit headers
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

class DataBuffer;
class DataGroup;

/*!
 * \class DataStore
 *
 * \brief DataStore is the main interface for creating and accessing
 *        DataBuffer objects.
 *
 * It maintains a collection of DataBuffer objects and owns the "root"
 * DataGroup, called "/". A DataGroup hierachy (a tree) is created by
 * creating (sub) DataGroups withing the root group and (sub) DataGroups
 * within child DataGroups.
 *
 */
class DataStore
{
public:

  /*!
   * \brief Default ctor initializes datastore object and creates root group.
   */
  DataStore();

  /*!
   * \brief Dtor destroys all contents of the datastore, including data held
   *        in buffers if owned by the buffers.
   */
  ~DataStore();

  /*!
   * \brief Return pointer to the root DataGroup.
   */
  DataGroup * getRoot()
  {
    return m_RootGroup;
  };


//@{
//!  @name DataBuffer methods

  /*!
   * \brief Return true if DataStore owns a DataBuffer with given index;
   *        else false.
   */
  bool hasBuffer( IndexType idx ) const
  {
    return ( 0 <= idx && static_cast<unsigned>(idx) < m_data_buffers.size() &&
             m_data_buffers[idx] != ATK_NULLPTR );
  }

  /*!
   * \brief Return (non-const) pointer to data buffer object with given index,
   *        or ATK_NULLPTR if none exists.
   */
  DataBuffer * getBuffer( IndexType idx ) const;

  /*!
   * \brief Create an undescribed data buffer object and return a pointer to it.
   *    When allocating, the type and number of elements must be provided
   *    in the allocate call.
   *
   *    The buffer object is assigned a unique index when created and the
   *    buffer object is owned by the data store.
   */
  DataBuffer * createBuffer();

  /*!
   * \brief Create a data buffer object with specified type and number of
   *    elements and return a pointer to it.
   *
   *    The buffer object is assigned a unique index when created and the
   *    buffer object is owned by the data store.
   */
  DataBuffer * createBuffer( TypeID type, SidreLength num_elems );

  /*!
   * \brief Remove data buffer with given index from the datastore and
   *        destroy it (including its data if data buffer owns it).
   *
   *   Note that buffer destruction detaches it from all groups and views
   *   it was associated with.
   */
  void destroyBuffer( IndexType idx );

  /*!
   * \brief Remove all data buffers from the datastore and destroy them
   *        (including data they own).
   *
   *   Note that buffer destruction detaches it from all groups and views
   *   it was associated with.
   */
  void destroyBuffers();

  /*!
   * \brief Remove data buffer with given index from the datastore, but leave
   *        it intact.
   *
   * \return pointer to DataBuffer object that was datached.
   */
  DataBuffer * detachBuffer( IndexType idx );

  /*!
   * \brief Return number of buffers in the datastore.
   */
  size_t getNumBuffers() const
  {
    return m_data_buffers.size() - m_free_buffer_ids.size();
  }

//@}

//@{
//!  @name DataBuffer iteration methods

  /*!
   * \brief Return first valid DataBuffer index (i.e., smallest index
   *        over all DataBuffers).
   *
   * sidre::InvalidIndex is returned if group has no buffers.
   */
  IndexType getFirstValidBufferIndex() const;

  /*!
   * \brief Return next valid DataBuffer index after given index (i.e.,
   *        smallest index over all buffer indices larger than given one).
   *
   * sidre::InvalidIndex is returned if there is no valid index greater
   * than given one.
   */
  IndexType getNextValidBufferIndex(IndexType idx) const;

//@}

  /*!
   * \brief Copy buffer descriptions and group tree, starting at root,
   *        to given Conduit node.
   */
  void info(Node& n) const;

  /*!
   * \brief Print JSON description of data buffers and group tree,
   *        starting at root, to stdout.
   */
  void print() const;

  /*!
   * \brief Print JSON description of data buffers and group tree,
   *        starting at root, to an ostream.
   */
  void print(std::ostream& os) const;


private:
  /*!
   *  Unimplemented ctors and copy-assignment operators.
   */
#ifdef USE_CXX11
  DataStore( const DataStore& ) = delete;
  DataStore( DataStore&& ) = delete;

  DataStore& operator=( const DataStore& ) = delete;
  DataStore& operator=( DataStore&& ) = delete;
#else
  DataStore( const DataStore& );
  DataStore& operator=( const DataStore& );
#endif

  /// Root data group, created when DataStore instance is created.
  DataGroup * m_RootGroup;

  /// Collection of DataBuffers holding data in DataStore instance.
  std::vector<DataBuffer *> m_data_buffers;

  /// Collection of unused unique buffer indices (they can be recycled).
  std::stack< IndexType > m_free_buffer_ids;

  /// Flag indicating whether SLIC logging environment was initialized in ctor.
  bool m_i_initialized_slic;
};


} /* end namespace sidre */
} /* end namespace asctoolkit */

#endif /* DATASTORE_HPP_ */

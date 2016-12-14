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

#include "hdf5.h"

// Other CS Toolkit headers
#include "common/config.hpp"    // defines ATK_USE_CXX11
#include "common/CommonTypes.hpp"
#include "slic/slic.hpp"

// Sidre project headers
#include "SidreTypes.hpp"

namespace asctoolkit
{
namespace sidre
{

class DataBuffer;
class DataGroup;

/*!
 * \class DataStore
 *
 * \brief DataStore is the main interface for creating and accessing
 *        Buffer objects.
 *
 * It maintains a collection of Buffer objects and owns the "root"
 * Group, called "/". A Group hierarchy (a tree) is created by
 * creating child Groups of Groups.
 */
class DataStore
{
public:

  /*!
   * \brief Default ctor initializes DataStore object and creates root Group.
   *
   * Ctor also initializes SLIC logging environment if it is not already
   * initialized.
   */
  DataStore();

  /*!
   * \brief Dtor destroys all contents of the DataStore, including data held
   *        in Buffers.
   */
  ~DataStore();

  /*!
   * \brief Return pointer to the root Group.
   */
  DataGroup * getRoot()
  {
    return m_RootGroup;
  };

  /*!
   * \brief Return pointer to the root DataGroup.
   */
  const DataGroup * getRoot() const
  {
    return m_RootGroup;
  };

//@{
//!  @name Methods to query, access, create, and destroy Buffers.

  /*!
   * \brief Return number of Buffers in the DataStore.
   */
  size_t getNumBuffers() const
  {
    return m_data_buffers.size() - m_free_buffer_ids.size();
  }

  /*!
   * \brief Return true if DataStore owns a Buffer with given index;
   *        else false.
   */
  bool hasBuffer( IndexType idx ) const
  {
    return ( 0 <= idx && static_cast<unsigned>(idx) < m_data_buffers.size() &&
             m_data_buffers[static_cast<unsigned>(idx)] != ATK_NULLPTR );
  }

  /*!
   * \brief Return (non-const) pointer to Buffer object with given index,
   *        or ATK_NULLPTR if none exists.
   */
  DataBuffer * getBuffer( IndexType idx ) const;

  /*!
   * \brief Create an undescribed Buffer object and return a pointer to it.
   *
   *        The Buffer must be described before it can be allocated.
   *
   *        The Buffer object is assigned a unique index when created and the
   *        Buffer object is owned by the DataStore object.
   */
  DataBuffer * createBuffer();

  /*!
   * \brief Create a Buffer object with specified type and number of
   *        elements and return a pointer to it.
   *
   *        See the DataBuffer::describe() method for valid data description.
   *
   *        The Buffer object is assigned a unique index when created and the
   *        Buffer object is owned by the DataStore object.
   */
  DataBuffer * createBuffer( TypeID type, SidreLength num_elems );

  /*!
   * \brief Remove Buffer from the DataStore and destroy it and
   *        its data.
   *
   *        Note that Buffer destruction detaches it from all Views to
   *        which it is attached.
   */
  void destroyBuffer( DataBuffer * buff );

  /*!
   * \brief Remove Buffer with given index from the DataStore and
   *        destroy it and its data.
   *
   *        Note that Buffer destruction detaches it from all Views to
   *        which it is attached.
   */
  void destroyBuffer( IndexType idx );

  /*!
   * \brief Remove all Buffers from the DataStore and destroy them
   *        and their data.
   *
   *        Note that Buffer destruction detaches it from all Views to
   *        which it is attached.
   */
  void destroyAllBuffers();

//@}

//@{
//!  @name Methods for iterating over Buffers in DataStore

  /*!
   * \brief Return first valid Buffer index.
   *
   *        sidre::InvalidIndex is returned if Group has no Buffers.
   */
  IndexType getFirstValidBufferIndex() const;

  /*!
   * \brief Return next valid Buffer index after given index.
   *
   *        sidre::InvalidIndex is returned if there is no valid next index.
   */
  IndexType getNextValidBufferIndex(IndexType idx) const;

//@}

  /*!
   * \brief Copy DataStore Group hierarchy (starting at root) and Buffer
   *        descriptions to given Conduit node.
   */
  void copyToConduitNode(Node& n) const;


  /*!
   * \brief Copy DataStore native layout (starting at root) to given Conduit node.
   *
   * The native layout is a Conduit Node hierarchy that maps the Conduit Node data
   * externally to the Sidre View data so that it can be filled in from the data
   * in the file (independent of file format) and can be accessed as a Conduit tree.
   */
  void createNativeLayout(Node& n) const;

  /*!
   * \brief Copy DataStore native layout (starting at root) to given Conduit node.
   *
   * The native layout is a Conduit Node hierarchy that maps the Conduit Node data
   * externally to the Sidre View data so that it can be filled in from the data
   * in the file (independent of file format) and can be accessed as a Conduit tree.
   */
  void createExternalLayout(Node& n) const;


  /*!
   * \brief Print JSON description of DataStore Group hierarchy (starting at
   *        root) and Buffer descriptions to std::cout.
   */
  void print() const;

  /*!
   * \brief Print JSON description of DataStore Group hierarchy (starting at
   *        root) and Buffer descriptions to given output stream.
   */
  void print(std::ostream& os) const;


  /// Developer notes:
  /// We should reduce these functions when SPIO is fully available ( in both serial and parallel ).
  /// We only need one or two simple save functions.  Try to keep this class simple and move the I/O
  /// interfaces to SPIO.

  /*!
   * \brief Save the DataStore to a new file.
   * \see DataGroup::save() for supported protocols
   */
  void save( const std::string& file_path,
             const std::string& protocol ) const;

  /*!
   * \brief Save the DataStore to an existing hdf5 handle.
   */
  void save( const hid_t& h5_id) const;

  /*!
   * \brief Load the DataStore from a file.
   */
  void load( const std::string& file_path,
             const std::string& protocol);

  /*!
   * \brief Load the DataStore from an hdf5 handle.
   */
  void load( const hid_t& h5_id);

  /*!
   * \brief Load the DataStore external data from a file.
   */
  void loadExternalData( const std::string& file_path,
                         const std::string& protocol);

  /*!
   * \brief Load the DataStore external data from an hdf5 handle.
   */
  void loadExternalData( const hid_t& h5_id);

private:
  /*!
   *  Unimplemented ctors and copy-assignment operators.
   */
#ifdef ATK_USE_CXX11
  DataStore( const DataStore& ) = delete;
  DataStore( DataStore&& ) = delete;

  DataStore& operator=( const DataStore& ) = delete;
  DataStore& operator=( DataStore&& ) = delete;
#else
  DataStore( const DataStore& );
  DataStore& operator=( const DataStore& );
#endif

  /// Root Group, created when DataStore object is created.
  DataGroup * m_RootGroup;

  /// Collection of Buffers in DataStore instance.
  std::vector<DataBuffer *> m_data_buffers;

  /// Collection of unused unique Buffer indices (they can be recycled).
  std::stack< IndexType > m_free_buffer_ids;

  /// Flag indicating whether SLIC logging environment was initialized in ctor.
  bool m_need_to_finalize_slic;
};


} /* end namespace sidre */
} /* end namespace asctoolkit */

#endif /* DATASTORE_HPP_ */

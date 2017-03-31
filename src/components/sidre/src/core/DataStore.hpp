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
 * \file DataStore.hpp
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

// Other axom headers
#include "common/AxomMacros.hpp"
#include "common/CommonTypes.hpp"
#include "slic/slic.hpp"

// Sidre project headers
#include "SidreTypes.hpp"

namespace axom
{
namespace sidre
{

class Buffer;
class Group;

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
  Group * getRoot()
  {
    return m_RootGroup;
  };

  /*!
   * \brief Return pointer to the root Group.
   */
  const Group * getRoot() const
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
             m_data_buffers[static_cast<unsigned>(idx)] != AXOM_NULLPTR );
  }

  /*!
   * \brief Return (non-const) pointer to Buffer object with given index,
   *        or AXOM_NULLPTR if none exists.
   */
  Buffer * getBuffer( IndexType idx ) const;

  /*!
   * \brief Create an undescribed Buffer object and return a pointer to it.
   *
   *        The Buffer must be described before it can be allocated.
   *
   *        The Buffer object is assigned a unique index when created and the
   *        Buffer object is owned by the DataStore object.
   */
  Buffer * createBuffer();

  /*!
   * \brief Create a Buffer object with specified type and number of
   *        elements and return a pointer to it.
   *
   *        See the Buffer::describe() method for valid data description.
   *
   *        The Buffer object is assigned a unique index when created and the
   *        Buffer object is owned by the DataStore object.
   */
  Buffer * createBuffer( TypeID type, SidreLength num_elems );

  /*!
   * \brief Remove Buffer from the DataStore and destroy it and
   *        its data.
   *
   *        Note that Buffer destruction detaches it from all Views to
   *        which it is attached.
   */
  void destroyBuffer( Buffer * buff );

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
   * \brief Print JSON description of DataStore Group hierarchy (starting at
   *        root) and Buffer descriptions to std::cout.
   */
  void print() const;

  /*!
   * \brief Print JSON description of DataStore Group hierarchy (starting at
   *        root) and Buffer descriptions to given output stream.
   */
  void print(std::ostream& os) const;

private:
  DISABLE_COPY_AND_ASSIGNMENT(DataStore);
  DISABLE_MOVE_AND_ASSIGNMENT(DataStore);

  /// Root Group, created when DataStore object is created.
  Group * m_RootGroup;

  /// Collection of Buffers in DataStore instance.
  std::vector<Buffer *> m_data_buffers;

  /// Collection of unused unique Buffer indices (they can be recycled).
  std::stack< IndexType > m_free_buffer_ids;

  /// Flag indicating whether SLIC logging environment was initialized in ctor.
  bool m_need_to_finalize_slic;
};


} /* end namespace sidre */
} /* end namespace axom */

#endif /* DATASTORE_HPP_ */

/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
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

// Other axom headers
#include "axom/Macros.hpp"
#include "axom/Types.hpp"
#include "slic/slic.hpp"

// Sidre project headers
#include "Attribute.hpp"
#include "SidreTypes.hpp"

namespace axom
{
namespace sidre
{

class Buffer;
class Group;
template <typename TYPE> class MapCollection;

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
  Group* getRoot()
  {
    return m_RootGroup;
  };

  /*!
   * \brief Return pointer to the root Group.
   */
  const Group* getRoot() const
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
  Buffer* getBuffer( IndexType idx ) const;

  /*!
   * \brief Create an undescribed Buffer object and return a pointer to it.
   *
   *        The Buffer must be described before it can be allocated.
   *
   *        The Buffer object is assigned a unique index when created and the
   *        Buffer object is owned by the DataStore object.
   */
  Buffer* createBuffer();

  /*!
   * \brief Create a Buffer object with specified type and number of
   *        elements and return a pointer to it.
   *
   *        See the Buffer::describe() method for valid data description.
   *
   *        The Buffer object is assigned a unique index when created and the
   *        Buffer object is owned by the DataStore object.
   */
  Buffer* createBuffer( TypeID type, SidreLength num_elems );

  /*!
   * \brief Remove Buffer from the DataStore and destroy it and
   *        its data.
   *
   *        Note that Buffer destruction detaches it from all Views to
   *        which it is attached.
   */
  void destroyBuffer( Buffer* buff );

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

//@{
//!  @name Methods to query, access, create, and destroy Attributes.

  /*!
   * \brief Return number of Attributes in the DataStore.
   */
  SidreLength getNumAttributes() const;

  /*!
   * \brief Create a Attribute object with a default scalar value.
   *
   *        The Attribute object is assigned a unique index when created and the
   *        Attribute object is owned by the DataStore object.
   */
  template<typename ScalarType>
  Attribute* createAttributeScalar( const std::string & name,
                                    ScalarType default_value)
  {
    Attribute* new_attribute = createAttributeEmpty(name);
    if ( new_attribute != AXOM_NULLPTR )
    {
      new_attribute->setDefaultScalar(default_value);
    }
    return new_attribute;
  }

  /*!
   * \brief Create a Attribute object with a default string value.
   *
   *        The Attribute object is assigned a unique index when created and the
   *        Attribute object is owned by the DataStore object.
   */
  Attribute* createAttributeString( const std::string & name,
                                    const std::string & default_value)
  {
    Attribute* new_attribute = createAttributeEmpty(name);
    if ( new_attribute != AXOM_NULLPTR )
    {
      new_attribute->setDefaultString(default_value);
    }
    return new_attribute;
  }

  /*!
   * \brief Return true if DataStore has created attribute name; else false.
   */
  bool hasAttribute( const std::string& name ) const;

  /*!
   * \brief Return true if DataStore has created attribute with index; else
   * false.
   */
  bool hasAttribute( IndexType idx ) const;

  /*!
   * \brief Remove Attribute from the DataStore and destroy it and
   *        its data.
   *
   * XXX    Note that Attribute destruction detaches it from all Views to
   *        which it is attached.
   */
  void destroyAttribute( const std::string & name );

  /*!
   * \brief Remove Attribute with given index from the DataStore and
   *        destroy it and its data.
   *
   *        Note that Attribute destruction detaches it from all Views to
   *        which it is attached.
   */
  void destroyAttribute( IndexType idx );

  /*!
   * \brief Remove Attribute from the DataStore and destroy it and
   *        its data.
   *
   * XXX    Note that Attribute destruction detaches it from all Views to
   *        which it is attached.
   */
  void destroyAttribute( Attribute* attr );

  /*!
   * \brief Remove all Attributes from the DataStore and destroy them
   *        and their data.
   *
   * XXX    Note that Attribute destruction detaches it from all Views to
   *        which it is attached.
   */
  void destroyAllAttributes();

//@}

//@{
//!  @name Attribute access and iteration methods.

  /*!
   * \brief Return pointer to non-const Attribute with given index.
   *
   * If no such Attribute exists, AXOM_NULLPTR is returned.
   */
  Attribute* getAttribute( IndexType idx );

  /*!
   * \brief Return pointer to const Attribute with given index.
   *
   * If no such Attribute exists, AXOM_NULLPTR is returned.
   */
  const Attribute* getAttribute( IndexType idx ) const;

  /*!
   * \brief Return pointer to non-const Attribute with given name.
   *
   * If no such Attribute exists, AXOM_NULLPTR is returned.
   */
  Attribute* getAttribute( const std::string& name );

  /*!
   * \brief Return pointer to const Attribute with given name.
   *
   * If no such Attribute exists, AXOM_NULLPTR is returned.
   */
  const Attribute* getAttribute( const std::string& name ) const;

  /*!
   * \brief Return first valid Attribute index in DataStore object
   *        (i.e., smallest index over all Attributes).
   *
   * sidre::InvalidIndex is returned if DataStore has no Attributes.
   */
  IndexType getFirstValidAttributeIndex() const;

  /*!
   * \brief Return next valid Attribute index in DataStore object after given
   * index
   *        (i.e., smallest index over all Attribute indices larger than given
   * one).
   *
   * sidre::InvalidIndex is returned if there is no valid index greater
   * than given one.
   */
  IndexType getNextValidAttributeIndex(IndexType idx) const;

  /*!
   * \brief Copy Attribute and default value to Conduit node.
   *        Return true if attributes were copied.
   */
  bool saveAttributeLayout(Node& node) const;

  /*!
   * \brief Create attributes from name/value pairs in node["attribute"].
   */
  void loadAttributeLayout(Node& node);

//@}
//----------------


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

//@{
//!  @name Private View declaration methods.
//!        (callable only by Group and View methods).

  /*!
   * \brief Create an Attribute and insert it into the DataStore.
   *        The attribute will be untyped.
   */
  Attribute* createAttributeEmpty(const std::string & name);

//@}

  /// Root Group, created when DataStore object is created.
  Group* m_RootGroup;

  /// Collection of Buffers in DataStore instance.
  std::vector<Buffer*> m_data_buffers;

  /// Collection of unused unique Buffer indices (they can be recycled).
  std::stack< IndexType > m_free_buffer_ids;

  ///////////////////////////////////////////////////////////////////
  //
  typedef MapCollection<Attribute> AttributeCollection;
  ///////////////////////////////////////////////////////////////////

  /// Collection of Attributes
  AttributeCollection* m_attribute_coll;

  /// Flag indicating whether SLIC logging environment was initialized in ctor.
  bool m_need_to_finalize_slic;
};


} /* end namespace sidre */
} /* end namespace axom */

#endif /* DATASTORE_HPP_ */

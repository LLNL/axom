// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file DataStore.hpp
 *
 * \brief   Header file containing definition of DataStore class.
 *
 ******************************************************************************
 */

#ifndef SIDRE_DATASTORE_HPP_
#define SIDRE_DATASTORE_HPP_

// Standard C++ headers
#include <vector>
#include <stack>

// Other axom headers
#include "axom/config.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/Types.hpp"
#include "axom/slic/interface/slic.hpp"

// Sidre project headers
#include "Attribute.hpp"
#include "SidreTypes.hpp"

namespace axom
{
namespace sidre
{
class Buffer;
class Group;
template <typename TYPE>
class MapCollection;

/*!
 * \class DataStore
 *
 * \brief DataStore is the main interface for creating and accessing
 *        Buffer objects.
 *
 * It maintains a collection of Buffer objects and owns the "root"
 * Group.  The initial name of the root Group is the empty string: a
 * code uses the getRoot() method to retrieve the root Group.  A Group
 * hierarchy (a tree) is created by creating child Groups of Groups.
 */
class DataStore
{
public:
  /*!
   * \brief Default ctor initializes DataStore object and creates a root Group.
   *
   * The ctor also initializes SLIC logging environment if it is not already
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
  Group* getRoot() { return m_RootGroup; };

  /*!
   * \brief Return pointer to the root Group.
   */
  const Group* getRoot() const { return m_RootGroup; };

  //@{
  //!  @name Methods to query, access, create, and destroy Buffers.

  /*!
   * \brief Return number of Buffers in the DataStore.
   */
  IndexType getNumBuffers() const
  {
    return static_cast<IndexType>(m_data_buffers.size() -
                                  m_free_buffer_ids.size());
  }

  /*!
   * \brief Return true if DataStore owns a Buffer with given index;
   *        else false.
   */
  bool hasBuffer(IndexType idx) const
  {
    return (0 <= idx && static_cast<unsigned>(idx) < m_data_buffers.size() &&
            m_data_buffers[static_cast<unsigned>(idx)] != nullptr);
  }

  /*!
   * \brief Return (non-const) pointer to Buffer object with the given
   *        index, or nullptr if none exists.
   */
  Buffer* getBuffer(IndexType idx) const;

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
  Buffer* createBuffer(TypeID type, IndexType num_elems);

  /*!
   * \brief Remove Buffer from the DataStore and destroy it and
   *        its data.
   *
   *        Note that Buffer destruction detaches it from all Views to
   *        which it is attached.
   */
  void destroyBuffer(Buffer* buff);

  /*!
   * \brief Remove Buffer with given index from the DataStore and
   *        destroy it and its data.
   *
   *        Note that Buffer destruction detaches it from all Views to
   *        which it is attached.
   */
  void destroyBuffer(IndexType idx);

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
  //!  @name Methods for iterating over Buffers in the DataStore.
  //!
  //! Using these methods, a code can get the first Buffer index and each
  //! succeeding index.  This allows Buffer iteration using the same
  //! constructs in C++, C, and Fortran.  Example:
  //!
  //!      for (sidre::IndexType idx = ds->getFirstValidBufferIndex();
  //!           sidre::indexIsValid(idx);
  //!           idx = ds->getNextValidBufferIndex(idx))
  //!      {
  //!          Buffer * buf = ds->getBuffer(idx);
  //!
  //!          /// code here using buf
  //!      }

  /*!
   * \brief Return first valid Buffer index.
   *
   * sidre::InvalidIndex is returned if Group has no Buffers.
   *
   * \sa axom::sidre::indexIsValid()
   */
  IndexType getFirstValidBufferIndex() const;

  /*!
   * \brief Return next valid Buffer index after given index.
   *
   * sidre::InvalidIndex is returned if there is no valid next index.
   *
   * \sa axom::sidre::indexIsValid()
   */
  IndexType getNextValidBufferIndex(IndexType idx) const;

  //@}

  //@{
  //!  @name Methods to query, access, create, and destroy Attributes.

  /*!
   * \brief Return number of Attributes in the DataStore.
   */
  IndexType getNumAttributes() const;

  /*!
   * \brief Create a Attribute object with a default scalar value.
   *
   *        The Attribute object is assigned a unique index when created and the
   *        Attribute object is owned by the DataStore object.
   */
  template <typename ScalarType>
  Attribute* createAttributeScalar(const std::string& name,
                                   ScalarType default_value)
  {
    Attribute* new_attribute = createAttributeEmpty(name);
    if(new_attribute != nullptr)
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
  Attribute* createAttributeString(const std::string& name,
                                   const std::string& default_value)
  {
    Attribute* new_attribute = createAttributeEmpty(name);
    if(new_attribute != nullptr)
    {
      new_attribute->setDefaultString(default_value);
    }
    return new_attribute;
  }

  /*!
   * \brief Return true if DataStore has created attribute name; else false.
   */
  bool hasAttribute(const std::string& name) const;

  /*!
   * \brief Return true if DataStore has created attribute with index; else
   * false.
   */
  bool hasAttribute(IndexType idx) const;

  /*!
   * \brief Remove Attribute from the DataStore and destroy it and
   *        its data.
   *
   * \note Destruction of an Attribute detaches it from all Views to
   *       which it is attached.
   */
  void destroyAttribute(const std::string& name);

  /*!
   * \brief Remove Attribute with given index from the DataStore and
   *        destroy it and its data.
   *
   * \note Destruction of an Attribute detaches it from all Views to
   *       which it is attached.
   */
  void destroyAttribute(IndexType idx);

  /*!
   * \brief Remove Attribute from the DataStore and destroy it and
   *        its data.
   *
   * \note Destruction of an Attribute detaches it from all Views to
   *       which it is attached.
   */
  void destroyAttribute(Attribute* attr);

  /*!
   * \brief Remove all Attributes from the DataStore and destroy them
   *        and their data.
   *
   * \note Destruction of an Attribute detaches it from all Views to
   *       which it is attached.
   */
  void destroyAllAttributes();

  //@}

  //@{
  //!  @name Attribute access methods.

  /*!
   * \brief Return pointer to non-const Attribute with given index.
   *
   * If no such Attribute exists, nullptr is returned.
   */
  Attribute* getAttribute(IndexType idx);

  /*!
   * \brief Return pointer to const Attribute with given index.
   *
   * If no such Attribute exists, nullptr is returned.
   */
  const Attribute* getAttribute(IndexType idx) const;

  /*!
   * \brief Return pointer to non-const Attribute with given name.
   *
   * If no such Attribute exists, nullptr is returned.
   */
  Attribute* getAttribute(const std::string& name);

  /*!
   * \brief Return pointer to const Attribute with given name.
   *
   * If no such Attribute exists, nullptr is returned.
   */
  const Attribute* getAttribute(const std::string& name) const;

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

  //@{
  //!  @name Methods for iterating over Attributes in the DataStore.
  //!
  //! Using these methods, a code can get the first Attribute index and each
  //! succeeding index.  This allows Attribute iteration using the same
  //! constructs in C++, C, and Fortran.  Example:
  //!
  //!      for (sidre::IndexType idx = ds->getFirstValidAttributeIndex();
  //!           sidre::indexIsValid(idx);
  //!           idx = ds->getNextValidAttributeIndex(idx))
  //!      {
  //!          Attribute * attr = ds->getAttribute(idx);
  //!
  //!          /// code here using attr
  //!      }

  /*!
   * \brief Return first valid Attribute index in DataStore object
   *        (i.e., smallest index over all Attributes).
   *
   * sidre::InvalidIndex is returned if DataStore has no Attributes.
   *
   * \sa axom::sidre::indexIsValid()
   */
  IndexType getFirstValidAttributeIndex() const;

  /*!
   * \brief Return next valid Attribute index in DataStore object after given
   * index (i.e., smallest index over all Attribute indices larger than given
   * one).
   *
   * sidre::InvalidIndex is returned if there is no valid index greater
   * than given one.
   *
   * \sa axom::sidre::indexIsValid()
   */
  IndexType getNextValidAttributeIndex(IndexType idx) const;

  //@}

  /*!
   * \brief Generate a Conduit Blueprint index based on a mesh in stored in
   *        this DataStore.
   *
   * If this DataStore contains data for a mesh that adheres to Conduit's
   * Blueprint format, this method generates a Blueprint index and stores
   * it at a specified location within this DataStore.  It uses as input a
   * path to a representative domain from that mesh.
   *
   * The domain must be held in a Group stored in this DataStore.  The location
   * of this Group is specified by the domain_path argument.  The generated
   * Blueprint index will be stored in a newly-created Group that will be
   * located at the path specified by the index_path argument.
   *
   * \param domain_path      path to a domain stored in the Blueprint format
   * \param mesh_name        name for the mesh to be described
   * \param index_path       path where the Blueprint index will be written
   *                         within this DataStore
   * \param num_domains      number of domains in the mesh
   *
   * \return true if the Blueprint index is successfully generated.
   */
  bool generateBlueprintIndex(const std::string& domain_path,
                              const std::string& mesh_name,
                              const std::string& index_path,
                              int num_domains);

#ifdef AXOM_USE_MPI
  /*!
   * \brief Generate a Conduit Blueprint index from a distributed mesh
   *        stored in this Datastore
   *
   * \param comm            communicator for the mesh distribution
   * \param domain_path     path where domains are located
   * \param mesh_name       name for the mesh to be described
   * \param index_path      path where index written
   */
  bool generateBlueprintIndex(MPI_Comm comm,
                              const std::string& domain_path,
                              const std::string& mesh_name,
                              const std::string& index_path);
#endif

  //----------------

  /*!
   * \brief Print JSON description of the DataStore Group hierarchy (starting
   *        at root) and Buffer descriptions to std::cout.
   */
  void print() const;

  /*!
   * \brief Print JSON description of the DataStore Group hierarchy (starting
   *        at root) and Buffer descriptions to given output stream.
   */
  void print(std::ostream& os) const;

  /*!
   * \brief Wires Conduit `info`, `warning`, and `error` messages to
   *        SLIC style handling.
   */
  static void setConduitSLICMessageHandlers();

  /*!
   * \brief Wires Conduit `info`, `warning`, and `error` messages to
   *        Conduit's default handling.
   */
  static void setConduitDefaultMessageHandlers();

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
  Attribute* createAttributeEmpty(const std::string& name);

  //@}

  /// Root Group, created when DataStore object is created.
  Group* m_RootGroup;

  /// Collection of Buffers in DataStore instance.
  std::vector<Buffer*> m_data_buffers;

  /// Collection of unused unique Buffer indices (they can be recycled).
  std::stack<IndexType> m_free_buffer_ids;

  ///////////////////////////////////////////////////////////////////
  //
  using AttributeCollection = MapCollection<Attribute>;
  ///////////////////////////////////////////////////////////////////

  /// Collection of Attributes
  AttributeCollection* m_attribute_coll;

  /// Flag indicating whether SLIC logging environment was initialized in ctor.
  bool m_need_to_finalize_slic;
};

} /* end namespace sidre */
} /* end namespace axom */

#endif /* SIDRE_DATASTORE_HPP_ */

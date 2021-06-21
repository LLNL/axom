// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Standard C++ headers
#include <fstream>

#include "conduit_blueprint.hpp"
#include "conduit_utils.hpp"  // for setting conduit's message logging handlers

// Associated header file
#include "DataStore.hpp"

// Other axom headers
#include "axom/slic/interface/slic.hpp"
#include "axom/slic/streams/GenericOutputStream.hpp"

// Sidre component headers
#include "MapCollection.hpp"
#include "Buffer.hpp"
#include "Group.hpp"
#include "Attribute.hpp"

#ifdef AXOM_USE_MPI
  #include "conduit_blueprint_mpi.hpp"
#endif

namespace axom
{
namespace sidre
{
/*
 *************************************************************************
 *
 * Callback function used to map Conduit errors to SLIC errors.
 *
 *************************************************************************
 */
void DataStoreConduitErrorHandler(const std::string& message,
                                  const std::string& fileName,
                                  int line)
{
  axom::slic::logErrorMessage(message, fileName, line);
}

/*
 *************************************************************************
 *
 * Callback function used to map Conduit warnings to SLIC warnings.
 *
 *************************************************************************
 */
void DataStoreConduitWarningHandler(const std::string& message,
                                    const std::string& fileName,
                                    int line)
{
  axom::slic::logWarningMessage(message, fileName, line);
}

/*
 *************************************************************************
 *
 * Callback function used to map Conduit info messages to SLIC info
 * messages.
 *
 *************************************************************************
 */
void DataStoreConduitInfoHandler(const std::string& message,
                                 const std::string& fileName,
                                 int line)
{
  axom::slic::logMessage(axom::slic::message::Info, message, fileName, line);
}

/*
 *************************************************************************
 *
 * DataStore ctor creates root Group.
 *
 *************************************************************************
 */
DataStore::DataStore()
  : m_RootGroup(nullptr)
  , m_attribute_coll(new AttributeCollection())
  , m_need_to_finalize_slic(false)
{
  if(!axom::slic::isInitialized())
  {
    axom::slic::initialize();

    std::string format =
      std::string("\n***********************************\n") +
      std::string("LEVEL=<LEVEL>\n") + std::string("MESSAGE=<MESSAGE>\n") +
      std::string("FILE=<FILE>\n") + std::string("LINE=<LINE>\n") +
      std::string("***********************************\n");

    axom::slic::setLoggingMsgLevel(axom::slic::message::Debug);
    axom::slic::addStreamToAllMsgLevels(
      new axom::slic::GenericOutputStream(&std::cout, format));

    m_need_to_finalize_slic = true;
  }

  // Provide SLIC error handler function to Conduit to log
  // internal Conduit errors.
  conduit::utils::set_error_handler(DataStoreConduitErrorHandler);
  conduit::utils::set_warning_handler(DataStoreConduitWarningHandler);
  conduit::utils::set_info_handler(DataStoreConduitInfoHandler);

  m_RootGroup = new Group("", this, false);
  m_RootGroup->m_parent = m_RootGroup;
};

/*
 *************************************************************************
 *
 * DataStore dtor destroys all contents.
 *
 *************************************************************************
 */
DataStore::~DataStore()
{
  // clean up Groups and Views before we destroy Buffers
  delete m_RootGroup;
  destroyAllBuffers();
  destroyAllAttributes();
  delete m_attribute_coll;

  if(m_need_to_finalize_slic)
  {
    axom::slic::finalize();
  }
}

/*
 *************************************************************************
 *
 * Re-wire Conduit Message Handlers to SLIC.
 *
 *************************************************************************
 */
void DataStore::setConduitSLICMessageHandlers()
{
  // Provide SLIC message handler functions to Conduit to log
  // internal Conduit info, warning, error messages.
  conduit::utils::set_error_handler(DataStoreConduitErrorHandler);
  conduit::utils::set_warning_handler(DataStoreConduitWarningHandler);
  conduit::utils::set_info_handler(DataStoreConduitInfoHandler);
}

/*
 *************************************************************************
 *
 * Restore Conduit Message Handlers to Conduit Defaults.
 *
 *************************************************************************
 */
void DataStore::setConduitDefaultMessageHandlers()
{
  // restore default handlers
  conduit::utils::set_info_handler(conduit::utils::default_info_handler);
  conduit::utils::set_warning_handler(conduit::utils::default_warning_handler);
  conduit::utils::set_error_handler(conduit::utils::default_error_handler);
}

/*
 *************************************************************************
 *
 * Return non-const pointer to Buffer with given index or null ptr.
 *
 *************************************************************************
 */
Buffer* DataStore::getBuffer(IndexType idx) const
{
  if(!hasBuffer(idx))
  {
    SLIC_CHECK_MSG(hasBuffer(idx),
                   "DataStore has no Buffer with index == " << idx);
    return nullptr;
  }

  return m_data_buffers[idx];
}

/*
 *************************************************************************
 *
 * Create new Buffer and assign unique id.
 *
 *************************************************************************
 */
Buffer* DataStore::createBuffer()
{
  // TODO: implement pool, look for free nodes.  Allocate in blocks.
  IndexType newIndex;
  if(m_free_buffer_ids.empty())
  {
    newIndex = m_data_buffers.size();
    m_data_buffers.push_back(nullptr);
  }
  else
  {
    newIndex = m_free_buffer_ids.top();
    m_free_buffer_ids.pop();
  }

  Buffer* const obj = new(std::nothrow) Buffer(newIndex);
  m_data_buffers[newIndex] = obj;

  return obj;
}

/*
 *************************************************************************
 *
 * Create new Buffer and assign unique id.
 *
 *************************************************************************
 */
Buffer* DataStore::createBuffer(TypeID type, IndexType num_elems)
{
  Buffer* buffer = createBuffer();

  if(buffer != nullptr)
  {
    buffer->describe(type, num_elems);
  }

  return buffer;
}

/*
 *************************************************************************
 *
 * Remove Buffer from the DataStore and destroy it, recover its
 * id for reuse.
 *
 *************************************************************************
 */
void DataStore::destroyBuffer(Buffer* buff)
{
  if(buff != nullptr)
  {
    buff->detachFromAllViews();
    IndexType idx = buff->getIndex();
    delete buff;
    SLIC_ASSERT(m_data_buffers[idx] != nullptr);
    m_data_buffers[idx] = nullptr;
    m_free_buffer_ids.push(idx);
  }
}

/*
 *************************************************************************
 *
 * Remove Buffer with given index from the DataStore and destroy it,
 * recover its id for reuse.
 *
 *************************************************************************
 */
void DataStore::destroyBuffer(IndexType idx) { destroyBuffer(getBuffer(idx)); }

/*
 *************************************************************************
 *
 * Destroy all Buffers in DataStore and reclaim indices.
 *
 *************************************************************************
 */
void DataStore::destroyAllBuffers()
{
  IndexType bidx = getFirstValidBufferIndex();
  while(indexIsValid(bidx))
  {
    destroyBuffer(bidx);
    bidx = getNextValidBufferIndex(bidx);
  }
}

/*
 *************************************************************************
 *
 * Return first valid Buffer index, or InvalidIndex if there is none.
 *
 *************************************************************************
 */
IndexType DataStore::getFirstValidBufferIndex() const
{
  return getNextValidBufferIndex(-1);
}

/*
 *************************************************************************
 *
 * Return first valid Buffer index, or InvalidIndex if there is none.
 *
 *************************************************************************
 */
IndexType DataStore::getNextValidBufferIndex(IndexType idx) const
{
  idx++;
  while(static_cast<unsigned>(idx) < m_data_buffers.size() &&
        m_data_buffers[idx] == nullptr)
  {
    idx++;
  }
  return ((static_cast<unsigned>(idx) < m_data_buffers.size()) ? idx
                                                               : InvalidIndex);
}

/*
 *************************************************************************
 *
 * Return number of Attributes in the DataStore.
 *
 *************************************************************************
 */
IndexType DataStore::getNumAttributes() const
{
  return m_attribute_coll->getNumItems();
}

/*
 *************************************************************************
 *
 * PRIVATE method to create an untyped Attribute object.
 *
 *************************************************************************
 */
Attribute* DataStore::createAttributeEmpty(const std::string& name)
{
  if(name.empty() || hasAttribute(name))
  {
    SLIC_CHECK(!name.empty());
    SLIC_CHECK_MSG(hasAttribute(name),
                   "Cannot create Attribute with name '"
                     << name
                     << " since it already has an Attribute with that name");
    return nullptr;
  }

  Attribute* new_attribute = new(std::nothrow) Attribute(name);
  if(new_attribute == nullptr)
  {
    return nullptr;
  }

  new_attribute->m_index = m_attribute_coll->insertItem(new_attribute, name);
  return new_attribute;
}

/*
 *************************************************************************
 *
 * Return true if this DataStore has an Attribute with given name; else false.
 *
 *************************************************************************
 */
bool DataStore::hasAttribute(const std::string& name) const
{
  return m_attribute_coll->hasItem(name);
}

/*
 *************************************************************************
 *
 * Return true if this DataStore has an Attribute with given index; else false.
 *
 *************************************************************************
 */
bool DataStore::hasAttribute(IndexType idx) const
{
  return m_attribute_coll->hasItem(idx);
}

/*
 *************************************************************************
 *
 * Destroy Attribute with given name.
 *
 *************************************************************************
 */
void DataStore::destroyAttribute(const std::string& name)
{
  Attribute* attr = getAttribute(name);
  destroyAttribute(attr);
}

/*
 *************************************************************************
 *
 * Destroy Attribute with given index.
 *
 *************************************************************************
 */
void DataStore::destroyAttribute(IndexType idx)
{
  Attribute* attr = m_attribute_coll->removeItem(idx);
  if(attr != nullptr)
  {
    delete attr;
  }
}

/*
 *************************************************************************
 *
 * Destroy Attribute.
 *
 *************************************************************************
 */
void DataStore::destroyAttribute(Attribute* attr)
{
  SLIC_ASSERT(attr != nullptr);

  destroyAttribute(attr->getIndex());
}

/*
 *************************************************************************
 *
 * Destroy all Attributes in DataStore and reclaim indices.
 *
 *************************************************************************
 */
void DataStore::destroyAllAttributes()
{
  IndexType bidx = getFirstValidAttributeIndex();
  while(indexIsValid(bidx))
  {
    destroyAttribute(bidx);
    bidx = getNextValidAttributeIndex(bidx);
  }
}

/*
 *************************************************************************
 *
 * Return pointer to non-const Attribute with given index.
 *
 * If no such Attribute exists, nullptr is returned.
 *
 *************************************************************************
 */
Attribute* DataStore::getAttribute(IndexType idx)
{
  SLIC_CHECK_MSG(hasAttribute(idx),
                 "DataStore has no Attribute with index " << idx);

  return m_attribute_coll->getItem(idx);
}

/*
 *************************************************************************
 *
 * Return pointer to const Attribute with given index.
 *
 * If no such Attribute exists, nullptr is returned.
 *
 *************************************************************************
 */
const Attribute* DataStore::getAttribute(IndexType idx) const
{
  SLIC_CHECK_MSG(hasAttribute(idx),
                 "DataStore has no Attribute with index " << idx);

  return m_attribute_coll->getItem(idx);
}

/*
 *************************************************************************
 *
 * Return pointer to non-const Attribute with given name.
 *
 * If no such Attribute exists, nullptr is returned.
 *
 *************************************************************************
 */
Attribute* DataStore::getAttribute(const std::string& name)
{
  SLIC_CHECK_MSG(hasAttribute(name),
                 "DataStore has no Attribute with name " << name);

  return m_attribute_coll->getItem(name);
}

/*
 *************************************************************************
 *
 * Return pointer to const Attribute with given name.
 *
 * If no such Attribute exists, nullptr is returned.
 *
 *************************************************************************
 */
const Attribute* DataStore::getAttribute(const std::string& name) const
{
  SLIC_CHECK_MSG(hasAttribute(name),
                 "DataStore has no Attribute with name " << name);

  return m_attribute_coll->getItem(name);
}

/*
 *************************************************************************
 *
 * \brief Return first valid Attribute index in DataStore object
 *        (i.e., smallest index over all Attributes).
 *
 * sidre::InvalidIndex is returned if DataStore has no Attributes.
 *************************************************************************
 */
IndexType DataStore::getFirstValidAttributeIndex() const
{
  return m_attribute_coll->getFirstValidIndex();
}

/*
 *************************************************************************
 *
 * \brief Return next valid Attribute index in DataStore object after given
 * index
 *        (i.e., smallest index over all Attribute indices larger than given
 * one).
 *
 * sidre::InvalidIndex is returned if there is no valid index greater
 * than given one.
 *************************************************************************
 */
IndexType DataStore::getNextValidAttributeIndex(IndexType idx) const
{
  return m_attribute_coll->getNextValidIndex(idx);
}

/*
 *************************************************************************
 *
 * Copy Attribute and default value to Conduit node.
 * Return true if attributes were copied.
 *
 *************************************************************************
 */
bool DataStore::saveAttributeLayout(Node& node) const
{
  // Adds a conduit node for this group if it has external views,
  // or if any of its children groups has an external view

  node.set(DataType::object());

  bool hasAttributes = false;

  IndexType aidx = getFirstValidAttributeIndex();
  while(indexIsValid(aidx))
  {
    const Attribute* attr = getAttribute(aidx);

    node[attr->getName()] = attr->getDefaultNodeRef();

    aidx = getNextValidAttributeIndex(aidx);
    hasAttributes = true;
  }

  return hasAttributes;
}

/*
 *************************************************************************
 *
 * Create attributes from name/value pairs in node["attribute"].  If no
 * attributes are available, do nothing. The values are used as the
 * default value of the name attribute.
 *
 *************************************************************************
 */
void DataStore::loadAttributeLayout(Node& node)
{
  if(node.has_path("attribute"))
  {
    conduit::NodeIterator attrs_itr = node["attribute"].children();
    while(attrs_itr.has_next())
    {
      Node& n_attr = attrs_itr.next();
      std::string attr_name = attrs_itr.name();

      Attribute* attr = createAttributeEmpty(attr_name);
      attr->setDefaultNodeRef(n_attr);
    }
  }
}

bool DataStore::generateBlueprintIndex(const std::string& domain_path,
                                       const std::string& mesh_name,
                                       const std::string& index_path,
                                       int num_domains)
{
  Group* domain =
    (domain_path == "/") ? getRoot() : getRoot()->getGroup(domain_path);

  conduit::Node mesh_node;
  domain->createNativeLayout(mesh_node);

  Group* bpindex = getRoot()->hasGroup(index_path)
    ? getRoot()->getGroup(index_path)
    : getRoot()->createGroup(index_path);

  bool success = false;
  conduit::Node info;
  if(conduit::blueprint::verify("mesh", mesh_node, info))
  {
    conduit::Node index;
    conduit::blueprint::mesh::generate_index(mesh_node,
                                             mesh_name,
                                             num_domains,
                                             index);

    bpindex->importConduitTree(index);

    success = true;
  }
  else
  {
    SLIC_DEBUG("Blueprint verify failed. Node info: " << info.to_string());
  }

  return success;
}

#ifdef AXOM_USE_MPI
bool DataStore::generateBlueprintIndex(MPI_Comm comm,
                                       const std::string& domain_path,
                                       const std::string& mesh_name,
                                       const std::string& index_path)
{
  Group* domain;
  if(domain_path == "/")
  {
    domain = getRoot();
  }
  else if(getRoot()->hasGroup(domain_path))
  {
    domain = getRoot()->getGroup(domain_path);
  }
  else
  {
    //There could be ranks with no domains--the MPI blueprint
    //functions can handle this.
    domain = nullptr;
  }

  conduit::Node mesh_node;
  if(domain)
  {
    domain->createNativeLayout(mesh_node);
  }

  Group* bpindex = getRoot()->hasGroup(index_path)
    ? getRoot()->getGroup(index_path)
    : getRoot()->createGroup(index_path);

  bool success = false;
  conduit::Node info;
  if(conduit::blueprint::mpi::verify("mesh", mesh_node, info, comm))
  {
    conduit::Node index;
    conduit::blueprint::mpi::mesh::generate_index(mesh_node, mesh_name, index, comm);

    Node& domain_rank_map = index["state/partition_map/datagroup"];
    conduit::blueprint::mpi::mesh::generate_domain_to_rank_map(mesh_node,
                                                               domain_rank_map,
                                                               comm);

    bpindex->importConduitTree(index);

    success = true;
  }
  else
  {
    SLIC_DEBUG("Blueprint verify failed. Node info: " << info.to_string());
  }

  return success;
}
#endif

/*
 *************************************************************************
 *
 * Print JSON description of Buffers and Group tree, starting at root,
 * to stdout.
 *
 *************************************************************************
 */
void DataStore::print() const { print(std::cout); }

/*
 *************************************************************************
 *
 * Print JSON description of Buffers and Group tree, starting at root,
 * to an ostream.
 *
 *************************************************************************
 */
void DataStore::print(std::ostream& os) const
{
  Node n;
  m_RootGroup->copyToConduitNode(n);
  n.to_json_stream(os);
}

} /* end namespace sidre */
} /* end namespace axom */

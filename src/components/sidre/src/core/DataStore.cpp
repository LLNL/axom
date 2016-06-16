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
 * \brief   Implementation file for DataStore class.
 *
 ******************************************************************************
 */

// Standard C++ headers
#include <fstream>

#include "conduit.hpp"
#include "relay.hpp"

// Associated header file
#include "DataStore.hpp"

// Sidre project headers
#include "DataBuffer.hpp"
#include "DataGroup.hpp"

// Other CS Toolkit headers
#include "slic/slic.hpp"
#include "slic/GenericOutputStream.hpp"

namespace asctoolkit
{
namespace sidre
{

/*
 *************************************************************************
 *
 * Function to map SLIC method for logging error message to Conduit
 * error handler.
 *
 *************************************************************************
 */
void DataStoreConduitErrorHandler( const std::string& message,
                                   const std::string& fileName,
                                   int line )
{
  slic::logErrorMessage( message, fileName, line );
}

/*
 *************************************************************************
 *
 * DataStore ctor creates root Group.
 *
 *************************************************************************
 */
DataStore::DataStore()
  : m_RootGroup(ATK_NULLPTR), m_need_to_finalize_slic(false)
{

  if ( !slic::isInitialized() )
  {
    slic::initialize();

    std::string format =
      std::string("\n***********************************\n")+
      std::string( "LEVEL=<LEVEL>\n" ) +
      std::string( "MESSAGE=<MESSAGE>\n" ) +
      std::string( "FILE=<FILE>\n" ) +
      std::string( "LINE=<LINE>\n" ) +
      std::string("***********************************\n");

    slic::setLoggingMsgLevel( slic::message::Debug );
    slic::addStreamToAllMsgLevels( new slic::GenericOutputStream(&std::cout,
                                                                 format) );

    m_need_to_finalize_slic = true;
  }

  // Provide SLIC error handler function to Conduit to log
  // internal Conduit errors.
  conduit::utils::set_error_handler( DataStoreConduitErrorHandler );

  m_RootGroup = new DataGroup("/", this);

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

  if ( m_need_to_finalize_slic )
  {
    slic::finalize();
  }
}


/*
 *************************************************************************
 *
 * Return non-const pointer to Buffer with given index or null ptr.
 *
 *************************************************************************
 */
DataBuffer * DataStore::getBuffer( IndexType idx ) const
{
  if ( !hasBuffer(idx) )
  {
    SLIC_CHECK_MSG(hasBuffer(idx),
                   "DataStore has no Buffer with index == " << idx);
    return ATK_NULLPTR;
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
DataBuffer * DataStore::createBuffer()
{
  // TODO: implement pool, look for free nodes.  Allocate in blocks.
  IndexType newIndex;
  if( m_free_buffer_ids.empty() )
  {
    newIndex = m_data_buffers.size();
    m_data_buffers.push_back( ATK_NULLPTR );
  }
  else
  {
    newIndex = m_free_buffer_ids.top();
    m_free_buffer_ids.pop();
  }

  DataBuffer * const obj = new(std::nothrow) DataBuffer( newIndex );
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
DataBuffer * DataStore::createBuffer( TypeID type, SidreLength num_elems )
{
  DataBuffer * buffer = createBuffer();

  if (buffer != ATK_NULLPTR)
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
void DataStore::destroyBuffer( DataBuffer * buff )
{
  if ( buff != ATK_NULLPTR )
  {
    buff->detachFromAllViews();
    IndexType idx = buff->getIndex();
    delete buff;
    SLIC_ASSERT( m_data_buffers[idx] != ATK_NULLPTR);
    m_data_buffers[idx] = ATK_NULLPTR;
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
void DataStore::destroyBuffer( IndexType idx )
{
  destroyBuffer(getBuffer(idx));
}

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
  while ( indexIsValid(bidx) )
  {
    destroyBuffer( bidx );
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
  while ( static_cast<unsigned>(idx) < m_data_buffers.size() &&
          m_data_buffers[idx] == ATK_NULLPTR )
  {
    idx++;
  }
  return ((static_cast<unsigned>(idx) < m_data_buffers.size()) ? idx
          : InvalidIndex);
}

/*
 *************************************************************************
 *
 * Copy Buffer descriptions and Group tree, starting at root, to given
 * Conduit node.
 *
 *************************************************************************
 */
void DataStore::copyToConduitNode(Node& n) const
{
  m_RootGroup->copyToConduitNode(n["DataStore/root"]);

  IndexType bidx = getFirstValidBufferIndex();
  while ( indexIsValid(bidx) )
  {
    Node& b = n["DataStore/buffers"].append();
    m_data_buffers[bidx]->copyToConduitNode(b);

    bidx = getNextValidBufferIndex(bidx);
  }
}


/*
 *************************************************************************
 *
 * Copy DataStore native layout, starting at root, to given Conduit node.
 *
 *************************************************************************
 */
void DataStore::createNativeLayout(Node& n) const
{
  m_RootGroup->createNativeLayout(n);
}


/*
 *************************************************************************
 *
 * Copy DataStore native external layout, starting at root, to given Conduit node.
 *
 *************************************************************************
 */
void DataStore::createExternalLayout(Node& n) const
{
  m_RootGroup->createExternalLayout(n);
}


/*
 *************************************************************************
 *
 * Print JSON description of Buffers and Group tree, starting at root,
 * to stdout.
 *
 *************************************************************************
 */
void DataStore::print() const
{
  print(std::cout);
}

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
  copyToConduitNode(n);
  n.to_json_stream(os);
}

/*************************************************************************/
// see ATK-735 - Add ability to control saving buffers/externals with save

void DataStore::save(const std::string& file_path,
                     const std::string& protocol,
                     const DataGroup * group) const
{

  SLIC_ERROR_IF(group != ATK_NULLPTR && group->getDataStore() != this,
                "Cannot call save method on Group not owned by this DataStore.");

  Node data_holder;
  data_holder.set_dtype(DataType::object());
  exportTo( group, data_holder);

  Node external_holder;
  external_holder.set_dtype(DataType::object());
  createExternalLayout(external_holder);

  Node two_part;
  two_part["sidre"] = data_holder;
  two_part["external"] = external_holder;

  if (protocol == "conduit")
  {
    conduit::relay::io::save(two_part, file_path);
  }
  else if (protocol == "conduit_hdf5")
  {
    conduit::relay::io::hdf5_write( two_part, file_path );
  }
  else if (protocol == "text")
  {
    std::ofstream output_file( file_path.c_str() );
    SLIC_ERROR_IF(!output_file.is_open(),
                  "Unable to create file " << file_path);
    if (output_file)
    {
      output_file << two_part.to_json();
    }
  }
  else
  {
    SLIC_ERROR("Invalid protocol " << protocol << " for file load.");
  }
}

/*************************************************************************/

void DataStore::save(const hid_t& h5_file_id,
                     const DataGroup * group) const
{
  SLIC_ERROR_IF(
    group != ATK_NULLPTR && group->getDataStore() != this,
    "Must call save function on Group that resides in this DataStore.");

  Node data_holder;
  data_holder.set_dtype(DataType::object());
  exportTo( group, data_holder);

  Node external_holder;
  external_holder.set_dtype(DataType::object());
  createExternalLayout(external_holder);

  Node two_part;
  two_part.set_dtype(DataType::object());
  two_part["sidre"] = data_holder;
  two_part["external"] = external_holder;

  conduit::relay::io::hdf5_write(two_part, h5_file_id);
}

/*************************************************************************/

/*
 *************************************************************************
 *
 * Load Group (including Views and child Groups) from a file
 *
 *************************************************************************
 */
void DataStore::load(const std::string& file_path,
                     const std::string& protocol,
                     DataGroup * group)
{
  SLIC_ERROR_IF(
    group != ATK_NULLPTR && group->getDataStore() != this,
    "Must call load function on Group that resides in this DataStore.");

  Node node;

  if (protocol == "conduit")
  {
    SLIC_ERROR("Invalid protocol " << protocol << " for file load.");
    conduit::relay::io::load(file_path, node);
  }
  else if (protocol == "conduit_hdf5")
  {
    conduit::relay::io::hdf5_read( file_path + ":sidre", node);
  }
  else
  {
    SLIC_ERROR("Invalid protocol " << protocol << " for file load.");
  }

  importFrom( group, node );

}

/*
 *************************************************************************
 *
 * Load Group (including Views and child Groups) from an hdf5 file
 *
 *************************************************************************
 */
void DataStore::load(const hid_t& h5_file_id,
                     DataGroup * group)
{
  SLIC_ERROR_IF(
    group != ATK_NULLPTR && group->getDataStore() != this,
    "Must call load function on Group that resides in this DataStore.");

  Node node;
  conduit::relay::io::hdf5_read(h5_file_id, "sidre", node);
  // for debugging call: n.print();
  importFrom( group, node );
}

/*
 *************************************************************************
 *
 * Load External Data from a file
 *
 *************************************************************************
 */
void DataStore::loadExternalData(const std::string& file_path,
                                 const std::string& protocol,
                                 DataGroup * group)
{
  SLIC_ERROR_IF(
    group != ATK_NULLPTR && group->getDataStore() != this,
    "Must call load function on Group that resides in this DataStore.");

  Node external_holder;
  external_holder.set_dtype(DataType::object());
  createExternalLayout(external_holder);

  if (protocol == "conduit")
  {
    SLIC_ERROR("Invalid protocol " << protocol << " for file load.");
    conduit::relay::io::load(file_path, external_holder);
  }
  else if (protocol == "conduit_hdf5")
  {
    conduit::relay::io::hdf5_read( file_path + ":external", external_holder);
  }
  else
  {
    SLIC_ERROR("Invalid protocol " << protocol << " for file load.");
  }
}

/*
 *************************************************************************
 *
 * Load External Data from an hdf5 file
 *
 *************************************************************************
 */
void DataStore::loadExternalData(const hid_t& h5_file_id,
                                 DataGroup * group)
{
  SLIC_ERROR_IF(
    group != ATK_NULLPTR && group->getDataStore() != this,
    "Must call load function on Group that resides in this DataStore.");

  Node external_holder;
  external_holder.set_dtype(DataType::object());
  createExternalLayout(external_holder);

  conduit::relay::io::hdf5_read(h5_file_id, "external", external_holder);
  // for debugging call: n.print();
}


/*
 *************************************************************************
 *
 * Serialize tree identified by a Group into a conduit node.  Include
 * any Buffers attached to Views in that tree.
 *
 * If Group is not specified, the DataStore root group will be used.
 *
 *************************************************************************
 */
void DataStore::exportTo(const DataGroup * group,
                         conduit::Node& data_holder) const
{
  if (group == ATK_NULLPTR)
  {
    group = getRoot();
  }

  // TODO - This implementation will change in the future.  We want to write
  // out some separate set of conduit nodes:
  // #1 A set of nodes representing the Group and Views (hierarchy), with
  // the data descriptions ( schemas ).
  // #2 A set of nodes for our data ( Buffers, external data, etc ).
  // On a load, we want to be able to create our DataStore tree first,
  // then call allocate ourself, then have conduit load the data directly
  // into our allocated memory areas.  Conduit can do this, as long as the
  // conduit node set is compatible with what's in the file.
  std::set<IndexType> buffer_indices;

  // Tell Group to add itself and all sub-Groups and Views to node.
  // Any Buffers referenced by those Views will be tracked in the
  // buffer_indices
  group->exportTo(data_holder, buffer_indices);

  if (!buffer_indices.empty())
  {
    // Now, add all the referenced buffers to the node.
    Node & bnode = data_holder["buffers"];
    for (std::set<IndexType>::iterator s_it = buffer_indices.begin() ;
         s_it != buffer_indices.end() ; ++s_it)
    {
      // Use a dictionary layout here instead of conduit list.
      // Conduit IO HDF5 doesn't support conduit list objects.
      std::ostringstream oss;
      oss << "buffer_id_" << *s_it;
      Node& buffer_holder = bnode.fetch( oss.str() );
      getBuffer( *s_it )->exportTo(buffer_holder);
    }
  }
}

/*
 *************************************************************************
 *
 * Imports tree from a conduit node into DataStore.  Includes
 * any Buffers attached to Views in that tree.
 *
 * If Group is not specified, the DataStore root Group will be used.
 *
 *************************************************************************
 */

void DataStore::importFrom(DataGroup * group,
                           conduit::Node& data_holder)
{
  // TODO - May want to put in a little meta-data into these files like a 'version'
  // or tag identifying the data.  We don't want someone giving us a file that
  // doesn't have our full multiView->buffer connectivity in there.

  if (group == ATK_NULLPTR)
  {
    group = getRoot();
  }

  group->destroyGroups();
  group->destroyViews();

  // First - Import Buffers into the DataStore.
  std::map<IndexType, IndexType> buffer_indices_map;

  // Added CON-132 ticket asking if has_path can just return false if node is empty or not an object type.
  if (data_holder.dtype().is_object() && data_holder.has_path("buffers"))
  {
    conduit::NodeIterator buffs_itr = data_holder["buffers"].children();
    while (buffs_itr.has_next())
    {
      Node& buffer_data_holder = buffs_itr.next();
      IndexType old_buffer_id = buffer_data_holder["id"].as_int32();

      DataBuffer * buffer = createBuffer();

      // track change of old Buffer id to new Buffer id
      buffer_indices_map[ old_buffer_id ] = buffer->getIndex();

      // populate the new Buffer's state
      buffer->importFrom(buffer_data_holder);
    }
  }

  // Next - import tree of Groups, sub-Groups, Views into the DataStore.
  // Use the mapping of old to new Buffer ids to connect the Views to the
  // right Buffers.
  group->importFrom(data_holder, buffer_indices_map);

}

} /* end namespace sidre */
} /* end namespace asctoolkit */

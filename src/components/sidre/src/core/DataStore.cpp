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

// Associated header file
#include "DataStore.hpp"

// Other toolkit component headers
#include "common/CommonTypes.hpp"

// SiDRe project headers
#include "DataBuffer.hpp"
#include "DataGroup.hpp"
#include "SidreTypes.hpp"

// Other CS Toolkit headers
#include "slic/slic.hpp"
#include "slic/GenericOutputStream.hpp"

#include "conduit.hpp"

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
 * Datastore ctor creates root group.
 *
 *************************************************************************
 */
DataStore::DataStore()
  : m_i_initialized_slic(false)
{

  // Initialize SLIC loggin environement, if not initialized already.
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

    m_i_initialized_slic = true;
  }

  // Provide SLIC error handler function to Conduit to log
  // internal Conduit errors.
  conduit::utils::set_error_handler( DataStoreConduitErrorHandler );

  m_RootGroup = new DataGroup("/", this);

};


/*
 *************************************************************************
 *
 * Datastore dtor destroys all contents.
 *
 *************************************************************************
 */
DataStore::~DataStore()
{
  // clean up groups and views before we destroy buffers
  delete m_RootGroup;
  destroyBuffers();

  if ( m_i_initialized_slic )
  {
    slic::finalize();
  }
}


/*
 *************************************************************************
 *
 * Return non-cost pointer to buffer with given index or null ptr.
 *
 *************************************************************************
 */
DataBuffer * DataStore::getBuffer( IndexType idx ) const
{

  if ( !hasBuffer(idx) )
  {
    SLIC_CHECK_MSG(hasBuffer(idx),
                   "Datastore has no buffer with index == " << idx);
    return ATK_NULLPTR;
  }

  return m_data_buffers[idx];
}

/*
 *************************************************************************
 *
 * Create new data buffer and assign unique id.
 *
 *************************************************************************
 */
DataBuffer * DataStore::createBuffer()
{
  // TODO: implement pool, look for free nodes.  Allocate in blocks.
  IndexType newIndex = m_data_buffers.size();
  m_data_buffers.push_back( ATK_NULLPTR );
  if( !m_free_buffer_ids.empty() )
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
 * Create new data buffer and assign unique id.
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
 * Remove data buffer with given index from the datastore and destroy it,
 * recover its id for reuse.
 *
 *************************************************************************
 */
void DataStore::destroyBuffer( IndexType idx )
{
  if ( !hasBuffer(idx) || m_data_buffers[idx] == ATK_NULLPTR ||
       m_data_buffers[idx]->getNumViews() != 0 )
  {
    SLIC_CHECK_MSG( hasBuffer(idx), "No buffer found with index " << idx);
    SLIC_CHECK_MSG( m_data_buffers[idx] != ATK_NULLPTR,
                    "Datastore has NULL pointer for buffer index " << idx);
    SLIC_CHECK_MSG( m_data_buffers[idx] != ATK_NULLPTR &&
                    m_data_buffers[idx]->getNumViews() != 0,
                    "Unable to delete buffer, it has " <<
                    m_data_buffers[idx]->getNumViews() << " still attached.");
    return;
  }

  delete m_data_buffers[idx];
  m_data_buffers[idx] = ATK_NULLPTR;
  m_free_buffer_ids.push(idx);
}


/*
 *************************************************************************
 *
 * Destroy all buffers in datastore and reclaim indices.
 *
 *************************************************************************
 */
void DataStore::destroyBuffers()
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
 * Remove data buffer with given index from the datastore leaving it intact,
 * and return a pointer to it. Its index is recovered for reuse.
 *
 *************************************************************************
 */
DataBuffer * DataStore::detachBuffer( IndexType idx )
{
  SLIC_CHECK_MSG(hasBuffer(idx), "no buffer exists with index == " << idx);

  DataBuffer * rval = ATK_NULLPTR;

  if ( hasBuffer(idx) )
  {
    rval = m_data_buffers[idx];
    m_data_buffers[idx] = ATK_NULLPTR;
    m_free_buffer_ids.push(idx);
  }

  return rval;
}

/*
 *************************************************************************
 *
 * Return first valid buffer index, or InvalidIndex if there is none.
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
 * Return first valid buffer index, or InvalidIndex if there is none.
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
 * Copy buffer descriptions and group tree, starting at root, to given
 * Conduit node.
 *
 *************************************************************************
 */
void DataStore::info(Node& n) const
{
  m_RootGroup->info(n["DataStore/root"]);

  IndexType bidx = getFirstValidBufferIndex();
  while ( indexIsValid(bidx) )
  {
    Node& b = n["DataStore/buffers"].append();
    m_data_buffers[bidx]->info(b);

    bidx = getNextValidBufferIndex(bidx);
  }
}

/*
 *************************************************************************
 *
 * Print JSON description of data buffers and group tree, starting at root,
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
 * Print JSON description of data buffers and group tree, starting at root,
 * to an ostream.
 *
 *************************************************************************
 */
void DataStore::print(std::ostream& os) const
{
  Node n;
  info(n);
  n.to_json_stream(os);
}

/*************************************************************************/

void DataStore::save(const std::string& obase,
                     const std::string& protocol,
                     const DataGroup* group) const
{
  if (protocol == "conduit")
  {
    Node data_holder;
    exportTo( group, data_holder, protocol );

    // for debugging call: n.print();
    conduit::io::save(data_holder, obase);
  }
}

/*************************************************************************/

void DataStore::save(const std::string& obase,
                     const std::string& protocol,
                     const hid_t& h5_file_id,
                     const DataGroup * group) const
{
  if (protocol == "conduit_hdf5")
  {
    Node data_holder;
    exportTo(group, data_holder, protocol);

    std::string file_path;
    std::string hdf5_path;
    conduit::utils::split_string(obase,
                                 std::string(":"),
                                 file_path,
                                 hdf5_path);

    conduit::io::hdf5_write(data_holder, h5_file_id, hdf5_path);
  }
}

/*************************************************************************/

/*
 *************************************************************************
 *
 * Load data group (including data views and child groups) from a file
 * set named "obase" into this group object.
 *
 * Note: Only valid protocol is "conduit".
 *
 *************************************************************************
 */
void DataStore::load(const std::string& obase,
                     const std::string& protocol,
                     DataGroup * group)
{
  if (protocol == "conduit")
  {
    Node node;
    conduit::io::load(obase, node);
    // for debugging call: n.print();
    importFrom( group, node );
  }
}

/*
 *************************************************************************
 *
 * Load data group (including data views and child groups) from an hdf5 file
 * named "obase" into this group object.
 *
 * Note: Only valid protocol is "conduit_hdf5".
 *
 *************************************************************************
 */
void DataStore::load(const std::string& obase,
                     const std::string& protocol,
                     const hid_t& h5_file_id,
                     DataGroup * group)
{
  if (protocol == "conduit_hdf5")
  {
    std::string file_path;
    std::string hdf5_path;
    conduit::utils::split_string(obase,
                                 std::string(":"),
                                 file_path,
                                 hdf5_path);
    Node node;
    conduit::io::hdf5_read(h5_file_id, hdf5_path, node);
    // for debugging call: n.print();
    importFrom( group, node );
  }
}


/*
 *************************************************************************
 *
 * Serialize tree identified by a group into a conduit node.  Include
 * any buffers attached to views in that tree.
 *
 * If group is not specified, the datastore root group will be used.
 *
 *************************************************************************
 */
void DataStore::exportTo(const DataGroup * group,
               conduit::Node& data_holder,
               const std::string& protocol) const
{
  if (group == ATK_NULLPTR)
  {
    group = getRoot();
  }

  //TODO - Use the protocol to choose whether or not we want to save out our
  // full multiview->buffer connectivity (for restarts), or whether we want
  // to collapse those and have each view have a copy of the buffer data
  // ( for exporting to other consumers ).

  std::set<IndexType> buffer_indices;

  // Tell group to add itself and all sub-groups and views to node.
  // Any buffers referenced by those views will be tracked in the
  // buffer_indices
  group->exportTo(data_holder, buffer_indices);

  // Now, add all those referenced buffers to the node.

  // TODO - Conduit lists are supported by the conduit protocol, but not
  //        the conduit_hdf5 protocol.  May need to re-implement this
  //        to use a dictionary/map instead of list layout if conduit_hdf5
  //        doesn't add support for conduit lists.
  for (std::set<IndexType>::iterator s_it = buffer_indices.begin();
       s_it != buffer_indices.end(); ++s_it)
  {
    conduit::Node& buffer_holder = data_holder["buffers"].append();
    getBuffer( *s_it )->exportTo(buffer_holder);
  }
}

/*
 *************************************************************************
 *
 * Imports tree from a conduit node into datastore.  Includes
 * any buffers attached to views in that tree.
 *
 * If group is not specified, the datastore root group will be used.
 *
 *************************************************************************
 */

void DataStore::importFrom(DataGroup * group,
                           conduit::Node& data_holder)
{
  // TODO - May want to put in a little meta-data into these files like a 'version'
  // or tag identifying the data.  We don't want someone giving us a file that
  // doesn't have our full multiview->buffer connectivity in there.

  if (group == ATK_NULLPTR)
  {
    group = getRoot();
  }

  group->destroyGroups();
  group->destroyViews();

  std::cerr << " importing data into group " << group->getName() << std::endl;

  // First - Import buffers into the datastore.
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

      // track change of old buffer id to new buffer id
      buffer_indices_map[ old_buffer_id ] = buffer->getIndex();

      // populate the new buffer's state
      buffer->importFrom(buffer_data_holder);
    }
  }

  // Next - import tree of groups, sub-groups, views into the datastore.
  // Use the mapping of old to new buffer ids to connect the views to the
  // right buffers.
  group->importFrom(data_holder, buffer_indices_map);

}

} /* end namespace sidre */
} /* end namespace asctoolkit */

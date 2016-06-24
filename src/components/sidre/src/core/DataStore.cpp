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

/*
 *************************************************************************
 *
 * Save Group (including Views and child Groups) to a file
 *
 *************************************************************************
 */
void DataStore::save(const std::string& file_path,
                     const std::string& protocol) const
{
  getRoot()->save(file_path,protocol);
}

/*
 *************************************************************************
 *
 * Save Group (including Views and child Groups) to a hdf5 handle
 *
 *************************************************************************
 */
void DataStore::save(const hid_t& h5_id) const
{
  getRoot()->save(h5_id);
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
                     const std::string& protocol)
{
  getRoot()->load(file_path,protocol);
}

/*
 *************************************************************************
 *
 * Load Group (including Views and child Groups) from an hdf5 id
 *
 *************************************************************************
 */
void DataStore::load(const hid_t& h5_id)
{
  getRoot()->load(h5_id);
}

/*
 *************************************************************************
 *
 * Load External Data from a file
 *
 *************************************************************************
 */
void DataStore::loadExternalData(const std::string& file_path,
                                 const std::string& protocol)
{
  getRoot()->loadExternalData(file_path,protocol);
}

/*
 *************************************************************************
 *
 * Load External Data from an hdf5 handle
 *
 *************************************************************************
 */
void DataStore::loadExternalData(const hid_t& h5_id)
{
  getRoot()->loadExternalData(h5_id);
}

} /* end namespace sidre */
} /* end namespace asctoolkit */

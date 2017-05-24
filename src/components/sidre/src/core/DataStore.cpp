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
 * \file DataStore.cpp
 *
 * \brief   Implementation file for DataStore class.
 *
 ******************************************************************************
 */

// Standard C++ headers
#include <fstream>

#include "conduit_utils.hpp" // for setting conduit's message logging handlers

// Associated header file
#include "DataStore.hpp"

// Other axom headers
#include "slic/slic.hpp"
#include "slic/GenericOutputStream.hpp"

// Sidre component headers
#include "MapCollection.hpp"
#include "Buffer.hpp"
#include "Group.hpp"
#include "Attribute.hpp"


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
void DataStoreConduitErrorHandler( const std::string& message,
                                   const std::string& fileName,
                                   int line )
{
  axom::slic::logErrorMessage( message, fileName, line );
}


/*
 *************************************************************************
 *
 * Callback function used to map Conduit warnings to SLIC warnings.
 *
 *************************************************************************
 */
void DataStoreConduitWarningHandler( const std::string& message,
                                     const std::string& fileName,
                                     int line )
{
  axom::slic::logWarningMessage( message, fileName, line );
}

/*
 *************************************************************************
 *
 * Callback function used to map Conduit info messages to SLIC info
 * messages.
 *
 *************************************************************************
 */
void DataStoreConduitInfoHandler( const std::string& message,
                                  const std::string& fileName,
                                  int line )
{
  axom::slic::logMessage( axom::slic::message::Info, message, fileName, line );
}

/*
 *************************************************************************
 *
 * DataStore ctor creates root Group.
 *
 *************************************************************************
 */
DataStore::DataStore()
  : m_RootGroup(AXOM_NULLPTR),
    m_attribute_coll(new AttributeCollection()),
    m_need_to_finalize_slic(false)
{

  if ( !axom::slic::isInitialized() )
  {
    axom::slic::initialize();

    std::string format =
      std::string("\n***********************************\n")+
      std::string( "LEVEL=<LEVEL>\n" ) +
      std::string( "MESSAGE=<MESSAGE>\n" ) +
      std::string( "FILE=<FILE>\n" ) +
      std::string( "LINE=<LINE>\n" ) +
      std::string("***********************************\n");

    axom::slic::setLoggingMsgLevel( axom::slic::message::Debug );
    axom::slic::addStreamToAllMsgLevels(
      new axom::slic::GenericOutputStream(&std::cout,format) );

    m_need_to_finalize_slic = true;
  }

  // Provide SLIC error handler function to Conduit to log
  // internal Conduit errors.
  conduit::utils::set_error_handler( DataStoreConduitErrorHandler );
  conduit::utils::set_warning_handler( DataStoreConduitWarningHandler );
  conduit::utils::set_info_handler( DataStoreConduitInfoHandler );

  m_RootGroup = new Group("", this);
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
  delete m_attribute_coll;

  if ( m_need_to_finalize_slic )
  {
    axom::slic::finalize();
  }
}


/*
 *************************************************************************
 *
 * Return non-const pointer to Buffer with given index or null ptr.
 *
 *************************************************************************
 */
Buffer * DataStore::getBuffer( IndexType idx ) const
{
  if ( !hasBuffer(idx) )
  {
    SLIC_CHECK_MSG(hasBuffer(idx),
                   "DataStore has no Buffer with index == " << idx);
    return AXOM_NULLPTR;
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
Buffer * DataStore::createBuffer()
{
  // TODO: implement pool, look for free nodes.  Allocate in blocks.
  IndexType newIndex;
  if( m_free_buffer_ids.empty() )
  {
    newIndex = m_data_buffers.size();
    m_data_buffers.push_back( AXOM_NULLPTR );
  }
  else
  {
    newIndex = m_free_buffer_ids.top();
    m_free_buffer_ids.pop();
  }

  Buffer * const obj = new(std::nothrow) Buffer( newIndex );
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
Buffer * DataStore::createBuffer( TypeID type, SidreLength num_elems )
{
  Buffer * buffer = createBuffer();

  if (buffer != AXOM_NULLPTR)
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
void DataStore::destroyBuffer( Buffer * buff )
{
  if ( buff != AXOM_NULLPTR )
  {
    buff->detachFromAllViews();
    IndexType idx = buff->getIndex();
    delete buff;
    SLIC_ASSERT( m_data_buffers[idx] != AXOM_NULLPTR);
    m_data_buffers[idx] = AXOM_NULLPTR;
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
          m_data_buffers[idx] == AXOM_NULLPTR )
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
size_t DataStore::getNumAttributes() const
{
  return m_attribute_coll->getNumItems();
}

Attribute * DataStore::createAttribute( const std::string & name )
{
  // XXX look for existing attribute

  Attribute * new_attribute = new(std::nothrow) Attribute(name);
  if ( new_attribute == AXOM_NULLPTR )
  {
    return AXOM_NULLPTR;
  }

  m_attribute_coll->insertItem(new_attribute, name);
  IndexType idx = m_attribute_coll->getItemIndex(name);
  new_attribute->setIndex(idx);
  return new_attribute;
}

/*
 *************************************************************************
 *
 * Return true if this DataStore has an Attribute with given name; else false.
 *
 *************************************************************************
 */
bool DataStore::hasAttribute( const std::string & name ) const
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
bool DataStore::hasAttribute( IndexType idx ) const
{
  return m_attribute_coll->hasItem(idx);
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
  m_RootGroup->copyToConduitNode(n);
  n.to_json_stream(os);
}

} /* end namespace sidre */
} /* end namespace axom */

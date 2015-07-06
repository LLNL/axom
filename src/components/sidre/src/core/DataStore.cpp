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


namespace asctoolkit
{
namespace sidre
{


/*
 *************************************************************************
 *
 * Datastore ctor creates root group.
 *
 *************************************************************************
 */
DataStore::DataStore()
{
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
}


/*
 *************************************************************************
 *
 * Return non-cost pointer to buffer with given index or null ptr.
 *
 *************************************************************************
 */
DataBuffer * DataStore::getBuffer( IndexType idx )
{
  SLIC_CHECK_MSG(hasBuffer(idx), "no buffer exists with index == " << idx);

  if ( hasBuffer(idx) ) 
  { 
    return m_data_buffers[idx];
  }
  else 
  { 
    return ATK_NULLPTR;
  }
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
  DataBuffer * const obj = new DataBuffer( newIndex );

  m_data_buffers[newIndex] = obj;

  return obj;
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
  SLIC_CHECK_MSG(hasBuffer(idx), "no buffer exists with index == " << idx);

  if ( hasBuffer(idx) ) 
  {
    delete m_data_buffers[idx];
    m_data_buffers[idx] = ATK_NULLPTR;
    m_free_buffer_ids.push(idx);
  }
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
  /// TODO: after conduit update, use new ostream variant of to_json.
  std::ostringstream oss;
  n.json_to_stream(oss);
  os << oss.str();
}


} /* end namespace sidre */
} /* end namespace asctoolkit */

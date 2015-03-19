/*
 * DataGroup.cpp
 *
 *  Created on: Dec 2, 2014
 *      Author: settgast
 */

#include "DataStore.hpp"
#include "DataBuffer.hpp"

namespace DataStoreNS
{

  DataStore::DataStore() :
    m_RootGroup(nullptr, this),
    m_DataBuffers(),
    m_AvailableDataBuffers()
    {};

  DataStore::~DataStore()
  {
    for( dataBufferContainerType::iterator iter=m_DataBuffers.begin() ; iter!=m_DataBuffers.end() ; ++iter )
    {
      delete *iter;
    }
  }

  DataBuffer* DataStore::CreateDataBuffer()
  {
    // TODO: implement pool, look for free nodes.  Allocate in blocks.
    IDType newIndex = m_DataBuffers.size();
    m_DataBuffers.push_back( nullptr );
    if( m_AvailableDataBuffers.empty() )
    {
      newIndex = m_AvailableDataBuffers.top();
      m_AvailableDataBuffers.pop();
    }
    DataBuffer* const obj = new DataBuffer( newIndex );

    m_DataBuffers[newIndex] = obj ;

    return obj;
  }

  void DataStore::DeleteDataBuffer( const IDType id )
  {
    delete m_DataBuffers[id];
    m_AvailableDataBuffers.push(id);
  }

  DataBuffer* DataStore::DetatchDataBuffer( const IDType id )
  {
    DataBuffer* const rval = m_DataBuffers[id];
    m_DataBuffers[id] = nullptr;
    m_AvailableDataBuffers.push(id);

    return rval;
  }

} /* namespace DataStoreNS */

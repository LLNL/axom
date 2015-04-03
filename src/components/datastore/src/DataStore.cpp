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
    m_RootGroup("/", this),
    m_DataBuffers(),
    m_AvailableDataBuffers()
    {};

  DataStore::~DataStore()
  {
      DestroyBuffers();
  }

  DataBuffer* DataStore::CreateBuffer()
  {
    // TODO: implement pool, look for free nodes.  Allocate in blocks.
    IDType newIndex = m_DataBuffers.size();
    m_DataBuffers.push_back( nullptr );
    if( !m_AvailableDataBuffers.empty() )
    {
      newIndex = m_AvailableDataBuffers.top();
      m_AvailableDataBuffers.pop();
    }
    DataBuffer* const obj = new DataBuffer( newIndex );

    m_DataBuffers[newIndex] = obj ;

    return obj;
  }

  void DataStore::DestroyBuffer( const IDType id )
  {
    delete m_DataBuffers[id];
    m_DataBuffers[id] = nullptr;
    std::cout<<m_DataBuffers[id]<<std::endl;
    m_AvailableDataBuffers.push(id);
  }

  DataBuffer* DataStore::DetatchBuffer( const IDType id )
  {
    DataBuffer* const rval = m_DataBuffers[id];
    m_DataBuffers[id] = nullptr;
    m_AvailableDataBuffers.push(id);

    return rval;
  }
  
  void DataStore::DestroyBuffers()
  {
      for( dataBufferContainerType::iterator iter=m_DataBuffers.begin() ;                  iter!=m_DataBuffers.end() ; ++iter )
      {
        delete *iter;
      }
  }

  void DataStore::Print() const
  {
      Node n;
      Info(n);
      n.print();
  }

  void DataStore::Info(Node &n) const
  {
      m_RootGroup.Info(n["DataStore/root"]);
      for( dataBufferContainerType::const_iterator iter=m_DataBuffers.begin() ;
           iter!=m_DataBuffers.end() ;
           ++iter )
      {
          Node &b = n["DataStore/buffers"].append();
          if(*iter != nullptr)
          {
              (*iter)->Info(b);
          }
      }
  }


} /* namespace DataStoreNS */

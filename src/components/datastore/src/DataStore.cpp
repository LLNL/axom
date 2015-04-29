/*
 * DataGroup.cpp
 *
 *  Created on: Dec 2, 2014
 *      Author: settgast
 */

#include "DataStore.hpp"
#include "DataBuffer.hpp"

namespace sidre
{

  DataStore::DataStore() :
    m_DataBuffers(),
    m_AvailableDataBuffers()
    {
        m_RootGroup = new DataGroup("/",this);
    };

  DataStore::~DataStore()
  {
      // clean up views before we destroy buffers
      delete m_RootGroup;
      destroyBuffers();
  }

  DataBuffer* DataStore::createBuffer()
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

  void DataStore::destroyBuffer( const IDType id )
  {
    delete m_DataBuffers[id];
    m_DataBuffers[id] = nullptr;
    m_AvailableDataBuffers.push(id);
  }

  DataBuffer* DataStore::detachBuffer( const IDType id )
  {
    DataBuffer* const rval = m_DataBuffers[id];
    m_DataBuffers[id] = nullptr;
    m_AvailableDataBuffers.push(id);

    return rval;
  }
  
  void DataStore::destroyBuffers()
  {
      for( std::vector<DataBuffer*>::iterator iter=m_DataBuffers.begin() ;
           iter!=m_DataBuffers.end() ; ++iter )
      {
        delete *iter;
      }
  }

  void DataStore::print() const
  {
      Node n;
      info(n);
      n.print();
  }

  void DataStore::info(Node &n) const
  {
      m_RootGroup->info(n["DataStore/root"]);
      for( std::vector<DataBuffer*>::const_iterator iter=m_DataBuffers.begin() ;
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


} /* namespace sidre */

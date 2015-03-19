/*
 * DataGroup.cpp
 *
 *  Created on: Dec 2, 2014
 *      Author: settgast
 */

#include "DataStore.hpp"
#include "DataObject.hpp"

namespace DataStoreNS
{

  DataStore::~DataStore()
  {
    for( dataObjectContainerType::iterator iter=m_DataBuffers.begin() ; iter!=m_DataBuffers.end() ; ++iter )
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

  void DataStore::DeleteDataBuffer( DataBuffer*& obj )
  {
/*
    DataObject::GroupContainerType *gset = obj->GetDataGroups();
    for( DataObject::GroupContainerType::iterator it = gset->begin() ;
         it != gset->end() ; ++it )
    {
      DataGroup *grp = *it;
      grp->RemoveDataView( obj );
    }

    m_AvailableDataObjects.push( obj );
    return;
    */
  }

} /* namespace DataStoreNS */

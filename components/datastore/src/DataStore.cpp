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
  for (size_t i=0; i < m_DataObjects.size(); i++) {
    delete m_DataObjects[i];
  }
}

DataObject *DataStore::CreateDataObject()
{
  // TODO: implement pool, look for free nodes.  Allocate in blocks.
  DataObject *obj;
  if (m_AvailableDataObjects.empty()) {
    obj = new DataObject(m_IDCounter);
    ++m_IDCounter;
    m_DataObjects.push_back(obj);
  } else {
    obj = m_AvailableDataObjects.top();
    m_AvailableDataObjects.pop();
  }

  return obj;
}

void DataStore::DeleteDataObject( DataObject *obj )
{
  DataObject::GroupContainerType *gset = obj->GetDataGroups();
  for (DataObject::GroupContainerType::iterator it = gset->begin(); it != gset->end(); ++it) {
    DataGroup *grp = *it;
    grp->RemoveDataObject(obj);
  }

  m_AvailableDataObjects.push(obj);
  return;
}


} /* namespace DataStoreNS */

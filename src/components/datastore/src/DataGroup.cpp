/*
 * DataGroup.cpp
 *
 *  Created on: Dec 2, 2014
 *      Author: settgast
 */

#include "DataGroup.hpp"
#include "DataStore.hpp"


namespace DataStoreNS
{


bool DataGroup::HasName(const std::string& name) {
  // XXX must perform two lookups: objects and groups
  DataGroup::lookupType::iterator it = m_DataObjectLookup.find(name);
  if (it != m_DataObjectLookup.end())
    return true;
  DataGroup::lookupGroup::iterator itg = m_childGroups.find(name);
  if (itg != m_childGroups.end())
    return true;
  return false;
}

DataObject *DataGroup::AddDataObject(const std::string& name, DataObject *obj)
{
  if (HasName(name)) {
      throw std::exception();
  }
  if (obj == nullptr ) {
    obj = m_datastore->CreateDataObject();
  }
  m_DataObjectLookup[name] = m_DataObjects.size();  // map name to index
  m_DataObjects.push_back(obj);
  obj->AttachGroup(this);
  // XXX how does user get index? Search IndexDataObject(name) or return from here?
  return obj;
}

DataObject *DataGroup::RemoveDataObject(const std::string& name)
{
  DataGroup::lookupType::iterator it = m_DataObjectLookup.find(name);
  if (it != m_DataObjectLookup.end()) {
    IDType indx = it->second;
    DataObject *obj = m_DataObjects[indx];
    obj->DetachGroup(this);
    m_DataObjectLookup.erase(it);
    m_DataObjects[indx] = nullptr;    // XXX remove from m_DataObjects
    // XXX this makes a hole in the table, but perserved existing indexes.
  } else {
    throw std::exception();
  }
}

DataObject *DataGroup::RemoveDataObject(DataObject *obj)
{
  std::string const & name = NameDataObject(obj);
  RemoveDataObject(name);
}

DataGroup *DataGroup::CreateDataGroup(const std::string& name)
{
  if (HasName(name)) {
      throw std::exception();
  }
  DataGroup *grp = new DataGroup(this, this->m_datastore);
  m_childGroups[name] = grp;
  return grp;
}

void DataGroup::ClearDataObjects()
{
  m_DataObjects.clear();
  m_DataObjectLookup.clear();
}

std::string const & DataGroup::NameDataObject(DataObject *obj)
{
  DataGroup::lookupType::iterator it;

  // XXX brute force search for obj
  for (it=m_DataObjectLookup.begin(); it != m_DataObjectLookup.end(); ++it) {
    IDType indx = it->second;
    DataObject *obj1 = m_DataObjects[indx];
    if (obj1 == obj) {
      return it->first;
    }
  }
  throw std::exception();
}

DataGroup::~DataGroup()
{
    //TODO: Real Cleanup ..
}

} /* namespace DataStore */

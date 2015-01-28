/*
 * DatastoreInterface.cpp
 *
 *  Created on: Dec 2, 2014
 *      Author: settgast
 */




#include "DatastoreInterface.hpp"

namespace DataStoreNS
{

namespace
{
  std::unordered_map<std::string,DataStore> DataStores;
}


DataStore* CreateDataStore( const std::string& name )
{
  DataStores.insert(std::make_pair( name,DataStore() ) );
  return GetDataStore(name);
}

DataStore* GetDataStore( const std::string& name )
{
  return &(DataStores.at(name));
}


DataObject* CreateDataObject( DataGroup* const dataGroup, const std::string& name )
{
  dataGroup->CreateDataObject(name);
}





}

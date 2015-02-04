/*
 * DatastoreInterface.cpp
 *
 *  Created on: Dec 2, 2014
 *      Author: settgast
 */




#include "DatastoreInterface.hpp"

namespace DataStore
{

namespace
{
  std::unordered_map<std::string,DataGroup> DataStores;
}


DataGroup* CreateDataStore( const std::string& name )
{
  DataStores.insert(std::make_pair( name,DataGroup(name, std::string(), nullptr ) ) );
  return GetDataStore(name);
}

DataGroup* GetDataStore( const std::string& name )
{
  return &(DataStores.at(name));
}


DataObject* CreateDataObject( DataGroup* const dataGroup, const std::string& name )
{
  dataGroup->CreateDataObject(name);
}





}

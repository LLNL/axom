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



DataGroup::DataGroup( const std::string& name,
                      const std::string& path,
                      DataGroup* const parent,
                      DataStore* const dataStore ):
    m_parent(parent),
    m_dataStore(dataStore),
    m_DataObjects(),
    m_DataObjectLookup(),
    m_childGroups(),
    m_name(name),
    m_dataShape()
{}


DataGroup::DataGroup( const DataGroup& source ):
    m_parent(source.m_parent),
    m_dataStore(source.m_dataStore),
    m_DataObjects(source.m_DataObjects),
    m_DataObjectLookup(source.m_DataObjectLookup),
    m_childGroups(source.m_childGroups),
    m_name(source.m_name),
    m_dataShape(source.m_dataShape)
{
}

DataGroup::DataGroup( DataGroup&& source ):
    m_parent(std::move(source.m_parent)),
    m_dataStore(std::move(source.m_dataStore)),
    m_DataObjects(std::move(source.m_DataObjects)),
    m_DataObjectLookup(std::move(source.m_DataObjectLookup)),
    m_childGroups(std::move(source.m_childGroups)),
    m_name(std::move(source.m_name)),
    m_dataShape(std::move(source.m_dataShape))
{
  source.m_dataStore = nullptr;
}

DataGroup::~DataGroup()
{
  for( auto obj : m_childGroups )
  {
    delete obj.second;
  }

}








DataObject* DataGroup::CreateDataObject( const std::string& name )
{
  if( m_DataObjectLookup.count(name) == 0 )
  {
    DataObject* newObj = m_dataStore->CreateDataObject(name,this);
    m_DataObjectLookup.insert( std::make_pair(name, m_DataObjects.size() ) );
    m_DataObjects.push_back( newObj );
    m_DataObjects.back()->SetDataShape( this->GetDataShape() );
  }
  else
  {
    throw std::exception();
  }

  return m_DataObjects.back();
}

void DataGroup::RemoveDataObject( const std::string& name, const bool removeFromDataStore )
{
  auto iterLookup = m_DataObjectLookup.find(name);
  if( iterLookup != m_DataObjectLookup.end() )
  {
    const auto index = iterLookup->second;
    DataObject* dobj = m_DataObjects.at(index);
    m_DataObjects.erase(m_DataObjects.begin()+index);
    m_DataObjectLookup.erase(iterLookup);

    if( removeFromDataStore )
    {
//      m_dataStore->RemoveDataObject();
    }
  }


}


DataGroup* DataGroup::CreateDataGroup( const std::string& name )
{
  auto iter = m_childGroups.find(name);
  if( iter == m_childGroups.end() )
  {
    m_childGroups.insert( std::make_pair(name, new DataGroup(name,"", this, this->m_dataStore) ) );
    iter =  m_childGroups.find(name);
  }

  return iter->second;
}



DataGroup* DataGroup::GetDataGroup( const Attribute& Attribute )
{
  DataGroup* rGroup = CreateDataGroup("att");
  for( auto i : m_DataObjects )
  {
    if( i->HasAttribute( Attribute.Name() ) )
    {
      rGroup->m_DataObjects.push_back(i);
    }
  }
  return rGroup;
}


} /* namespace DataStore */

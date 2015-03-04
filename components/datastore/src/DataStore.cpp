/*
 * DataGroup.cpp
 *
 *  Created on: Dec 2, 2014
 *      Author: settgast
 */

#include "DataStore.hpp"
#include "DataObject.hpp"
#include "DataGroup.hpp"

namespace DataStoreNS
{



DataStore::DataStore():
    m_DataObjects(),
//    m_DataObjectLookup(),
    m_childGroups(),
    m_DataGroupLookup()

{
}



DataStore::DataStore( const DataStore& source ):
    m_DataObjects(source.m_DataObjects),
//    m_DataObjectLookup(source.m_DataObjectLookup),
    m_childGroups(source.m_childGroups),
    m_DataGroupLookup(source.m_DataGroupLookup)
{
}

DataStore::DataStore( DataStore&& source ):
    m_DataObjects(std::move(source.m_DataObjects)),
//    m_DataObjectLookup(std::move(source.m_DataObjectLookup)),
    m_childGroups(std::move(source.m_childGroups)),
    m_DataGroupLookup(std::move(source.m_DataGroupLookup))
{}

DataStore::~DataStore()
{
  for( auto obj : m_DataObjects )
  {
    delete obj;
  }

  for( auto obj : m_childGroups )
  {
    delete obj.second;
  }

}








DataObject* DataStore::CreateDataObject( const std::string& name,
                                         const DataGroup* const parent )
{
//  if( m_DataObjectLookup.count(name) == 0 )
  {
//    m_DataObjectLookup.insert( std::make_pair(name, m_DataObjects.size() ) );
    m_DataObjects.push_back(  new DataObject( name, parent ) );
  }
/*  else
  {
    throw std::exception();
  }
*/
  return m_DataObjects.back();
}



void DataStore::RemoveDataObject( const std::size_t& index )
{
  DataObject* const dobj = m_DataObjects[index];
  for( auto group : dobj->GetGroups() )
  {
    group.second->RemoveDataObject( dobj->Name(), false );
  }

  m_DataObjects.erase(m_DataObjects.begin()+index);
}



void DataStore::RemoveDataObject( DataObject*& target_dobj )
{
  DataObject* dobj = nullptr;
  std::size_t index = -1;
  for( index = 0 ; index < m_DataObjects.size() ; ++index )
  {
    dobj = m_DataObjects[index];
    if( target_dobj == dobj )
      break;
  }

  for( auto group : target_dobj->GetGroups() )
  {
    group.second->RemoveDataObject( target_dobj->Name(), false );
  }

  m_DataObjects.erase(m_DataObjects.begin()+index);
  target_dobj = nullptr;
}



DataGroup* DataStore::CreateDataGroup( const std::string& name )
{
  auto iter = m_childGroups.find(name);
  if( iter == m_childGroups.end() )
  {
    m_childGroups.insert( std::make_pair(name, new DataGroup(name,"/", nullptr, this) ) );
    iter =  m_childGroups.find(name);
  }

  return iter->second;
}

/*
DataGroup* DataStore::GetDataGroup( const Attribute& Attribute )
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
*/



} /* namespace DataStoreNS */

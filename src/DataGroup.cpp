/*
 * DataGroup.cpp
 *
 *  Created on: Dec 2, 2014
 *      Author: settgast
 */

#include "DataGroup.hpp"

namespace DataStore
{




DataGroup::DataGroup( const std::string& name,
                      const std::string& path,
                      DataGroup* const parent ):
    DataObject(name,path, parent )
{
}


DataGroup::DataGroup( const DataGroup& source ):
    DataObject(source)
{

}

DataGroup::DataGroup( const DataGroup&& source ):
    DataObject(std::move(source))
{

}

DataGroup::~DataGroup()
{
  for( auto obj : m_DataObjects )
  {
    if( obj->GetParent() == this )
    {
      delete obj;
    }
  }

  for( auto obj : m_childGroups )
  {
    delete obj.second;
  }

}








DataObject* DataGroup::CreateDataObject( const std::string& name )
{
  if( m_DataObjectLookup.count(name) == 0 )
  {
    m_DataObjectLookup.insert( std::make_pair(name, m_DataObjects.size() ) );

//    m_DataObjects.push_back( std::move(std::unique_ptr< DataObject >( new DataObject( name,Path() ) ) ) );
    m_DataObjects.push_back(  new DataObject( name,Path(), this ) );
  }
  else
  {
    throw std::exception();
  }

  return m_DataObjects.back();
}



DataGroup* DataGroup::CreateDataGroup( const std::string& name )
{
  auto iter = m_childGroups.find(name);
  if( iter == m_childGroups.end() )
  {
    m_childGroups.insert( std::make_pair(name, new DataGroup(name,"", this) ) );
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

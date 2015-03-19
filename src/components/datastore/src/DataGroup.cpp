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

  bool DataGroup::HasName( const std::string& name )
  {
    // XXX must perform two lookups: Views and groups
    DataGroup::lookupType::iterator it = m_DataViewLookup.find( name );
    if( it != m_DataViewLookup.end() )
      return true;
    DataGroup::lookupGroup::iterator itg = m_childGroups.find( name );
    if( itg != m_childGroups.end() )
      return true;
    return false;
  }

  /*
  DataView *DataGroup::AddDataView( const std::string& name, DataView *obj )
  {
    if( HasName( name ) )
    {
      throw std::exception();
    }
    if( obj == nullptr )
    {
      obj = m_datastore->CreateDataBuffer();
    }
    m_DataViewLookup[name] = m_DataViews.size(); // map name to index
    m_DataViews.push_back( obj );
    obj->AttachGroup( this );
    // XXX how does user get index? Search IndexDataView(name) or return from here?
    return obj;
  }
*/
  DataView *DataGroup::RemoveDataView( const std::string& name )
  {
    DataGroup::lookupType::iterator it = m_DataViewLookup.find( name );
    if( it != m_DataViewLookup.end() )
    {
      IDType indx = it->second;
      DataView *obj = m_DataViews[indx];
      m_DataViewLookup.erase( it );
      m_DataViews[indx] = nullptr; // XXX remove from m_DataViews
      // XXX this makes a hole in the table, but perserved existing indexes.
    }
    else
    {
      throw std::exception();
    }
  }

  DataView *DataGroup::RemoveDataView( DataView *obj )
  {
    std::string const & name = NameDataView( obj );
    RemoveDataView( name );
  }

  DataGroup *DataGroup::CreateDataGroup( const std::string& name )
  {
    if( HasName( name ) )
    {
      throw std::exception();
    }
    DataGroup *grp = new DataGroup( this, this->m_datastore );
    m_childGroups[name] = grp;
    return grp;
  }

  void DataGroup::ClearDataViews()
  {
    m_DataViews.clear();
    m_DataViewLookup.clear();
  }

  std::string const & DataGroup::NameDataView( DataView *obj )
  {
    DataGroup::lookupType::iterator it;

    // XXX brute force search for obj
    for( it = m_DataViewLookup.begin(); it != m_DataViewLookup.end() ; ++it )
    {
      IDType indx = it->second;
      DataView *obj1 = m_DataViews[indx];
      if( obj1 == obj )
      {
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

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


  DataView *DataGroup::AttachDataView( const std::string& name, DataView * const obj )
  {
    if( HasName( name ) || obj==nullptr )
    {
      throw std::exception();
    }

    m_DataViewLookup[name] = m_DataViews.size(); // map name to index
    m_DataViews.push_back( obj );
    // XXX how does user get index? Search IndexDataView(name) or return from here?
    return obj;
  }

  DataView* DataGroup::CreateDataView( const std::string& name )
  {
    DataView* const view = new DataView( name, this, m_datastore );
    return AttachDataView(name, NULL);
  }


  DataView* DataGroup::DetatchDataView( const std::string& name )
  {
    DataView* view = nullptr;
    DataGroup::lookupType::iterator it = m_DataViewLookup.find( name );
    if( it != m_DataViewLookup.end() )
    {
      IDType id = it->second;
      view = m_DataViews[id];
      m_DataViewLookup.erase( it );
      m_DataViews[id] = nullptr;
    }
    else
    {
      throw std::exception();
    }

    return view;
  }


  void DataGroup::RemoveDataView( const std::string& name )
  {
    delete DetatchDataView(name);
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

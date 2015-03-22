/*
 * DataGroup.cpp
 *
 *  Created on: Dec 2, 2014
 *      Author: settgast
 */

#include "DataGroup.hpp"
#include "DataStore.hpp"
#include "DataView.hpp"

namespace DataStoreNS
{
    DataGroup::DataGroup(const std::string &name,
                         DataGroup *parent)
    :m_name(name),
     m_parent(parent),
     m_datastore(parent->GetDataStore())
    {}

    DataGroup::DataGroup(const std::string &name,
                         DataStore *datastore)
    :m_name(name),
     m_parent(datastore->GetRoot()),
     m_datastore(datastore)
    {}

    DataGroup::~DataGroup()
    {
        DestroyViews();
        DestroyGroups();
    }


    bool DataGroup::HasChild( const std::string& name )
    {
        return HasView(name) || HasGroup(name);
    }


    bool DataGroup::HasGroup( const std::string& name )
    {
        std::map<std::string,IDType>::iterator itr;
        itr = m_groupsNameMap.find( name );
        if( itr != m_groupsNameMap.end() )
            return true;
        return false;
    }

    bool DataGroup::HasView( const std::string& name )
    {
        std::map<std::string,IDType>::iterator itr;
        itr = m_viewsNameMap.find( name );
        if( itr != m_viewsNameMap.end() )
            return true;
        return false;
    }
    /// --- DataView Children --- ///

    DataView* DataGroup::CreateView( const std::string& name )
    {
        DataView* const view = new DataView( name, this, m_datastore );
        return AttachView(view);
    }

    DataView* DataGroup::CreateView( const std::string& name,
                                     DataBuffer *buff)
    {
        DataView* const view = new DataView( name, this, buff );
        return AttachView(view);
    }


    DataView *DataGroup::AttachView(DataView * const view)
    {
        if( HasChild( view->GetName()) || view==nullptr )
        {
          throw std::exception();
        }

        m_viewsNameMap[view->GetName()] = m_views.size(); // map name to index
        m_views.push_back( view );

//        for( std::map<std::string,IDType>::iterator iter=m_viewsNameMap.begin() ; iter!=m_viewsNameMap.end() ; ++iter )
//        { std::cout<<iter->first<<std::endl; }

        return view;
    }


    DataView* DataGroup::DetachView(const std::string& name )
    {
          DataView* view = nullptr;
          std::map<std::string,IDType>::iterator itr;
          itr = m_viewsNameMap.find( name );
          if( itr != m_viewsNameMap.end() )
          {
                IDType idx = itr->second;
                view = m_views[idx];
                m_viewsNameMap.erase( itr );
                m_views[idx] = nullptr; // remove?
          }
          else
          {
              throw std::exception();
          }
          return view;
    }

    DataView* DataGroup::DetachView(IDType idx)
    {
        DataView *view = m_views[idx];
        std::map<std::string,IDType>::iterator itr;
        itr = m_viewsNameMap.find(view->GetName());
        m_viewsNameMap.erase( itr );
        m_views[idx] = nullptr; // remove?
        return view;
    }


    DataView* DataGroup::DetachView(DataView *view)
    {
        /// TODO !
        return nullptr;
    }

    // remove vs destroy vs delete
    void DataGroup::DestroyView( const std::string& name )
    {
        delete DetachView(name);
    }

    void DataGroup::DestroyView( IDType idx )
    {
        delete DetachView(idx);
    }

    void DataGroup::DestroyView( DataView *view )
    {
        delete DetachView(view);
    }
    
    /// --- DataGroup Children --- ///

    DataGroup* DataGroup::CreateGroup( const std::string& name )
    {
        DataGroup*  grp = new DataGroup( name, this);
        return AttachGroup(grp);
    }

    DataGroup *DataGroup::AttachGroup(DataGroup * const grp)
    {
        if( HasChild( grp->GetName()) || grp==nullptr )
        {
          throw std::exception();
        }

        m_groupsNameMap[grp->GetName()] = m_groups.size(); // map name to index
        m_groups.push_back( grp );
        return grp;
    }


    DataGroup* DataGroup::DetachGroup(const std::string& name )
    {
          DataGroup* grp = nullptr;
          std::map<std::string,IDType>::iterator itr;
          itr = m_groupsNameMap.find( name );
          if( itr != m_groupsNameMap.end() )
          {
                IDType idx = itr->second;
                grp = m_groups[idx];
                m_groupsNameMap.erase( itr );
                m_groups[idx] = nullptr; // remove?
          }
          else
          {
              throw std::exception();
          }
          return grp;
    }

    DataGroup* DataGroup::DetachGroup(IDType idx)
    {
        DataGroup *grp = m_groups[idx];
        std::map<std::string,IDType>::iterator itr;
        itr = m_groupsNameMap.find(grp->GetName());
        m_groupsNameMap.erase( itr );
        m_groups[idx] = nullptr; // remove?
        return grp;
    }


    DataGroup* DataGroup::DetachGroup(DataGroup *grp)
    {
        /// TODO !
        return nullptr;
    }


    // remove vs destroy vs delete
    void DataGroup::DestroyGroup( const std::string& name )
    {
        delete DetachGroup(name);
    }

    void DataGroup::DestroyGroup( IDType idx )
    {
        delete DetachGroup(idx);
    }

    void DataGroup::DestroyGroup( DataGroup *view )
    {
        delete DetachGroup(view);
    }
    
    // real cleanup
    void DataGroup::DestroyGroups()
    {
    }
    
    
    void DataGroup::DestroyViews()
    {
    }

    void DataGroup::Print() const
    {
        Node n;
        Print(n);
        n.print();
    }

    
    void DataGroup::Print(Node &n) const
    {
        n["DataGroup/name"] = m_name;
        for(IDType i=0;i<this->CountViews();i++)
        {
            DataView const *view = this->GetView(i);
            Node &v = n["DataGroup/views"].fetch(view->GetName());
            view->Print(v);

        }
        for(IDType i=0;i<this->CountGroups();i++)
        {
            DataGroup const *grp =  this->GetGroup(i);
            Node &g = n["DataGroup/groups"].fetch(grp->GetName());
            grp->Print(g);
        }
    }
    
    
    void DataGroup::PrintTree( const int level ) const
    {
      for( int i=0 ; i<level ; ++i ) std::cout<<"    ";
      std::cout<<"DataGroup "<<this->GetName()<<std::endl;

      for( std::map<std::string,IDType>::const_iterator viewIter=m_viewsNameMap.begin() ;
           viewIter!=m_viewsNameMap.end() ;
           ++viewIter )
      {
        for( int i=0 ; i<level+1 ; ++i ) std::cout<<"    ";
        std::cout<<"DataView "<<viewIter->first<<std::endl;
      }


      for( std::map<std::string,IDType>::const_iterator groupIter=m_groupsNameMap.begin() ;
           groupIter!=m_groupsNameMap.end() ;
           ++groupIter )
      {
        IDType index = groupIter->second;
        m_groups[index]->PrintTree( level + 1 );
      }

    }


} /* namespace DataStore */

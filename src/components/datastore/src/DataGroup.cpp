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
        if(HasChild(name))
        {
            throw std::exception();
        }
        DataView* const view = new DataView( name, this);
        return AttachView(view);
    }

    DataView* DataGroup::CreateView( const std::string& name,
                                     DataBuffer *buff)
    {
        if(HasChild(name))
        {
            throw std::exception();
        }
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
          IDType idx;
          itr = m_viewsNameMap.find( name );
          if( itr != m_viewsNameMap.end() )
          {
                idx = itr->second;
                view = m_views[idx];
                m_viewsNameMap.erase( itr );
                m_views.erase(m_views.begin() + idx);
          }
          else
          {
              throw std::exception();
          }
          
          // any entry in m_viewsNameMap above idx needs to shift down by 1
          for(itr = m_viewsNameMap.begin();itr!= m_viewsNameMap.end();itr++)
          {
              if(itr->second > idx)
              {
                  itr->second--;
              }
          }
          return view;
    }

    DataView* DataGroup::DetachView(IDType idx)
    {
        DataView *view = m_views[idx];
        std::map<std::string,IDType>::iterator itr;
        itr = m_viewsNameMap.find(view->GetName());
        m_viewsNameMap.erase( itr );
        m_views.erase(m_views.begin() + idx);
        // any entry in m_viewsNameMap above idx needs to shift down by 1
        for(itr = m_viewsNameMap.begin();itr!= m_viewsNameMap.end();itr++)
        {
            if(itr->second > idx)
            {
                itr->second--;
            }
        }

        return view;
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
          IDType idx;
          if( itr != m_groupsNameMap.end() )
          {
                idx = itr->second;
                grp = m_groups[idx];
                m_groupsNameMap.erase( itr );
                m_groups[idx] = nullptr; // remove?
          }
          else
          {
              throw std::exception();
          }
          // any entry in m_groupsNameMap above idx needs to shift down by 1
          for(itr = m_groupsNameMap.begin();itr!= m_groupsNameMap.end();itr++)
          {
              if(itr->second > idx)
              {
                  itr->second--;
              }
          }
          return grp;
    }

    DataGroup* DataGroup::DetachGroup(IDType idx)
    {
        DataGroup *grp = m_groups[idx];
        std::map<std::string,IDType>::iterator itr;
        itr = m_groupsNameMap.find(grp->GetName());
        m_groupsNameMap.erase( itr );
        // any entry in m_groupsNameMap above idx needs to shift down by 1
        for(itr = m_groupsNameMap.begin();itr!= m_groupsNameMap.end();itr++)
        {
            if(itr->second > idx)
            {
                itr->second--;
            }
        }
        return grp;
    }


    void DataGroup::DestroyGroup( const std::string& name )
    {
        delete DetachGroup(name);
    }

    void DataGroup::DestroyGroup( IDType idx )
    {
        delete DetachGroup(idx);
    }

    // real cleanup
    void DataGroup::DestroyGroups()
    {
        // delete all groups
        size_t ngroups = CountGroups();
        
        for(size_t i=0;i<ngroups;i++)
        {
            DataGroup *grp = this->GetGroup(i);
            delete grp;
        }

        // clean up book keeping
        m_groups.clear();
        m_groupsNameMap.clear();
    }
    
    
    void DataGroup::DestroyViews()
    {
        // delete all views
        size_t nviews = CountViews();
        
        for(size_t i=0;i<nviews;i++)
        {
            DataView *view = this->GetView(i);
            delete view;
        }
        // clean up book keeping
        m_views.clear();
        m_viewsNameMap.clear();
        
    }

    void DataGroup::Print() const
    {
        Node n;
        Info(n);
        n.print();
    }

    
    void DataGroup::Info(Node &n) const
    {
        n["name"] = m_name;
        for(IDType i=0;i<this->CountViews();i++)
        {
            DataView const *view = this->GetView(i);
            Node &v = n["views"].fetch(view->GetName());
            view->Info(v);

        }
        for(IDType i=0;i<this->CountGroups();i++)
        {
            DataGroup const *grp =  this->GetGroup(i);
            Node &g = n["groups"].fetch(grp->GetName());
            grp->Info(g);
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

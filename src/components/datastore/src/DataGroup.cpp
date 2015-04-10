/*
 * DataGroup.cpp
 *
 *  Created on: Dec 2, 2014
 *      Author: settgast
 */

#include "DataGroup.hpp"

#include "DataStore.hpp"
#include "DataBuffer.hpp"
#include "DataView.hpp"

#include "Utilities.hpp"

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

    DataView* DataGroup::CreateViewAndBuffer( const std::string& name )
    {
        ASCTK_ASSERT_MSG( HasChild(name) == false, "name == " << name );

        DataBuffer *buff = this->GetDataStore()->CreateBuffer();
        DataView* const view = new DataView( name, this,buff);
        buff->AttachView(view);
        return AttachView(view);
    }
    
     DataView *DataGroup::CreateOpaqueView( const std::string& name,
                                            void *opaque)
     {
        ASCTK_ASSERT_MSG( HasChild(name) == false, "name == " << name );
         
        DataView* const view = new DataView(name, this,opaque);
        return AttachView(view);
     }

    DataView* DataGroup::CreateView( const std::string& name,
                                     DataBuffer *buff)
    {
        ASCTK_ASSERT_MSG( HasChild(name) == false, "name == " << name );
        ASCTK_ASSERT( buff != 0 );

        DataView* const view = new DataView( name, this, buff );
        return AttachView(view);
    }

    DataView *DataGroup::MoveView(DataView *view)
    {
        ASCTK_ASSERT( view != 0 );
        ASCTK_ASSERT_MSG( HasChild(view->GetName()) == false, \
                          "view->GetName() == " << view->GetName() );
        
        // remove this view from its current parent
        DataGroup *curr_grp = view->GetParent();
        
        curr_grp->DetachView(view->GetName());
        
        /// finally, attach to this group
        AttachView(view);
        
        return view;
    }


    // creates a copy of the given view for this group
    // Recall:copying the view does not imply copying the buffer.
    // returns the new view
    DataView *DataGroup::CopyView(DataView *view)
    {
        ASCTK_ASSERT( view != 0 );
        ASCTK_ASSERT_MSG( HasChild(view->GetName()) == false, \
                          "view->GetName() == " << view->GetName() );
        
        DataView *res = CreateView(view->GetName(),view->GetBuffer());
        res->Declare(view->GetDescriptor());
        if(view->Applied())
        {
            res->Apply();
        }
        return res;
    }

    DataView *DataGroup::AttachView(DataView * const view)
    {
        ASCTK_ASSERT( view != 0 );
        ASCTK_ASSERT_MSG( HasChild(view->GetName()) == false, \
                          "view->GetName() == " << view->GetName() );

        m_viewsNameMap[view->GetName()] = m_views.size(); // map name to index
        m_views.push_back( view );
        return view;
    }


    DataView* DataGroup::DetachView(const std::string& name )
    {
        DataView* view = nullptr;
        std::map<std::string,IDType>::iterator itr;
        IDType idx;
        itr = m_viewsNameMap.find( name );
        if ( itr == m_viewsNameMap.end() )
        {
           ASCTK_WARNING("No view with name " << name << " -- null return value"); 
        }
        else {
           idx = itr->second;
           view = m_views[idx];
           m_viewsNameMap.erase( itr );
           m_views.erase(m_views.begin() + idx);
        
          // any entry in m_viewsNameMap above idx needs to shift down by 1
          for (itr = m_viewsNameMap.begin();itr!= m_viewsNameMap.end();itr++)
          {
             if(itr->second > idx)
             {
                itr->second--;
             }
          }
        
          view->m_group = nullptr;
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
        view->m_group = nullptr;
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
    
    // remove vs destroy vs delete
    void DataGroup::DestroyViewAndBuffer( const std::string& name )
    {
        DataView* view = DetachView(name);
        DataBuffer * const buffer = view->GetBuffer();
        delete view;

        // there should be a better way?
        GetDataStore()->DestroyBuffer(buffer->GetUID());

    }

    void DataGroup::DestroyViewAndBuffer( IDType idx )
    {
        DataView* view = DetachView(idx);
        // there should be a better way?
        GetDataStore()->DestroyBuffer(view->GetBuffer()->GetUID());
        delete view;
    }
    
    
    /// --- DataGroup Children --- ///

    DataGroup* DataGroup::CreateGroup( const std::string& name )
    {
        DataGroup*  grp = new DataGroup( name, this);
        return AttachGroup(grp);
    }


    DataGroup *DataGroup::MoveGroup(DataGroup *grp)
    {
        ASCTK_ASSERT( grp != 0 );
        ASCTK_ASSERT_MSG( HasChild(grp->GetName()) == false, \
                          "grp->GetName() == " << grp->GetName() );
        
        // remove this grp from its current parent
        DataGroup *curr_grp = grp->GetParent();
        
        curr_grp->DetachGroup(grp->GetName());
        
        /// finally, attach to this group
        AttachGroup(grp);

        return grp;
    }


    // creates a copy of the given group for this group
    // Recall:copying the views does not imply copying the buffers.
    // returns the new group
    DataGroup *DataGroup::CopyGroup(DataGroup *grp)
    {
        ASCTK_ASSERT( grp != 0 );
        ASCTK_ASSERT_MSG( HasChild(grp->GetName()) == false, \
                          "grp->GetName() == " << grp->GetName() );
        
        DataGroup *res = CreateGroup(grp->GetName());
    
        // copy all groups
        size_t nchild_grps = grp->CountGroups();
        for(size_t i=0; i < nchild_grps; i++)
        {
            res->CopyGroup(grp->GetGroup(i));
        }

    
        size_t nchild_views = grp->CountViews();
        for(size_t i=0; i < nchild_views; i++)
        {
            res->CopyView(grp->GetView(i));
        }

        return res;
    }

    DataGroup *DataGroup::AttachGroup(DataGroup * const grp)
    {
        ASCTK_ASSERT( grp != 0 );
        ASCTK_ASSERT_MSG( HasChild(grp->GetName()) == false, \
                          "grp->GetName() == " << grp->GetName() );

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
       if ( itr == m_groupsNameMap.end() )
       {
          ASCTK_WARNING("No view with name " << name << " -- null return value"); 
       }
       else
       {
          idx = itr->second;
          grp = m_groups[idx];
          m_groupsNameMap.erase( itr );
          m_groups.erase(m_groups.begin() + idx);

          // any entry in m_groupsNameMap above idx needs to shift down by 1
          for(itr = m_groupsNameMap.begin();itr!= m_groupsNameMap.end();itr++)
          {
              if(itr->second > idx)
              {
                  itr->second--;
              }
          }
          grp->m_parent = nullptr;

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
        grp->m_parent = nullptr;
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

    void DataGroup::DestroyViewsAndBuffers()
    {
        // delete all views
        size_t nviews = CountViews();
        
        for(size_t i=0;i<nviews;i++)
        {
            DataView *view = this->GetView(i);
            GetDataStore()->DestroyBuffer(view->GetBuffer()->GetUID());
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

/*!
 ******************************************************************************
 *
 * \file
 *
 * \brief   Implementation file for DataGroup class.
 *
 ******************************************************************************
 */

// Associated header file
#include "DataGroup.hpp"

// SiDRe project headers
#include "DataBuffer.hpp"
#include "DataStore.hpp"
#include "DataView.hpp"
#include "Utilities.hpp"

using conduit::NodeIterator;


namespace sidre
{

////////////////////////////////////////////////////////////////////////
//
//  Methods for managing DataView objects in DataGroup 
//
////////////////////////////////////////////////////////////////////////

/*
*************************************************************************
*
* Return true if DataGroup owns a DataView with given name; else false.
*
*************************************************************************
*/
bool DataGroup::hasView( const std::string& name ) const
{
    std::map<std::string, IDType>::const_iterator itr;
    itr = m_viewsNameMap.find( name );
    if( itr != m_viewsNameMap.end() )
        return true;
    return false;
}

/*
*************************************************************************
*
* Return true if DataGroup owns a DataView with given index; else false.
*
*************************************************************************
*/
bool DataGroup::hasView( const IDType idx ) const
{
   return ( idx >= 0 && idx < m_views.size() );
}

/*
*************************************************************************
*
* Create view and buffer with given name and attach to group.
*
*************************************************************************
*/
DataView* DataGroup::createViewAndBuffer( const std::string& name )
{
    ATK_ASSERT_MSG( hasView(name) == false, "name == " << name );

    DataBuffer *buff = this->getDataStore()->createBuffer();
    DataView* const view = new DataView( name, this,buff);
    buff->attachView(view);

    return attachView(view);
}

/*
*************************************************************************
*
* Create opaque view and attach to group.
*
*************************************************************************
*/
DataView *DataGroup::createOpaqueView( const std::string& name,
                                       void* opaque_ptr)
{
   ATK_ASSERT_MSG( hasView(name) == false, "name == " << name );
    
   DataView* const view = new DataView(name, this, opaque_ptr);
   return attachView(view);
}

/*
*************************************************************************
*
* Create view associated with given buffer and attach to group.
*
*************************************************************************
*/
DataView* DataGroup::createView( const std::string& name,
                                 DataBuffer* buff)
{
    ATK_ASSERT_MSG( hasView(name) == false, "name == " << name );
    ATK_ASSERT( buff != 0 );

    DataView* const view = new DataView( name, this, buff );
    return attachView(view);
}

/*
*************************************************************************
*
* Detach view with given name from group and destroy view.
*
* Data buffer remains intact in DataStore.
*
*************************************************************************
*/
void DataGroup::destroyView( const std::string& name )
{
    delete detachView(name);
}

/*
*************************************************************************
*
* Detach view with given index from group and destroy view.
*
* Data buffer remains intact in DataStore.
*
*************************************************************************
*/
void DataGroup::destroyView( IDType idx )
{
    delete detachView(idx);
}

/*
*************************************************************************
*
* Detach all views from group and destroy them.
*
* Data buffers remain intact in DataStore.
*
*************************************************************************
*/
void DataGroup::destroyViews()
{
    size_t nviews = getNumberOfViews();
                     
    for (size_t i=0;i<nviews;i++)
    {
        DataView *view = this->getView(i);
        delete view;
    }

    // clean up book keeping
    m_views.clear(); 
    m_viewsNameMap.clear();
}

/*
*************************************************************************
*
* Detach view with given name from group and destroy view. 
*
* DataBuffer in DataStore is destroyed.
*
*************************************************************************
*/
void DataGroup::destroyViewAndBuffer( const std::string& name )
{
    DataView* view = detachView(name);
    DataBuffer * const buffer = view->getBuffer();
    delete view;

    // there should be a better way?
    getDataStore()->destroyBuffer(buffer->getUID());

}

/*
*************************************************************************
*
* Detach view with given index from group and destroy view.
*
* DataBuffer in DataStore is destroyed.
*
*************************************************************************
*/
void DataGroup::destroyViewAndBuffer( IDType idx )
{
    DataView* view = detachView(idx);
    // there should be a better way?
    getDataStore()->destroyBuffer(view->getBuffer()->getUID());
    delete view;
}

/*
*************************************************************************
*
* Detach all views from group and destroy them.
*
* DataBuffers in DataStore are destroyed.
*
*************************************************************************
*/
void DataGroup::destroyViewsAndBuffers()
{
    size_t nviews = getNumberOfViews();

    for (size_t i=0;i<nviews;i++)
    {
        DataView *view = this->getView(i);
        getDataStore()->destroyBuffer(view->getBuffer()->getUID());
        delete view;
    }

    // clean up book keeping
    m_views.clear();
    m_viewsNameMap.clear();

}

/*
*************************************************************************
*
* Remove given view from its owning group and attach to this group.
*
*************************************************************************
*/
DataView* DataGroup::moveView(DataView* view)
{
    ATK_ASSERT( view != 0 );
    ATK_ASSERT_MSG( hasView(view->getName()) == false, \
                    "view->GetName() == " << view->getName() );

    // remove this view from its current parent
    DataGroup *curr_grp = view->getOwningGroup();

    curr_grp->detachView(view->getName());

    /// finally, attach to this group
    attachView(view);

    return view;
}

/*
*************************************************************************
*
* Create a copy of given view and attach to this group.
*
* Copying a view does not perform a deep copy of its data buffer. 
*
*************************************************************************
*/
DataView* DataGroup::copyView(DataView* view)
{
    ATK_ASSERT( view != 0 );
    ATK_ASSERT_MSG( hasView(view->getName()) == false, \
                    "view->GetName() == " << view->getName() );

    DataView *res = createView(view->getName(),view->getBuffer());
    res->declare(view->getDescriptor());
    if(view->isApplied())
    {
        res->apply();
    }
    return res;
}


////////////////////////////////////////////////////////////////////////
//
//  Methods for managing child DataGroup objects in DataGroup 
//
////////////////////////////////////////////////////////////////////////

/*
*************************************************************************
*
* Return true if group has an (immediate) child group with given name;
* else false.
*
*************************************************************************
*/
bool DataGroup::hasGroup( const std::string& name ) const
{
    std::map<std::string, IDType>::const_iterator itr;
    itr = m_groupsNameMap.find( name );
    if( itr != m_groupsNameMap.end() )
        return true;
    return false;
}

/*
*************************************************************************
*
* Return true if group has an (immediate) child group with given index;
* else false.
*
*************************************************************************
*/
bool DataGroup::hasGroup( const IDType idx ) const
{
   return ( idx >= 0 && idx < m_groups.size() );
}

/*
*************************************************************************
*
* Create group with given name and make it a child of this group.
*
*************************************************************************
*/
DataGroup* DataGroup::createGroup( const std::string& name )
{
    DataGroup*  grp = new DataGroup( name, this);
    return attachGroup(grp);
}

/*
*************************************************************************
*
* Remove given group from its owning group and make it a child of this group.
*
*************************************************************************
*/
DataGroup* DataGroup::moveGroup(DataGroup* grp)
{
    ATK_ASSERT( grp != 0 );
    ATK_ASSERT_MSG( hasGroup(grp->getName()) == false, \
                    "grp->GetName() == " << grp->getName() );

    // remove this grp from its current parent
    DataGroup *curr_grp = grp->getParent();

    curr_grp->detachGroup(grp->getName());

    /// finally, attach to this group
    attachGroup(grp);

    return grp;
}

/*
*************************************************************************
*
* Create a copy of given group and make it a child of this group.
*
* Copying a group does not perform a deep copy of any of its buffers. 
*
*************************************************************************
*/
DataGroup *DataGroup::copyGroup(DataGroup *grp)
{
    ATK_ASSERT( grp != 0 );
    ATK_ASSERT_MSG( hasGroup(grp->getName()) == false, \
                    "grp->GetName() == " << grp->getName() );

    DataGroup *res = createGroup(grp->getName());

    // copy all groups
    size_t nchild_grps = grp->getNumberOfGroups();
    for(size_t i=0; i < nchild_grps; i++)
    {
        res->copyGroup(grp->getGroup(i));
    }


    size_t nchild_views = grp->getNumberOfViews();
    for(size_t i=0; i < nchild_views; i++)
    {
        res->copyView(grp->getView(i));
    }

    return res;
}

/*
*************************************************************************
*
* Detach child group with given name and destroy it.
*
*************************************************************************
*/
void DataGroup::destroyGroup( const std::string& name )
{
    delete detachGroup(name);
}

/*
*************************************************************************
*
* Detach child group with given index and destroy it.
*
*************************************************************************
*/
void DataGroup::destroyGroup( IDType idx )
{
    delete detachGroup(idx);
}

/*
*************************************************************************
*
* Detach all child groups and destroy them.
*
*************************************************************************
*/
void DataGroup::destroyGroups()
{
    size_t ngroups = getNumberOfGroups();

    for(size_t i=0;i<ngroups;i++)
    {
        DataGroup *grp = this->getGroup(i);
        delete grp;
    }

    // clean up book keeping
    m_groups.clear();
    m_groupsNameMap.clear();
}



    void DataGroup::info(Node &n) const
    {
        n["name"] = m_name;
        for(IDType i=0;i<this->getNumberOfViews();i++)
        {
            DataView const *view = this->getView(i);
            Node &v = n["views"].fetch(view->getName());
            view->info(v);

        }
        for(IDType i=0;i<this->getNumberOfGroups();i++)
        {
            DataGroup const *grp =  this->getGroup(i);
            Node &g = n["groups"].fetch(grp->getName());
            grp->info(g);
        }
    }

    void DataGroup::print() const
    {
        Node n;
        info(n);
        n.print();
    }

    void DataGroup::printTree( const int nlevels ) const
    {
      for( int i=0 ; i<nlevels ; ++i ) std::cout<<"    ";
      std::cout<<"DataGroup "<<this->getName()<<std::endl;

      for( std::map<std::string,IDType>::const_iterator viewIter=m_viewsNameMap.begin() ;
           viewIter!=m_viewsNameMap.end() ;
           ++viewIter )
      {
        for( int i=0 ; i<nlevels+1 ; ++i ) std::cout<<"    ";
        std::cout<<"DataView "<<viewIter->first<<std::endl;
      }


      for( std::map<std::string,IDType>::const_iterator groupIter=m_groupsNameMap.begin() ;
           groupIter!=m_groupsNameMap.end() ;
           ++groupIter )
      {
        IDType index = groupIter->second;
        m_groups[index]->printTree( nlevels + 1 );
      }

    }

    /// ---------------------------------------------------------------
    ///  Save + Restore Prototypes (ATK-39)
    /// ---------------------------------------------------------------
    /// saves "this", associated views and buffers to a file set.
    void DataGroup::save(const std::string &obase,
                         const std::string &protocol) const
    {
        if(protocol == "conduit")
        {
            Node n;
            copyToNode(n);
            // for debugging call: n.print();
            n.save(obase);
        }
    }

    /// restores as "this"
    void DataGroup::load(const std::string &obase,
                         const std::string &protocol)
    {
        if(protocol == "conduit")
        {
            destroyGroups();
            destroyViews();
            Node n;
            n.load(obase);
            // for debugging call: n.print();
            copyFromNode(n);
        }
    }

////////////////////////////////////////////////////////////////////////
//
// Private methods below.
//
////////////////////////////////////////////////////////////////////////



    DataGroup::DataGroup(const std::string &name,
                         DataGroup *parent)
    :m_name(name),
     m_parent(parent),
     m_datastore(parent->getDataStore())
    {}

    DataGroup::DataGroup(const std::string &name,
                         DataStore *datastore)
    :m_name(name),
     m_parent(datastore->getRoot()),
     m_datastore(datastore)
    {}

    DataGroup::~DataGroup()
    {
        destroyViews();
        destroyGroups();
    }


    /// --- DataView Children --- ///

    DataView *DataGroup::attachView(DataView * const view)
    {
        ATK_ASSERT( view != 0 );
        ATK_ASSERT_MSG( hasView(view->getName()) == false, \
                        "view->GetName() == " << view->getName() );

        m_viewsNameMap[view->getName()] = m_views.size(); // map name to index
        m_views.push_back( view );
        return view;
    }


    DataView* DataGroup::detachView(const std::string& name )
    {
        DataView* view = nullptr;
        std::map<std::string,IDType>::iterator itr;
        IDType idx;
        itr = m_viewsNameMap.find( name );
        if ( itr == m_viewsNameMap.end() )
        {
           ATK_WARNING("No view with name " << name << " -- null return value"); 
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
        
          view->m_owning_group = nullptr;
       } 
       return view;
    }

    DataView* DataGroup::detachView(IDType idx)
    {
        DataView *view = m_views[idx];
        std::map<std::string,IDType>::iterator itr;
        itr = m_viewsNameMap.find(view->getName());
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
        view->m_owning_group = nullptr;
        return view;
    }


    /// --- DataGroup Children --- ///

    DataGroup *DataGroup::attachGroup(DataGroup * const grp)
    {
        ATK_ASSERT( grp != 0 );
        ATK_ASSERT_MSG( hasGroup(grp->getName()) == false, \
                        "grp->GetName() == " << grp->getName() );

        m_groupsNameMap[grp->getName()] = m_groups.size(); // map name to index
        m_groups.push_back( grp );
        return grp;
    }


    DataGroup* DataGroup::detachGroup(const std::string& name )
    {
       DataGroup* grp = nullptr;
       std::map<std::string,IDType>::iterator itr;
       itr = m_groupsNameMap.find( name );
       IDType idx;
       if ( itr == m_groupsNameMap.end() )
       {
          ATK_WARNING("No view with name " << name << " -- null return value"); 
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

    DataGroup* DataGroup::detachGroup(IDType idx)
    {
        DataGroup *grp = m_groups[idx];
        std::map<std::string,IDType>::iterator itr;
        itr = m_groupsNameMap.find(grp->getName());
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


    void DataGroup::copyToNode(Node &n) const
    {
        std::vector<IDType> buffer_ids;
        copyToNode(n,buffer_ids);

        // save the buffers discovered by buffer_ids
        for(size_t i=0; i < buffer_ids.size(); i++)
        {
            Node &buff = n["buffers"].append();
            IDType buffer_id = buffer_ids[i];
            DataBuffer *ds_buff =  m_datastore->getBuffer(buffer_id);
            buff["id"].set(buffer_id);
            buff["descriptor"].set(ds_buff->getDescriptor().to_json());
            
            // only set our data if the buffer was initialized 
            if (ds_buff->getData() != NULL )
            {
                buff["data"].set_external(ds_buff->getNode());
            }
        }

    }

    void DataGroup::copyToNode(Node &n,
                               std::vector<IDType> &buffer_ids) const
    {
        for(IDType i=0; i < this->getNumberOfViews(); i++)
        {
            DataView const *view = this->getView(i);
            Node &n_view = n["views"].fetch(view->getName());
            n_view["descriptor"].set(view->getDescriptor().to_json());
            n_view["applied"].set(view->isApplied());
            // if we have a buffer, simply add the id to the list
            if(view->hasBuffer())
            {
                IDType buffer_id = view->getBuffer()->getUID();
                n_view["buffer_id"].set(buffer_id);
                buffer_ids.push_back(view->getBuffer()->getUID());
            }
        }
        
        for(IDType i=0; i < this->getNumberOfGroups(); i++)
        {
            DataGroup const *grp =  this->getGroup(i);
            Node &n_grp = n["groups"].fetch(grp->getName());
            grp->copyToNode(n_grp,buffer_ids);
        }
    }

    void DataGroup::copyFromNode(Node &n)
    {
         std::map<IDType,IDType> id_map;
         copyFromNode(n,id_map);
    }
    
    void DataGroup::copyFromNode(Node &n,
                                 std::map<IDType,IDType> &id_map)
    {
        /// for restore each group contains:
        /// buffers, views, and groups
        
        // create the buffers
        if(n.has_path("buffers"))
        {
            NodeIterator buffs_itr = n["buffers"].iterator();
            while(buffs_itr.has_next())
            {
                Node &n_buff = buffs_itr.next();
                IDType buffer_id = n_buff["id"].as_uint64();

                // create a new mapping and buffer if necessary
                if( id_map.find(buffer_id) == id_map.end())
                {
                    DataBuffer *ds_buff = this->getDataStore()->createBuffer();
                    // map "id" to whatever new id the data store gives us.
                    IDType buffer_ds_id = ds_buff->getUID();
                    id_map[buffer_id] = buffer_ds_id;
                    // setup the new data store buffer
                    Schema schema(n_buff["descriptor"].as_string());
                    ds_buff->declare(schema);
                    if(n_buff.has_path("data"))
                    {
                        ds_buff->allocate();
                        // copy the data from the node
                        ds_buff->getNode().update(n_buff["data"]);
                    }
                }
            }
        }

        // create the child views
        NodeIterator views_itr = n["views"].iterator();
        while(views_itr.has_next())
        {
            Node &n_view = views_itr.next();
            if(n_view.has_path("buffer_id"))
            {
                std::string view_name = views_itr.path();
                
                IDType buffer_id = n_view["buffer_id"].as_uint64();
                // get the mapped buffer id
                if( id_map.find(buffer_id) == id_map.end() )
                {
                    ATK_ERROR("Invalid buffer id mapping.");
                }
                
                buffer_id = id_map[buffer_id];
                DataBuffer *ds_buff = m_datastore->getBuffer(buffer_id);

                // create a new view with the buffer
                DataView   *ds_view = createView(view_name,ds_buff);
                // declare using the schema
                Schema schema(n_view["descriptor"].as_string());
                ds_view->declare(schema);
                // if the descriptor was applied, restore this state
                if(n_view["applied"].to_uint64() != 0)
                    ds_view->apply();
            }
            else
            {
                ATK_WARNING("DataGroup cannot restore opaque views.");
            }
        }

        // create the child groups
        NodeIterator grps_itr = n["groups"].iterator();
        while(grps_itr.has_next())
        {
            Node &n_grp = grps_itr.next();
            std::string grp_name = grps_itr.path();
            DataGroup *ds_grp = createGroup(grp_name);
            ds_grp->copyFromNode(n_grp,id_map);
        }
    }


} /* namespace sidre */

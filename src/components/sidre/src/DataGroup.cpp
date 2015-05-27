/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

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

// Other toolkit component headers
#include "common/CommonTypes.hpp"

// SiDRe project headers
#include "DataBuffer.hpp"
#include "DataStore.hpp"
#include "DataView.hpp"


namespace asctoolkit
{
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
* Create view and buffer with given name, use the data type to 
* allocate the buffer and initialize the view. Attach new view to group.
*
*************************************************************************
*/
DataView* DataGroup::createViewAndBuffer( const std::string& name,
                                          const DataType &dtype)
{
    ATK_ASSERT_MSG( hasView(name) == false, "name == " << name );

    DataBuffer *buff = this->getDataStore()->createBuffer();
    DataView* const view = new DataView( name, this,buff);
    buff->attachView(view);
    view->allocate(dtype);
    return attachView(view);
}


/*
*************************************************************************
*
* Create view and buffer with given name, use the data type to 
* allocate the buffer and initialize the view. Attach new view to group.
*
*************************************************************************
*/
DataView* DataGroup::createViewAndBuffer( const std::string& name,
                                          const Schema &schema)
{
    ATK_ASSERT_MSG( hasView(name) == false, "name == " << name );

    DataBuffer *buff = this->getDataStore()->createBuffer();
    DataView* const view = new DataView( name, this,buff);
    buff->attachView(view);
    view->allocate(schema);
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
* Create view associated with given buffer, apply given data type 
* and attach to group.
*
*************************************************************************
*/
DataView* DataGroup::createView( const std::string& name,
                                 DataBuffer* buff,
                                 const DataType &dtype)
{
    ATK_ASSERT_MSG( hasView(name) == false, "name == " << name );
    ATK_ASSERT( buff != 0 );

    DataView* const view = new DataView( name, this, buff );
    view->apply(dtype);
    return attachView(view);
}


/*
*************************************************************************
*
* Create view associated with given buffer, apply given schema 
* and attach to group.
*
*************************************************************************
*/
DataView* DataGroup::createView( const std::string& name,
                                 DataBuffer* buff,
                                 const Schema &schema)
{
    ATK_ASSERT_MSG( hasView(name) == false, "name == " << name );
    ATK_ASSERT( buff != 0 );

    DataView* const view = new DataView( name, this, buff );
    view->apply(schema);
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
* Create external view with given data type and attach to group.
*
*************************************************************************
*/
DataView *DataGroup::createExternalView( const std::string& name,
                                         void* external_data,
                                         const DataType& dtype)
{
   ATK_ASSERT_MSG( hasView(name) == false, "name == " << name );
   
   DataView* const view = new DataView(name, this, 
                                       external_data, dtype);
   return attachView(view);
}

/*
*************************************************************************
*
* Create external view with given schema and attach to group.
*
*************************************************************************
*/
DataView *DataGroup::createExternalView( const std::string& name,
                                         void* external_data,
                                         const Schema& schema)
{
   ATK_ASSERT_MSG( hasView(name) == false, "name == " << name );
  
   DataView* const view = new DataView(name, this, 
                                       external_data, schema);
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
    size_t nviews = getNumViews();
                     
    for (size_t i=0; i<nviews; ++i)
    {
        DataView *view = this->getView(i);
        delete view;
    }

    m_view_coll.removeAllItems();
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
    size_t nviews = getNumViews();

    for (size_t i=0; i<nviews; ++i)
    {
        DataView *view = this->getView(i);
        getDataStore()->destroyBuffer(view->getBuffer()->getUID());
        delete view;
    }

    m_view_coll.removeAllItems();
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
                    "view->getName() == " << view->getName() );

    // remove this view from its current parent
    DataGroup* curr_group = view->getOwningGroup();

    curr_group->detachView(view->getName());

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
                    "view->getName() == " << view->getName() );

    DataView *res = createView(view->getName(), view->getBuffer());
    res->declare(view->getSchema());
    if (view->isApplied())
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
* Create group with given name and make it a child of this group.
*
*************************************************************************
*/
DataGroup* DataGroup::createGroup( const std::string& name )
{
    DataGroup*  group = new DataGroup( name, this);
    return attachGroup(group);
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
    size_t ngroups = getNumGroups();

    for (size_t i=0; i<ngroups; ++i)
    {
        DataGroup* group = this->getGroup(i);
        delete group;
    }

    m_group_coll.removeAllItems();
}

/*
*************************************************************************
*
* Remove given group from its owning group and make it a child of this group.
*
*************************************************************************
*/
DataGroup* DataGroup::moveGroup(DataGroup* group)
{
    ATK_ASSERT( group != 0 );
    ATK_ASSERT_MSG( hasGroup(group->getName()) == false, \
                    "group->getName() == " << group->getName() );

    // remove this group from its current parent
    DataGroup *curr_group = group->getParent();

    curr_group->detachGroup(group->getName());

    /// finally, attach to this group
    attachGroup(group);

    return group;
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
DataGroup* DataGroup::copyGroup(DataGroup* group)
{
    ATK_ASSERT( group != 0 );
    ATK_ASSERT_MSG( hasGroup(group->getName()) == false, \
                    "group->getName() == " << group->getName() );

    DataGroup *res = createGroup(group->getName());

    // copy all groups
    size_t nchild_groups = group->getNumGroups();
    for (size_t i=0; i < nchild_groups; ++i)
    {
        res->copyGroup(group->getGroup(i));
    }


    size_t nchild_views = group->getNumViews();
    for (size_t i=0; i < nchild_views; ++i)
    {
        res->copyView(group->getView(i));
    }

    return res;
}

/*
*************************************************************************
*
* Copy data group description to given Conduit node.
*
*************************************************************************
*/
void DataGroup::info(Node& n) const
{
    n["name"] = m_name;
    for (size_t i=0; i < this->getNumViews(); ++i)
    {
        DataView const *view = this->getView(i);
        Node& v = n["views"].fetch(view->getName());
        view->info(v);

    }
    for (size_t i=0; i<this->getNumGroups(); ++i)
    {
        DataGroup const* group =  this->getGroup(i);
        Node& g = n["groups"].fetch(group->getName());
        group->info(g);
    }
}

/*
*************************************************************************
*
* Print JSON description of data group to stdout.
*
*************************************************************************
*/
void DataGroup::print() const
{
    Node n;
    info(n);
    n.print();
}

/*
*************************************************************************
*
* Print given number of levels of group (sub) tree starting at this 
* group to stdout.
*
*************************************************************************
*/
void DataGroup::printTree( const int nlevels ) const
{
  for ( int i=0 ; i<nlevels ; ++i ) std::cout<<"    ";
  std::cout<<"DataGroup "<<this->getName()<<std::endl;

  size_t nviews = getNumViews();
  for (size_t idx = 0; idx < nviews; ++idx)
  {
     const DataView* view = getView(idx);

     for ( int i=0 ; i<nlevels+1 ; ++i ) std::cout<<"    ";
     std::cout<< "DataView " << view->getName() << std::endl;
  } 
   
  size_t ngroups = getNumGroups();
  for (size_t idx = 0; idx < ngroups; ++idx)
  {
     const DataGroup* group = getGroup(idx);

     group->printTree( nlevels + 1 ); 
  }

}


/*
*************************************************************************
*
* Save this group object (including data views and child groups) to a 
* file set named "obase".
*
* Note: Only valid protocol is "conduit".
*
*************************************************************************
*/
void DataGroup::save(const std::string& obase,
                     const std::string& protocol) const
{
    if(protocol == "conduit")
    {
        Node n;
        copyToNode(n);
        // for debugging call: n.print();
        n.save(obase);
    }
}

/*
*************************************************************************
*
* Load data group (including data views and child groups) from a file 
* set named "obase" into this group object.
*
* Note: Only valid protocol is "conduit".
*
*************************************************************************
*/
void DataGroup::load(const std::string& obase,
                     const std::string& protocol)
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
// Private methods below
//
////////////////////////////////////////////////////////////////////////


/*
*************************************************************************
*
* PRIVATE ctor makes group with given name and make it a child of parent.
*
*************************************************************************
*/
DataGroup::DataGroup(const std::string& name,
                     DataGroup* parent) 
: m_name(name),
  m_parent(parent),
  m_datastore(parent->getDataStore())
{ 
}

/*
*************************************************************************
*
* PRIVATE ctor makes group with given name and make it a child of 
* root group in datastore.
*
*************************************************************************
*/
DataGroup::DataGroup(const std::string& name,
                     DataStore* datastore) 
: m_name(name),
  m_parent(datastore->getRoot()),
  m_datastore(datastore)
{
}

/*
*************************************************************************
*
* PRIVATE dtor destroys group and all its contents.
*
*************************************************************************
*/
DataGroup::~DataGroup()
{
   destroyViews();
   destroyGroups();
}


/*
*************************************************************************
*
* PRIVATE method to attach given view to group.
*
*************************************************************************
*/
DataView* DataGroup::attachView(DataView* view)
{
    ATK_ASSERT( view != 0 );
    ATK_ASSERT_MSG( hasView(view->getName()) == false, \
                    "view->getName() == " << view->getName() );

    m_view_coll.insertItem(view, view->getName());

    return view;
}

/*
*************************************************************************
*
* PRIVATE method to detach given with given name from group.
*
*************************************************************************
*/
DataView* DataGroup::detachView(const std::string& name )
{
   DataView* view = m_view_coll.removeItem(name);
   if (view) {
      view->m_owning_group = ATK_NULLPTR;
   }

   return view;
}

/*
*************************************************************************
*
* PRIVATE method to detach view with given index from group.
*
*************************************************************************
*/
DataView* DataGroup::detachView(IDType idx)
{
   DataView* view = m_view_coll.removeItem(idx);
   if (view) {
      view->m_owning_group = ATK_NULLPTR;
   }

   return view;
}


/*
*************************************************************************
*
* PRIVATE method to make given group a child of this group.
*
*************************************************************************
*/
DataGroup* DataGroup::attachGroup(DataGroup* group)
{
    ATK_ASSERT( group != 0 );
    ATK_ASSERT_MSG( hasGroup(group->getName()) == false, \
                    "group->getName() == " << group->getName() );

    m_group_coll.insertItem(group, group->getName());  
   
    return group;
}

/*
*************************************************************************
*
* PRIVATE method to detach child group with given name from group.
*
*************************************************************************
*/
DataGroup* DataGroup::detachGroup(const std::string& name )
{
   DataGroup* group = m_group_coll.removeItem(name);
   if (group) {
      group->m_parent = ATK_NULLPTR;
   }

   return group;
}

/*
*************************************************************************
*
* PRIVATE method to detach child group with given index from group.
*
*************************************************************************
*/
DataGroup* DataGroup::detachGroup(IDType idx)
{
   DataGroup* group = m_group_coll.removeItem(idx);
   if (group) {
      group->m_parent = ATK_NULLPTR;
   }

   return group;
}


/*
*************************************************************************
*
* PRIVATE method to copy group to given Conduit node.
*
*************************************************************************
*/
void DataGroup::copyToNode(Node& n) const
{
    std::vector<IDType> buffer_ids;
    copyToNode(n,buffer_ids);

    // save the buffers discovered by buffer_ids
    for (size_t i=0; i < buffer_ids.size(); i++)
    {
        Node& buff = n["buffers"].append();
        IDType buffer_id = buffer_ids[i];
        DataBuffer *ds_buff =  m_datastore->getBuffer(buffer_id);
        buff["id"].set(buffer_id);
        buff["schema"].set(ds_buff->getSchema().to_json());
        
        // only set our data if the buffer was initialized 
        if (ds_buff->getData() != NULL )
        {
            buff["data"].set_external(ds_buff->getNode());
        }
    }

}

/*
*************************************************************************
*
* PRIVATE method to copy from given Conduit node to this group.
*
*************************************************************************
*/
void DataGroup::copyFromNode(Node& n)
{
     std::map<IDType, IDType> id_map;
     copyFromNode(n, id_map);
}
    
/*
*************************************************************************
*
* PRIVATE method to copy from group to given Conduit node using
* given vector of ids to maintain correct association of data buffers
* to data views.
*
*************************************************************************
*/
void DataGroup::copyToNode(Node& n,
                           std::vector<IDType>& buffer_ids) const
{
    for (size_t i=0; i < this->getNumViews(); ++i)
    {
        DataView const* view = this->getView(i);
        Node& n_view = n["views"].fetch(view->getName());
        n_view["schema"].set(view->getSchema().to_json());
        n_view["is_applied"].set(view->isApplied());
        n_view["is_external"].set(view->isExternal());

        // if we have a buffer, simply add the id to the list
        if (view->hasBuffer())
        {
            IDType buffer_id = view->getBuffer()->getUID();
            n_view["buffer_id"].set(buffer_id);
            buffer_ids.push_back(view->getBuffer()->getUID());
        }
    }
    
    for (size_t i=0; i < this->getNumGroups(); ++i)
    {
        DataGroup const* group =  this->getGroup(i);
        Node& n_group = n["groups"].fetch(group->getName());
        group->copyToNode(n_group, buffer_ids);
    }
}

/*
*************************************************************************
*
* PRIVATE method to copy from given Conduit node to this group using
* given map of ids to indicate association of buffer ids in node to
* those in datastore.
*
*************************************************************************
*/
void DataGroup::copyFromNode(Node& n,
                             std::map<IDType, IDType>& id_map)
{
    /// for restore each group contains:
    /// buffers, views, and groups
    
    // create the buffers
    if (n.has_path("buffers"))
    {
        conduit::NodeIterator buffs_itr = n["buffers"].iterator();
        while (buffs_itr.has_next())
        {
            Node& n_buff = buffs_itr.next();
            IDType buffer_id = n_buff["id"].as_int32();

            // create a new mapping and buffer if necessary
            if ( id_map.find(buffer_id) == id_map.end())
            {
                DataBuffer* ds_buff = this->getDataStore()->createBuffer();
                // map "id" to whatever new id the data store gives us.
                IDType buffer_ds_id = ds_buff->getUID();
                id_map[buffer_id] = buffer_ds_id;
                // setup the new data store buffer
                Schema schema(n_buff["schema"].as_string());
                ds_buff->declare(schema);
                if (n_buff.has_path("data"))
                {
                    ds_buff->allocate();
                    // copy the data from the node
                    ds_buff->getNode().update(n_buff["data"]);
                }
            }
        }
    }

    // create the child views
    conduit::NodeIterator views_itr = n["views"].iterator();
    while (views_itr.has_next())
    {
        Node& n_view = views_itr.next();
        if (n_view.has_path("buffer_id"))
        {
            std::string view_name = views_itr.path();
            
            IDType buffer_id = n_view["buffer_id"].as_int32();
            // get the mapped buffer id
            if( id_map.find(buffer_id) == id_map.end() )
            {
                ATK_ERROR("Invalid buffer id mapping.");
            }
            
            buffer_id = id_map[buffer_id];
            DataBuffer* ds_buff = m_datastore->getBuffer(buffer_id);

            // create a new view with the buffer
            DataView* ds_view = createView(view_name,ds_buff);
            // declare using the schema
            Schema schema(n_view["schema"].as_string());
            ds_view->declare(schema);
            // if the schema was applied, restore this state
            if (n_view["is_applied"].to_uint64() != 0)
                ds_view->apply();
        }
        else
        {
            ATK_WARNING("DataGroup cannot restore opaque views.");
        }
    }

    // create the child groups
    conduit::NodeIterator groups_itr = n["groups"].iterator();
    while (groups_itr.has_next())
    {
        Node& n_group = groups_itr.next();
        std::string group_name = groups_itr.path();
        DataGroup* ds_group = createGroup(group_name);
        ds_group->copyFromNode(n_group, id_map);
    }
}


} /* end namespace sidre */
} /* end namespace asctoolkit */

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
//  Methods for manipulating DataView objects in DataGroup
//
////////////////////////////////////////////////////////////////////////

/*
 *************************************************************************
 *
 * Create view with given name, data type, and number of elements. 
 * Attach view to this group. Create associated buffer, allocate the 
 * data, and attach view to new buffer. 
 *
 *************************************************************************
 */
DataView * DataGroup::createViewAndAllocate( const std::string& name,
                                             TypeID type, 
                                             SidreLength num_elems )
{
  SLIC_ASSERT( !name.empty() );
  SLIC_ASSERT_MSG( hasView(name) == false, "name == " << name );
  SLIC_ASSERT(num_elems >= 0);

  if ( name.empty() || hasView(name) || num_elems < 0 ) 
  {
    return ATK_NULLPTR;
  }
  else 
  {
    DataView * const view = createViewAndBuffer(name);
    view->allocate(type, num_elems);
    return view;
  }
}

/*
 *************************************************************************
 *
 * Create view with given name and Conduit DataType. Attach view to
 * this group. Create associated buffer, allocate the data, and attach
 * view to new buffer.
 *
 *************************************************************************
 */
DataView * DataGroup::createViewAndAllocate( const std::string& name,
                                             const DataType& dtype)
{
  SLIC_ASSERT( !name.empty() );
  SLIC_ASSERT_MSG( hasView(name) == false, "name == " << name );

  if ( name.empty() || hasView(name) ) 
  {
    return ATK_NULLPTR;
  }
  else 
  {
    DataView * const view = createViewAndBuffer(name);
    view->allocate(dtype);
    return view;
  }
}

/*
 *************************************************************************
 *
 * Create view with given name and Conduit Schema. Attach view to
 * this group. Create associated buffer, allocate the data, and attach
 * view to new buffer.
 *
 *************************************************************************
 */
DataView * DataGroup::createViewAndAllocate( const std::string& name,
                                             const Schema& schema)
{
  SLIC_ASSERT( !name.empty() );
  SLIC_ASSERT_MSG( hasView(name) == false, "name == " << name );

  if ( name.empty() || hasView(name) ) 
  {
    return ATK_NULLPTR;
  }
  else 
  {
    DataView * const view = createViewAndBuffer(name);
    view->allocate(schema);
    return view;
  }
}

/*
 *************************************************************************
 *
 * Create empty view and attach view to group.
 *
 *************************************************************************
 */
DataView * DataGroup::createView( const std::string& name )
{
  SLIC_ASSERT( !name.empty() );
  SLIC_ASSERT_MSG( hasView(name) == false, "name == " << name );

  if ( name.empty() || hasView(name) )
  {
    return ATK_NULLPTR;
  }
  else
  {
    DataView * const view = new DataView( name, this);
    return attachView(view);
  }
}

/*
 *************************************************************************
 *
 * Create view with declaration and attach it to group.
 *
 *************************************************************************
 */
DataView * DataGroup::createView( const std::string& name,
                                  TypeID type,
                                  SidreLength num_elems )
{
  SLIC_ASSERT( !name.empty() );
  SLIC_ASSERT_MSG( hasView(name) == false, "name == " << name );
  SLIC_ASSERT_MSG(num_elems >= 0, "Must define view with number of elems >=0 ");

  if ( name.empty() || hasView(name) || num_elems < 0 )
  {
    return ATK_NULLPTR;
  } 
  else
  {
    DataView * const view = new DataView( name, this);
    view->declare(type, num_elems);
    return attachView(view);
  }
}

/*
 *************************************************************************
 *
 * Create view with declaration and attach it to group.
 *
 *************************************************************************
 */
DataView * DataGroup::createView( const std::string& name,
                                  const DataType& dtype )
{
  SLIC_ASSERT( !name.empty() );
  SLIC_ASSERT_MSG( hasView(name) == false, "name == " << name );

  if ( name.empty() || hasView(name) )
  {
    return ATK_NULLPTR;
  }
  else
  {
    DataView * const view = new DataView( name, this);
    view->declare(dtype);
    return attachView(view);
  }
}

/*
 *************************************************************************
 *
 * Create view with declaration and attach it to group.
 *
 *************************************************************************
 */
DataView * DataGroup::createView( const std::string& name,
                                  const Schema& schema )
{
  SLIC_ASSERT( !name.empty() );   
  SLIC_ASSERT_MSG( hasView(name) == false, "name == " << name );
  
  if ( name.empty() || hasView(name) )
  {
    return ATK_NULLPTR;
  }
  else
  {
    DataView * const view = new DataView( name, this);
    view->declare(schema);
    return attachView(view);
  }
}


/*
 *************************************************************************
 *
 * Create view into given buffer and attach view to group.
 *
 *************************************************************************
 */
DataView * DataGroup::createView( const std::string& name,
                                  DataBuffer * buff)
{
  SLIC_ASSERT( !name.empty() );
  SLIC_ASSERT_MSG( hasView(name) == false, "name == " << name );
  SLIC_CHECK( buff != ATK_NULLPTR );

  if ( name.empty() || hasView(name) )
  {
    return ATK_NULLPTR;
  }
  else 
  {
    DataView * view = new DataView( name, this );
    if ( buff != ATK_NULLPTR )
    {
       view->attachBuffer( buff );
    }
    return attachView(view);
  }
}

/*
 *************************************************************************
 *
 * Create external data view and attach to group.
 *
 *************************************************************************
 */
DataView * DataGroup::createView( const std::string& name,
                                  void * external_ptr )
{
  SLIC_ASSERT( !name.empty() );
  SLIC_ASSERT_MSG( hasView(name) == false, "name == " << name );
  SLIC_ASSERT_MSG( external_ptr != ATK_NULLPTR ,
                   "Cannot create external view with null data pointer" );
  
  if ( name.empty() || hasView(name) || external_ptr == ATK_NULLPTR )
  {
    return ATK_NULLPTR;
  }
  else
  {
    DataView * view = new DataView(name, this);
    view->setExternalDataPtr(external_ptr);

    return attachView(view);
  }
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
  SLIC_CHECK_MSG( hasView(name) == true, "name == " << name );

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
void DataGroup::destroyView( IndexType idx )
{
  SLIC_CHECK_MSG( hasView(idx) == true, "idx == " << idx );

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
  IndexType vidx = getFirstValidViewIndex();
  while ( indexIsValid(vidx) )
  {
    DataView * view = this->getView(vidx);
    delete view; 

    vidx = getNextValidViewIndex(vidx);
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
void DataGroup::destroyViewAndData( const std::string& name )
{
  SLIC_CHECK_MSG( hasView(name) == true, "name == " << name );

  DataView * view = detachView(name);
  if ( view != ATK_NULLPTR ) 
  { 
//
// RDH -- Should this check to see if there is only one view attached to 
//        the buffer before destroying it?
//
    DataBuffer * const buffer = view->getBuffer();
    delete view;

    if ( buffer != ATK_NULLPTR ) {
       getDataStore()->destroyBuffer(buffer->getIndex());
    } 
  }
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
void DataGroup::destroyViewAndData( IndexType idx )
{
  SLIC_CHECK_MSG( hasView(idx) == true, "idx == " << idx );

  DataView * view = detachView(idx);
  if ( view != ATK_NULLPTR )
  {
    // RDH TODO -- there should be a better way?
    DataBuffer * const buffer = view->getBuffer();
    delete view;

    if ( buffer != ATK_NULLPTR ) {
       getDataStore()->destroyBuffer(buffer->getIndex());
    } 
  }
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
void DataGroup::destroyViewsAndData()
{
  IndexType vidx = getFirstValidViewIndex();
  while ( indexIsValid(vidx) )
  {
    DataView * view = this->getView(vidx);

    // RDH TODO -- there should be a better way?
    DataBuffer * const buffer = view->getBuffer();
    delete view;

    if ( buffer != ATK_NULLPTR ) {
       getDataStore()->destroyBuffer(buffer->getIndex());
    } 

    vidx = getNextValidViewIndex(vidx);
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
DataView * DataGroup::moveView(DataView * view)
{
  SLIC_CHECK_MSG( view != ATK_NULLPTR,
                  "Attempting to move view, but given null ptr" );
  SLIC_CHECK_MSG( hasView(view->getName()) == false,
                  "Attempting to move view, but destination group already has a view named " << view->getName() );

  if ( view == ATK_NULLPTR || hasView(view->getName()) )
  {
    return ATK_NULLPTR;
  }
  else
  {  
    // remove this view from its current parent
    DataGroup * curr_group = view->getOwningGroup();

    curr_group->detachView(view->getName());

    /// finally, attach to this group
    attachView(view);

    return view;
  }
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
DataView * DataGroup::copyView(DataView * view)
{
  SLIC_CHECK_MSG( view != ATK_NULLPTR,
                  "Attempting to copy view, but given null ptr" );
  SLIC_CHECK_MSG( hasView(view->getName()) == false,
                  "Attempting to copy view, but destination group already has a view named " << view->getName() );

  if ( view == ATK_NULLPTR || hasView(view->getName()) )
  {
    return ATK_NULLPTR;
  }
  else
  {  
    DataView * res = createView(view->getName(), view->getBuffer());
    res->declare(view->getSchema());
    if (view->isApplied())
    {
      res->apply();
    }
    return res;
  }
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
DataGroup * DataGroup::createGroup( const std::string& name )
{
  SLIC_ASSERT( !name.empty() );
  SLIC_ASSERT_MSG( hasGroup(name) == false, "name == " << name );

  if ( name.empty() || hasGroup(name) )
  {
    return ATK_NULLPTR;
  }
  else
  {
    DataGroup * group = new DataGroup( name, this);
    return attachGroup(group);
  }
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
  SLIC_CHECK_MSG( hasGroup(name) == true, "name == " << name );

  delete detachGroup(name);
}

/*
 *************************************************************************
 *
 * Detach child group with given index and destroy it.
 *
 *************************************************************************
 */
void DataGroup::destroyGroup( IndexType idx )
{
  SLIC_CHECK_MSG( hasGroup(idx) == true, "idx == " << idx );

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
  IndexType gidx = getFirstValidGroupIndex();
  while ( indexIsValid(gidx) ) 
  {
    DataGroup * group = this->getGroup(gidx);
    delete group;

    gidx = getNextValidGroupIndex(gidx);
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
DataGroup * DataGroup::moveGroup(DataGroup * group)
{
  SLIC_CHECK_MSG( group != ATK_NULLPTR,
                  "Attempting to move group, but given null ptr" ); 
  SLIC_CHECK_MSG( hasGroup(group->getName()) == false,
                  "Attempting to move group, but destination group already has a group named " << group->getName() ); 

  if ( group == ATK_NULLPTR || hasGroup(group->getName()) )
  {
    return ATK_NULLPTR;
  }
  else 
  {
    // remove this group from its current parent
    DataGroup * curr_group = group->getParent();

    curr_group->detachGroup(group->getName());

    /// finally, attach to this group
    attachGroup(group);

    return group;
  }
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
DataGroup * DataGroup::copyGroup(DataGroup * group)
{
  SLIC_CHECK_MSG( group != ATK_NULLPTR,
                  "Attempting to move group, but given null ptr" );
  SLIC_CHECK_MSG( hasGroup(group->getName()) == false,
                  "Attempting to move group, but destination group already has a group named " << group->getName() );            

  if ( group == ATK_NULLPTR || hasGroup(group->getName()) )
  {
    return ATK_NULLPTR;
  }
  else 
  {
    DataGroup * res = createGroup(group->getName());

    // copy subgroups to new group
    IndexType gidx = group->getFirstValidGroupIndex();
    while ( indexIsValid(gidx) )
    {
      res->copyGroup(group->getGroup(gidx));
      gidx = group->getNextValidGroupIndex(gidx);
    }

    // copy views to new group
    IndexType vidx = group->getFirstValidViewIndex();
    while ( indexIsValid(vidx) )
    {
      res->copyView(group->getView(vidx));
      vidx = group->getNextValidViewIndex(vidx);
    }

    return res;
  } 
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

  IndexType vidx = getFirstValidViewIndex();
  while ( indexIsValid(vidx) )
  {
    const DataView * view = getView(vidx);
    Node& v = n["views"].fetch(view->getName());
    view->info(v);

    vidx = getNextValidViewIndex(vidx);
  }

  IndexType gidx = getFirstValidGroupIndex();
  while ( indexIsValid(gidx) )
  {
    const DataGroup * group =  getGroup(gidx);
    Node& g = n["groups"].fetch(group->getName());
    group->info(g);

    gidx = getNextValidGroupIndex(gidx);
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
  print(std::cout);
}

/*
 *************************************************************************
 *
 * Print JSON description of data group to an ostream.
 *
 *************************************************************************
 */
void DataGroup::print(std::ostream& os) const
{
  Node n;
  info(n);
  n.to_json_stream(os);
}

/*
 *************************************************************************
 *
 * Print given number of levels of group (sub) tree starting at this
 * group to stdout.
 *
 *************************************************************************
 */
void DataGroup::printTree( const int nlevels, 
                           std::ostream& os ) const
{
  for ( int i=0 ; i<nlevels ; ++i )
  {
    os <<"    ";
  }
  os <<"DataGroup "<<this->getName()<<std::endl;

  IndexType vidx = getFirstValidViewIndex();
  while ( indexIsValid(vidx) )
  {
    const DataView * view = getView(vidx);

    for ( int i=0 ; i<nlevels+1 ; ++i )
    {
      os <<"    ";
    }
    os << "DataView " << view->getName() << std::endl;

    vidx = getNextValidViewIndex(vidx);
  }

  IndexType gidx = getFirstValidGroupIndex();
  while ( indexIsValid(gidx) )
  {
    const DataGroup * group =  getGroup(gidx);

    group->printTree( nlevels + 1, os );

    gidx = getNextValidGroupIndex(gidx);
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
  if (protocol == "conduit")
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
  if (protocol == "conduit")
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
                     DataGroup * parent)
  : m_name(name),
  m_parent(parent),
  m_datastore(parent->getDataStore())
{}

/*
 *************************************************************************
 *
 * PRIVATE ctor makes group with given name and make it a child of
 * root group in datastore.
 *
 *************************************************************************
 */
DataGroup::DataGroup(const std::string& name,
                     DataStore * datastore)
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
 * PRIVATE method to create view and buffer and attach view to buffer.
 *
 *************************************************************************
 */
DataView * DataGroup::createViewAndBuffer( const std::string& name )
{
  DataView * view = new DataView( name, this );
  DataBuffer * buff = this->getDataStore()->createBuffer();
  view->attachBuffer( buff );
  return attachView(view);
}

/*
 *************************************************************************
 *
 * PRIVATE method to attach given view to group.
 *
 *************************************************************************
 */
DataView * DataGroup::attachView(DataView * view)
{
  if ( view == ATK_NULLPTR || hasView(view->getName()) )
  {
    return ATK_NULLPTR;
  } 
  else
  {
    m_view_coll.insertItem(view, view->getName());
    return view;
  }
}

/*
 *************************************************************************
 *
 * PRIVATE method to detach given with given name from group.
 *
 *************************************************************************
 */
DataView * DataGroup::detachView(const std::string& name )
{
  DataView * view = m_view_coll.removeItem(name);
  if ( view != ATK_NULLPTR )
  {
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
DataView * DataGroup::detachView(IndexType idx)
{
  DataView * view = m_view_coll.removeItem(idx);
  if ( view != ATK_NULLPTR )
  {
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
DataGroup * DataGroup::attachGroup(DataGroup * group)
{
  if ( group == ATK_NULLPTR || hasGroup(group->getName()) )
  {
    return ATK_NULLPTR;
  } 
  else
  {
    m_group_coll.insertItem(group, group->getName());
    return group;
  }
}

/*
 *************************************************************************
 *
 * PRIVATE method to detach child group with given name from group.
 *
 *************************************************************************
 */
DataGroup * DataGroup::detachGroup(const std::string& name )
{
  DataGroup * group = m_group_coll.removeItem(name);
  if ( group != ATK_NULLPTR )
  {
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
DataGroup * DataGroup::detachGroup(IndexType idx)
{
  DataGroup * group = m_group_coll.removeItem(idx);
  if ( group != ATK_NULLPTR )
  {
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
  std::vector<IndexType> buffer_ids;
  copyToNode(n,buffer_ids);

  // save the buffers discovered by buffer_ids
  for (size_t i=0 ; i < buffer_ids.size() ; i++)
  {
    Node& buff = n["buffers"].append();
    IndexType buffer_id = buffer_ids[i];
    DataBuffer * ds_buff =  m_datastore->getBuffer(buffer_id);
    buff["id"].set(buffer_id);
    DataType dtype = conduit::DataType::default_dtype(ds_buff->getTypeID());
    dtype.set_number_of_elements(ds_buff->getNumElements());
    buff["schema"].set(dtype.to_json());

    // only set our data if the buffer was initialized
    if (ds_buff->getVoidPtr() != ATK_NULLPTR )
    {
      buff["data"].set_external(dtype, ds_buff->getVoidPtr());
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
  std::map<IndexType, IndexType> id_map;
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
                           std::vector<IndexType>& buffer_ids) const
{
  IndexType vidx = getFirstValidViewIndex();
  while ( indexIsValid(vidx) )
  {
    const DataView * view = getView(vidx);
    Node& n_view = n["views"].fetch(view->getName());
    n_view["schema"].set(view->getSchema().to_json());
    n_view["is_applied"].set(view->isApplied());

    // if we have a buffer, simply add the index to the list
    if (view->hasBuffer())
    {
      IndexType buffer_id = view->getBuffer()->getIndex();
      n_view["buffer_id"].set(buffer_id);
      buffer_ids.push_back(view->getBuffer()->getIndex());
    }

    vidx = getNextValidViewIndex(vidx);
  }

  IndexType gidx = getFirstValidGroupIndex();
  while ( indexIsValid(gidx) )
  {
    const DataGroup * group =  getGroup(gidx);
    Node& n_group = n["groups"].fetch(group->getName());
    group->copyToNode(n_group, buffer_ids);

    gidx = getNextValidGroupIndex(gidx);
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
                             std::map<IndexType, IndexType>& id_map)
{
  /// for restore each group contains:
  /// buffers, views, and groups

  // create the buffers
  if (n.has_path("buffers"))
  {
    conduit::NodeIterator buffs_itr = n["buffers"].children();
    while (buffs_itr.has_next())
    {
      Node& n_buff = buffs_itr.next();
      IndexType buffer_id = n_buff["id"].as_int32();

      // create a new mapping and buffer if necessary
      if ( id_map.find(buffer_id) == id_map.end())
      {
        DataBuffer * ds_buff = this->getDataStore()->createBuffer();
        // map "id" to whatever new index the data store gives us.
        IndexType buffer_ds_id = ds_buff->getIndex();
        id_map[buffer_id] = buffer_ds_id;
        // setup the new data store buffer
        Schema schema(n_buff["schema"].as_string());
	TypeID type = static_cast<TypeID>(schema.dtype().id());
	SidreLength num_elems = schema.dtype().number_of_elements();
        ds_buff->declare(type, num_elems);
        if (n_buff.has_path("data"))
        {
          ds_buff->allocate();
          // copy the data from the node
	  void * data = n_buff["data"].element_ptr(0);
	  size_t nbytes = n_buff["data"].total_bytes();
	  ds_buff->update(data, nbytes);
        }
      }
    }
  }

  // create the child views
  conduit::NodeIterator views_itr = n["views"].children();
  while (views_itr.has_next())
  {
    Node& n_view = views_itr.next();
    if (n_view.has_path("buffer_id"))
    {
      std::string view_name = views_itr.path();

      IndexType buffer_id = n_view["buffer_id"].as_int32();
      // get the mapped buffer id

      SLIC_ASSERT_MSG( id_map.find(buffer_id) != id_map.end(),
                       "Invalid buffer index mapping." );

      buffer_id = id_map[buffer_id];
      DataBuffer * ds_buff = m_datastore->getBuffer(buffer_id);

      // create a new view with the buffer
      DataView * ds_view = createView(view_name, ds_buff);
      // declare using the schema
      Schema schema(n_view["schema"].as_string());
      ds_view->declare(schema);
      // if the schema was applied, restore this state
      if (n_view["is_applied"].to_uint64() != 0)
        ds_view->apply();
    }
    else
    {
      SLIC_ASSERT_MSG( 0, "DataGroup cannot restore opaque views." );
    }
  }

  // create the child groups
  conduit::NodeIterator groups_itr = n["groups"].children();
  while (groups_itr.has_next())
  {
    Node& n_group = groups_itr.next();
    std::string group_name = groups_itr.path();
    DataGroup * ds_group = createGroup(group_name);
    ds_group->copyFromNode(n_group, id_map);
  }
}


} /* end namespace sidre */
} /* end namespace asctoolkit */

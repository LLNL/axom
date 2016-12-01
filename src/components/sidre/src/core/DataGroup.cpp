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

#include "conduit_relay.hpp"
#include "conduit_relay_hdf5.hpp"

// Other toolkit component headers

// SiDRe project headers
#include "DataBuffer.hpp"
#include "DataStore.hpp"
#include "SidreUtilities.hpp"

namespace asctoolkit
{
namespace sidre
{

// Initialization of static path delimiter character for methods that
// support path syntax.
const char DataGroup::s_path_delimiter = '/';

////////////////////////////////////////////////////////////////////////
//
// Basic query and accessor methods.
//
////////////////////////////////////////////////////////////////////////

/*
 *************************************************************************
 *
 * Return path of Group object, not including its name.
 *
 *************************************************************************
 */
std::string DataGroup::getPath() const
{
  const DataGroup * root = getDataStore()->getRoot();
  const DataGroup * curr = getParent();
  std::string thePath = curr->getName();
  curr = curr->getParent();

  while (curr != root)
  {
    thePath = curr->getName() + s_path_delimiter + thePath;
    curr = curr->getParent();
  }

  return thePath;
}

////////////////////////////////////////////////////////////////////////
//
// View query methods.
//
////////////////////////////////////////////////////////////////////////

/*
 *************************************************************************
 *
 * Return true if Group owns a View with given name or path; else false.
 *
 *************************************************************************
 */
bool DataGroup::hasView( const std::string& path ) const
{
  std::string intpath(path);
  const DataGroup * group = walkPath( intpath );

  if (group == ATK_NULLPTR)
  {
    return false;
  }

  return group->hasChildView(intpath);
}


////////////////////////////////////////////////////////////////////////
//
// View access methods.
//
////////////////////////////////////////////////////////////////////////

/*
 *************************************************************************
 *
 * Return pointer to non-const View with given name or path if it exists.
 *
 *************************************************************************
 */
DataView * DataGroup::getView( const std::string& path )
{
  std::string intpath(path);
  bool create_groups_in_path = false;
  DataGroup * group = walkPath( intpath, create_groups_in_path );

  if ( group == ATK_NULLPTR )
  {
    SLIC_CHECK_MSG( group != ATK_NULLPTR,
		    "Non-existent group in path " << path );
    return ATK_NULLPTR;
  }

  SLIC_CHECK_MSG( !intpath.empty() && group->hasChildView(intpath),
                  "Group " << getName() <<
                  " has no View with name '" << intpath << "'");

  return group->m_view_coll.getItem(intpath);
}

/*
 *************************************************************************
 *
 * Return pointer to const View with given name or path if it exists.
 *
 *************************************************************************
 */
const DataView * DataGroup::getView( const std::string& path ) const
{
  std::string intpath(path);
  const DataGroup * group = walkPath( intpath );

  if (group == ATK_NULLPTR)
  {
    SLIC_CHECK_MSG( group != ATK_NULLPTR,
		    "Non-existent group in path " << path );
    return ATK_NULLPTR;
  }

  SLIC_CHECK_MSG( !intpath.empty() && group->hasChildView(intpath),
		  "Group " << getName() <<
		  " has no View with name '" << intpath << "'");
  
  return group->m_view_coll.getItem(intpath);
}


////////////////////////////////////////////////////////////////////////
//
//  Methods to create a View that has no associated data.
//
////////////////////////////////////////////////////////////////////////

/*
 *************************************************************************
 *
 * Create empty View (i.e., no data description) with given name or path
 * in this Group.
 *
 *************************************************************************
 */
DataView * DataGroup::createView( const std::string& path )
{
  std::string intpath(path);
  bool create_groups_in_path = true;
  DataGroup * group = walkPath( intpath, create_groups_in_path );

  if ( group == ATK_NULLPTR )
  {
    SLIC_CHECK_MSG( group != ATK_NULLPTR,
		    "Could not find or create path " << path <<
		    " since it appears there is already a view with that name" );
    return ATK_NULLPTR;
  }
  else if ( intpath.empty() || group->hasChildView(intpath) || group->hasChildGroup(intpath) )
  {
    SLIC_CHECK( !intpath.empty() );
    SLIC_CHECK_MSG( !group->hasChildView(intpath),
                    "Cannot create View with name '" << intpath <<
                    "' in Group '" << getName() <<
                    " since it already has a View with that name" );
    SLIC_CHECK_MSG( !group->hasChildGroup(intpath),
                    "Cannot create View with name '" << intpath <<
                    "' in Group '" << getName() <<
                    " since it already has a Group with that name" );
    return ATK_NULLPTR;
  }

  DataView * view = new(std::nothrow) DataView(intpath);
  if ( view != ATK_NULLPTR )
  {
    group->attachView(view);
  }
  return view;
}

/*
 *************************************************************************
 *
 * Create described View (type and # elems) with given name or path in
 * this Group.
 *
 *************************************************************************
 */
DataView * DataGroup::createView( const std::string& path,
                                  TypeID type,
                                  SidreLength num_elems )
{
  if ( type == NO_TYPE_ID || num_elems < 0 )
  {
    SLIC_CHECK_MSG(type != NO_TYPE_ID,
                   "Cannot create View with name '" << path <<
                   "' in Group '" << getName() << " without a valid type" );
    SLIC_CHECK_MSG(num_elems >= 0,
                   "Cannot create View with name '" << path <<
                   "' in Group '" << getName() << " with # elems < 0" );
    return ATK_NULLPTR;
  }

  DataView * view = createView(path);
  if (view != ATK_NULLPTR)
  {
    view->describe(type, num_elems);
  }
  return view;
}

/*
 *************************************************************************
 *
 * Create described View (type and shape) with given name or path in this
 * Group.
 *
 *************************************************************************
 */
DataView * DataGroup::createView( const std::string& path,
                                  TypeID type,
                                  int ndims,
                                  SidreLength * shape )
{
  if ( type == NO_TYPE_ID || ndims < 0 || shape == ATK_NULLPTR )
  {
    SLIC_CHECK_MSG(type != NO_TYPE_ID,
                   "Cannot create View with name '" << path <<
                   "' in Group '" << getName() << " without a valid type" );
    SLIC_CHECK_MSG(ndims >= 0,
                   "Cannot create View with name '" << path <<
                   "' in Group '" << getName() << " with ndims < 0" );
    SLIC_CHECK_MSG(shape != ATK_NULLPTR,
                   "Cannot create View with name '" << path <<
                   "' in Group '" << getName() << " with null shape ptr" );
    return ATK_NULLPTR;
  }

  DataView * view = createView(path);
  if (view != ATK_NULLPTR)
  {
    view->describe(type, ndims, shape);
  }
  return view;
}

/*
 *************************************************************************
 *
 * Create View with given name and data described by a conduit Datatype object.
 *
 *************************************************************************
 */
DataView * DataGroup::createView( const std::string& path,
                                  const DataType& dtype )
{
  DataView * view = createView(path);
  if (view != ATK_NULLPTR)
  {
    view->describe(dtype);
  }

  return view;
}


////////////////////////////////////////////////////////////////////////
//
//  Methods to create a View and attach a Buffer to it.
//
////////////////////////////////////////////////////////////////////////

/*
 *************************************************************************
 *
 * Create an undescribed View with given name or path in
 * this Group and attach Buffer to it.
 *
 *************************************************************************
 */
DataView * DataGroup::createView( const std::string& path,
                                  DataBuffer * buff )
{
  DataView * view = createView(path);
  if ( view != ATK_NULLPTR )
  {
    view->attachBuffer( buff );
  }
  return view;
}

/*
 *************************************************************************
 *
 * Create described View (type and # elems) with given name or path in
 * this Group and attach Buffer to it.
 *
 *************************************************************************
 */
DataView * DataGroup::createView( const std::string& path,
                                  TypeID type,
                                  SidreLength num_elems,
                                  DataBuffer * buff )
{
  DataView * view = createView(path, type, num_elems);
  if (view != ATK_NULLPTR)
  {
    view->attachBuffer(buff);
  }
  return view;
}

/*
 *************************************************************************
 *
 * Create described View (type and shape) with given name or path in
 * this Group and attach Buffer to it.
 *
 *************************************************************************
 */
DataView * DataGroup::createView( const std::string& path,
                                  TypeID type,
                                  int ndims,
                                  SidreLength * shape,
                                  DataBuffer * buff )
{
  DataView * view = createView(path, type, ndims, shape);
  if (view != ATK_NULLPTR)
  {
    view->attachBuffer(buff);
  }
  return view;
}

/*
 *************************************************************************
 *
 * Create described View (Conduit DataType) with given name or path in
 * this Group and attach Buffer to it.
 *
 *************************************************************************
 */
DataView * DataGroup::createView( const std::string& path,
                                  const DataType& dtype,
                                  DataBuffer * buff )
{
  DataView * view = createView(path, dtype);
  if (view != ATK_NULLPTR)
  {
    view->attachBuffer(buff);
  }

  return view;
}


////////////////////////////////////////////////////////////////////////
//
//  Methods to create a View and attach external data to it.
//
////////////////////////////////////////////////////////////////////////

/*
 *************************************************************************
 *
 * Create an undescribed View with given name or path in
 * this Group and attach external data ptr to it.
 *
 *************************************************************************
 */
DataView * DataGroup::createView( const std::string& path,
                                  void * external_ptr )
{
  DataView * view = createView(path);
  if ( view != ATK_NULLPTR )
  {
    view->setExternalDataPtr(external_ptr);
  }
  return view;
}

/*
 *************************************************************************
 *
 * Create described View (type and # elems) with given name or path in
 * this Group and attach external data ptr to it.
 *
 *************************************************************************
 */
DataView * DataGroup::createView( const std::string& path,
                                  TypeID type,
                                  SidreLength num_elems,
                                  void * external_ptr )
{
  DataView * view = createView(path, type, num_elems);
  if (view != ATK_NULLPTR)
  {
    view->setExternalDataPtr(external_ptr);
  }
  return view;
}

/*
 *************************************************************************
 *
 * Create described View (type and shape) with given name or path in
 * this Group and attach external data ptr to it.
 *
 *************************************************************************
 */
DataView * DataGroup::createView( const std::string& path,
                                  TypeID type,
                                  int ndims,
                                  SidreLength * shape,
                                  void * external_ptr )
{
  DataView * view = createView(path, type, ndims, shape);
  if (view != ATK_NULLPTR)
  {
    view->setExternalDataPtr(external_ptr);
  }
  return view;
}

/*
 *************************************************************************
 *
 * Create described View (Conduit DataType) with given name or path in
 * this Group and attach external data ptr to it.
 *
 *************************************************************************
 */
DataView * DataGroup::createView( const std::string& path,
                                  const DataType& dtype,
                                  void * external_ptr )
{
  DataView * view = createView(path, dtype);
  if (view != ATK_NULLPTR)
  {
    view->setExternalDataPtr(external_ptr);
  }
  return view;
}


////////////////////////////////////////////////////////////////////////
//
//  Methods to create a View and allocate its data.
//
////////////////////////////////////////////////////////////////////////

/*
 *************************************************************************
 *
 * Create described View (type and # elems) with given name or path in
 * this Group and allocate its data.
 *
 *************************************************************************
 */
DataView * DataGroup::createViewAndAllocate( const std::string& path,
                                             TypeID type,
                                             SidreLength num_elems )
{
  DataView * view = createView(path, type, num_elems);
  if ( view != ATK_NULLPTR )
  {
    view->allocate();
  }
  return view;
}

/*
 *************************************************************************
 *
 * Create described View (type and shape) with given name or path in
 * this Group and allocate its data.
 *
 *************************************************************************
 */
DataView * DataGroup::createViewAndAllocate( const std::string& path,
                                             TypeID type,
                                             int ndims,
                                             SidreLength * shape )
{
  DataView * view = createView(path, type, ndims, shape);
  if ( view != ATK_NULLPTR )
  {
    view->allocate();
  }
  return view;
}

/*
 *************************************************************************
 *
 * Create described View (Conduit DataType) with given name or path in
 * this Group and allocate its data.
 *
 *************************************************************************
 */
DataView * DataGroup::createViewAndAllocate( const std::string& path,
                                             const DataType& dtype)
{
  DataView * view = createView(path, dtype);
  if ( view != ATK_NULLPTR )
  {
    view->allocate();
  }
  return view;
}

/*
 *************************************************************************
 *
 * Create View with given name or path and set its data to given string.
 *
 *************************************************************************
 */
DataView * DataGroup::createViewString( const std::string& path,
                                        const std::string& value)
{
  DataView * view = createView(path);
  if (view != ATK_NULLPTR)
  {
    view->setString(value);
  }

  return view;
}


////////////////////////////////////////////////////////////////////////
//
//  Methods for destroying Views and their data.
//
////////////////////////////////////////////////////////////////////////

/*
 *************************************************************************
 *
 * Destroy View with given name or path and leave its data intact.
 *
 *************************************************************************
 */
void DataGroup::destroyView( const std::string& path )
{
  std::string intpath(path);
  bool create_groups_in_path = false;
  DataGroup * group = walkPath( intpath, create_groups_in_path );

  if ( group != ATK_NULLPTR )
  {
    DataView * view = group->detachView(intpath);
    if ( view != ATK_NULLPTR )
    {
      delete view;
    }
  }
}

/*
 *************************************************************************
 *
 * Destroy View with given index and leave its data intact.
 *
 *************************************************************************
 */
void DataGroup::destroyView( IndexType idx )
{
  DataView * view = detachView(idx);
  if ( view != ATK_NULLPTR )
  {
    delete view;
  }
}

/*
 *************************************************************************
 *
 * Destroy all Views in Group and leave their data intact.
 *
 *************************************************************************
 */
void DataGroup::destroyViews()
{
  IndexType vidx = getFirstValidViewIndex();
  while ( indexIsValid(vidx) )
  {
    DataView * view = detachView(vidx);
    if ( view != ATK_NULLPTR )
    {
      delete view;
    }

    vidx = getNextValidViewIndex(vidx);
  }

  m_view_coll.removeAllItems();
}

/*
 *************************************************************************
 *
 * Destroy View with given name or path and its data if it's the only View
 * associated with that data.
 *
 *************************************************************************
 */
void DataGroup::destroyViewAndData( const std::string& path )
{
  destroyViewAndData(getView(path));
}

/*
 *************************************************************************
 *
 * Destroy View with given index and its data if it's the only View
 * associated with that data.
 *
 *************************************************************************
 */
void DataGroup::destroyViewAndData( IndexType idx )
{
  destroyViewAndData(getView(idx));
}

/*
 *************************************************************************
 *
 * Destroy all Views in Group as well as the data for each View when it's
 * the only View associated with that data.
 *
 *************************************************************************
 */
void DataGroup::destroyViewsAndData()
{
  IndexType vidx = getFirstValidViewIndex();
  while ( indexIsValid(vidx) )
  {
    destroyViewAndData(vidx);
    vidx = getNextValidViewIndex(vidx);
  }

  m_view_coll.removeAllItems();
}


////////////////////////////////////////////////////////////////////////
//
//  Methods for moving and copying View objects from one Group to another.
//
////////////////////////////////////////////////////////////////////////

/*
 *************************************************************************
 *
 * Remove given View from its owning Group and attach to this Group.
 *
 *************************************************************************
 */
DataView * DataGroup::moveView(DataView * view)
{
  if ( view == ATK_NULLPTR )
  {
    SLIC_CHECK( view != ATK_NULLPTR );
    return ATK_NULLPTR;
  }

  DataGroup * curr_group = view->getOwningGroup();
  if (curr_group == this)
  {
    // this Group already owns the View
    return view;
  }
  else if (hasChildView(view->getName()))
  {
    SLIC_CHECK_MSG(!hasChildView(view->getName()),
                   "Group '" << getName() <<
                   "' already has a View named'" << view->getName() <<
                   "' so View move operation cannot happen");
    return ATK_NULLPTR;
  }

  curr_group->detachView(view);
  attachView(view);

  return view;
}

/*
 *************************************************************************
 *
 * Create a copy of given View and attach to this Group.
 *
 * Copying a View does not perform a deep copy of its data Buffer.
 *
 *************************************************************************
 */
DataView * DataGroup::copyView(DataView * view)
{
  if ( view == ATK_NULLPTR || hasChildView(view->getName()) )
  {
    SLIC_CHECK( view != ATK_NULLPTR );
    SLIC_CHECK_MSG(!hasChildView(view->getName()),
                   "Group '" << getName() <<
                   "' already has a View named'" << view->getName() <<
                   "' so View copy operation cannot happen");

    return ATK_NULLPTR;
  }

  DataView * copy = createView(view->getName());
  view->copyView(copy);
  return copy;
}

////////////////////////////////////////////////////////////////////////
//
// Child Group query methods.
//
////////////////////////////////////////////////////////////////////////


/*
 ***********************************************************************
 *
 * Return true if this Group has a descendant Group with given name or path;
 * else false.
 *
 ***********************************************************************
 */
bool DataGroup::hasGroup( const std::string& path ) const
{
  std::string intpath(path);
  const DataGroup * group = walkPath( intpath );

  if ( group == ATK_NULLPTR )
  {
    return false;
  }
  else
  {
    return group->hasChildGroup(intpath);
  }
}

////////////////////////////////////////////////////////////////////////
//
// Child Group access methods.
//
////////////////////////////////////////////////////////////////////////

/*
 *************************************************************************
 *
 * Return pointer to non-const child Group with given name or path
 * if it exists.
 *
 *************************************************************************
 */
DataGroup * DataGroup::getGroup( const std::string& path )
{
  std::string intpath(path);
  bool create_groups_in_path = false;
  DataGroup * group = walkPath( intpath, create_groups_in_path );

  if (group == ATK_NULLPTR)
  {
    SLIC_CHECK_MSG( group != ATK_NULLPTR,
		    "Non-existent group in path " << path );
    return ATK_NULLPTR;
  }

  SLIC_CHECK_MSG( !intpath.empty() && group->hasChildGroup(intpath),
                  "Group " << getName() <<
                  " has no descendant Group with name '" << path << "'");

  return group->m_group_coll.getItem(intpath);
}

/*
 *************************************************************************
 *
 * Return pointer to const child Group with given name or path if it exists.
 *
 *************************************************************************
 */
const DataGroup * DataGroup::getGroup( const std::string& path ) const
{
  std::string intpath(path);
  const DataGroup * group = walkPath( intpath );

  if (group == ATK_NULLPTR)
  {
    SLIC_CHECK_MSG( group != ATK_NULLPTR,
		    "Non-existent group in path " << path );
    return ATK_NULLPTR;
  }

  SLIC_CHECK_MSG( !intpath.empty() && group->hasChildGroup(intpath),
                  "Group " << getName() <<
                  " has no descendant Group with name '" << path << "'");

  return group->m_group_coll.getItem(intpath);
}


////////////////////////////////////////////////////////////////////////
//
//  Methods for managing child DataGroup objects in DataGroup
//
////////////////////////////////////////////////////////////////////////


/*
 *************************************************************************
 *
 * Create Group with given name or path and make it a child of this Group.
 *
 *************************************************************************
 */
DataGroup * DataGroup::createGroup( const std::string& path )
{
  std::string intpath(path);
  bool create_groups_in_path = true;
  DataGroup * group = walkPath( intpath, create_groups_in_path );

  if ( group == ATK_NULLPTR )
  {
    SLIC_CHECK_MSG( group != ATK_NULLPTR,
		    "Could not find or create path " << path <<
		    " since it appears there is already a view with that name" );
    return ATK_NULLPTR;
  }
  else if ( intpath.empty() || group->hasChildGroup(intpath) || group->hasChildView(intpath) )
  {
    SLIC_CHECK( !intpath.empty() );
    SLIC_CHECK_MSG( !group->hasChildGroup(intpath),
                    "Cannot create Group with name '" << path <<
                    " in Group '" << getName() <<
                    " since it already has a Group with that name" );
    SLIC_CHECK_MSG( !group->hasChildView(intpath),
                    "Cannot create Group with name '" << path <<
                    " in Group '" << getName() <<
                    " since it already has a View with that name" );

    return ATK_NULLPTR;
  }

  DataGroup * new_group = new(std::nothrow) DataGroup(intpath, group);
  if ( new_group == ATK_NULLPTR )
  {
    return ATK_NULLPTR;
  }
  return group->attachGroup(new_group);
}

/*
 *************************************************************************
 *
 * Detach child Group with given name or path and destroy it.
 *
 *************************************************************************
 */
void DataGroup::destroyGroup( const std::string& path )
{
  std::string intpath(path);
  bool create_groups_in_path = false;
  DataGroup * group = walkPath( intpath, create_groups_in_path );

  if ( group != ATK_NULLPTR )
  {
    DataGroup * targetgroup = group->detachGroup(intpath);
    if ( targetgroup != ATK_NULLPTR )
    {
      delete targetgroup;
    }
  }
}

/*
 *************************************************************************
 *
 * Detach child Group with given index and destroy it.
 *
 *************************************************************************
 */
void DataGroup::destroyGroup( IndexType idx )
{
  DataGroup * group = detachGroup(idx);
  if ( group != ATK_NULLPTR )
  {
    delete group;
  }
}

/*
 *************************************************************************
 *
 * Detach all child Groups and destroy them.
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
 * Remove given Group from its owning Group and make it a child of this Group.
 *
 *************************************************************************
 */
DataGroup * DataGroup::moveGroup(DataGroup * group)
{
  if ( group == ATK_NULLPTR || hasChildGroup(group->getName()))
  {
    SLIC_CHECK( group != ATK_NULLPTR );
    SLIC_CHECK_MSG(!hasChildGroup(group->getName()),
                   "Group '" << getName() <<
                   "' already has a child Group named'" << group->getName() <<
                   "' so Group move operation cannot happen");

    return ATK_NULLPTR;
  }

  DataGroup * curr_group = group->getParent();
  curr_group->detachGroup(group->getName());
  attachGroup(group);
  return group;
}

/*
 *************************************************************************
 *
 * Create a copy of given Group and make it a child of this Group.
 *
 * Copying a Group does not perform a deep copy of any of its Buffers.
 *
 *************************************************************************
 */
DataGroup * DataGroup::copyGroup(DataGroup * group)
{
  if ( group == ATK_NULLPTR || hasChildGroup(group->getName()) )
  {
    SLIC_CHECK( group != ATK_NULLPTR );
    SLIC_CHECK_MSG(!hasChildGroup(group->getName()),
                   "Group '" << getName() <<
                   "' already has a child Group named'" << group->getName() <<
                   "' so Group copy operation cannot happen");

    return ATK_NULLPTR;
  }
  else
  {
    DataGroup * res = createGroup(group->getName());

    // copy child Groups to new Group
    IndexType gidx = group->getFirstValidGroupIndex();
    while ( indexIsValid(gidx) )
    {
      res->copyGroup(group->getGroup(gidx));
      gidx = group->getNextValidGroupIndex(gidx);
    }

    // copy Views to new Group
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
 * Copy Group native layout to given Conduit node.
 *
 *************************************************************************
 */
void DataGroup::createNativeLayout(Node& n) const
{
  n.set(DataType::object());

  // Dump the group's views
  IndexType vidx = getFirstValidViewIndex();
  while ( indexIsValid(vidx) )
  {
    const DataView * view = getView(vidx);

    // Check that the view's name is not also a child group name
    SLIC_CHECK_MSG( !hasChildGroup(view->getName())
                    , view->getName() << " is the name of a groups and a view");

    view->createNativeLayout( n[view->getName()] );
    vidx = getNextValidViewIndex(vidx);
  }

  // Recursively dump the child groups
  IndexType gidx = getFirstValidGroupIndex();
  while ( indexIsValid(gidx) )
  {
    const DataGroup * group =  getGroup(gidx);
    group->createNativeLayout(n[group->getName()]);
    gidx = getNextValidGroupIndex(gidx);
  }
}

/*
 *************************************************************************
 *
 * Copy Group native layout to given Conduit node.
 *
 *************************************************************************
 * see ATK-786 - Improvements to createNativeLayout and createExternalLayout
 */
void DataGroup::createExternalLayout(Node& n) const
{
  n.set(DataType::object());

  // Dump the group's views
  IndexType vidx = getFirstValidViewIndex();
  while ( indexIsValid(vidx) )
  {
    const DataView * view = getView(vidx);

    // Check that the view's name is not also a child group name
    SLIC_CHECK_MSG( !hasChildGroup(view->getName())
                    , view->getName() << " is the name of a groups and a view");

    if(view->isExternal() && view->isDescribed())
    {
        view->createNativeLayout(  n[view->getName()]  );
    }

    vidx = getNextValidViewIndex(vidx);
  }

  // Recursively dump the child groups
  IndexType gidx = getFirstValidGroupIndex();
  while ( indexIsValid(gidx) )
  {
    const DataGroup * group =  getGroup(gidx);
    group->createExternalLayout(n[group->getName()]);
    gidx = getNextValidGroupIndex(gidx);
  }
}

/*
 *************************************************************************
 *
 * Print JSON description of data Group to stdout.
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
 * Print JSON description of data Group to an ostream.
 *
 *************************************************************************
 */
void DataGroup::print(std::ostream& os) const
{
  Node n;
  copyToConduitNode(n);
  n.to_json_stream(os);
}

/*
 *************************************************************************
 *
 * Print given number of levels of Group sub-tree rooted at this
 * Group to stdout.
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
  os << "Group "<< this->getName() <<std::endl;

  IndexType vidx = getFirstValidViewIndex();
  while ( indexIsValid(vidx) )
  {
    const DataView * view = getView(vidx);

    for ( int i=0 ; i<nlevels+1 ; ++i )
    {
      os <<"    ";
    }
    os << "View " << view->getName() << std::endl;

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
 * Copy data Group description to given Conduit node.
 *
 *************************************************************************
 */
void DataGroup::copyToConduitNode(Node& n) const
{
  n["name"] = m_name;

  IndexType vidx = getFirstValidViewIndex();
  while ( indexIsValid(vidx) )
  {
    const DataView * view = getView(vidx);
    Node& v = n["views"].fetch(view->getName());
    view->copyToConduitNode(v);

    vidx = getNextValidViewIndex(vidx);
  }

  IndexType gidx = getFirstValidGroupIndex();
  while ( indexIsValid(gidx) )
  {
    const DataGroup * group =  getGroup(gidx);
    Node& g = n["groups"].fetch(group->getName());
    group->copyToConduitNode(g);

    gidx = getNextValidGroupIndex(gidx);
  }
}

/*
 *************************************************************************
 *
 * Test this Group for equivalence to another Group.
 *
 *************************************************************************
 */
bool DataGroup::isEquivalentTo(const DataGroup * other) const
{
  // Equality of names
  bool is_equiv = (m_name == other->m_name);

  // Sizes of collections of child items must be equal
  if (is_equiv)
  {
    is_equiv = (m_view_coll.getNumItems() == other->m_view_coll.getNumItems())
               && (m_group_coll.getNumItems() ==
                   other->m_group_coll.getNumItems());
  }

  // Test equivalence of Views
  if (is_equiv)
  {
    IndexType vidx = getFirstValidViewIndex();
    while ( is_equiv && indexIsValid(vidx) )
    {
      const DataView * view = getView(vidx);
      const std::string& name = view->getName();

      is_equiv = other->hasChildView( name )
              && view->isEquivalentTo( other->getView( name ) );

      vidx = getNextValidViewIndex(vidx);
    }
  }

  // Recursively call this method to test equivalence of child DataGroups
  if (is_equiv)
  {
    IndexType gidx = getFirstValidGroupIndex();
    while ( is_equiv && indexIsValid(gidx) )
    {
      const DataGroup * group =  getGroup(gidx);
      const std::string& name = group->getName();

      is_equiv = other->hasChildGroup( name )
              && group->isEquivalentTo( other->getGroup( name ));

      gidx = getNextValidGroupIndex(gidx);
    }
  }

  return is_equiv;
}

/*
 *************************************************************************
 *
 * Save Group (including Views and child Groups) to a file
 *
 *************************************************************************
 */

void DataGroup::save(const std::string& path,
                     const std::string& protocol) const
{

  if (protocol == "sidre_hdf5")
  {
    Node n;
    exportTo(n["sidre"]);
    createExternalLayout(n["sidre/external"]);
    conduit::relay::io::save(n, path, "hdf5");
  }
  else if (protocol == "sidre_conduit_json")
  {
    Node n;
    exportTo(n["sidre"]);
    createExternalLayout(n["sidre/external"]);
    conduit::relay::io::save(n, path, "conduit_json");
  }
  else if (protocol == "sidre_json")
  {
    Node n;
    exportTo(n["sidre"]);
    createExternalLayout(n["sidre/external"]);
    conduit::relay::io::save(n, path, "json");
  }
  else if (protocol == "conduit_hdf5" )
  {
      Node n;
      createNativeLayout(n);
      conduit::relay::io::save(n, path,"hdf5");
  }
  else if (protocol == "conduit_bin"  || 
           protocol == "conduit_json" ||
           protocol == "json")
  {
    Node n;
    createNativeLayout(n);
    conduit::relay::io::save(n, path, protocol);
  }
  else
  {
    SLIC_ERROR("Invalid protocol " << protocol << " for file load.");
  }
}

/*
 *************************************************************************
 *
 * Save Group (including Views and child Groups) to a hdf5 handle
 *
 *************************************************************************
 */
void DataGroup::save(const hid_t& h5_id,
                     const std::string& protocol) const
{
  // supported here:
  // "sidre_hdf5"
  // "conduit_hdf5"
  if(protocol == "sidre_hdf5")
  {
    Node n;
    exportTo(n["sidre"]);
    createExternalLayout(n["sidre/external"]);
    conduit::relay::io::hdf5_write(n,h5_id);
  }
  else if( protocol == "conduit_hdf5")
  {
    Node n;
    createNativeLayout(n);
    conduit::relay::io::hdf5_write(n, h5_id);
  }
  else
  {
    SLIC_ERROR("Invalid protocol " 
                << protocol 
                << " for save with hdf5 handle.");
  }
}


/*************************************************************************/

/*
 *************************************************************************
 *
 * Load Group (including Views and child Groups) from a file
 *
 *************************************************************************
 */
void DataGroup::load(const std::string& path,
                     const std::string& protocol)
{

  if (protocol == "sidre_hdf5")
  {
    Node n;
    conduit::relay::io::load(path,"hdf5", n);
    SLIC_ASSERT(n.has_path("sidre"));
    importFrom(n["sidre"]);
  }
  else if (protocol == "sidre_conduit_json")
  {
    Node n;
    conduit::relay::io::load(path,"conduit_json", n);
    SLIC_ASSERT(n.has_path("sidre"));
    importFrom(n["sidre"]);
  }
  else if (protocol == "sidre_json")
  {
    Node n;
    conduit::relay::io::load(path,"json", n);
    SLIC_ASSERT(n.has_path("sidre"));
    importFrom(n["sidre"]);
  }
  else if (protocol == "conduit_hdf5")
  {
    Node n;
    conduit::relay::io::load(path,"hdf5", n);
    importConduitTree(n);
  }
  else if (protocol == "conduit_bin"  || 
           protocol == "conduit_json" ||
           protocol == "json")
  {
    Node n;
    conduit::relay::io::load(path,protocol, n);
    importConduitTree(n);
  }
  else
  {
    SLIC_ERROR("Invalid protocol " << protocol << " for file load.");
  }
}

/*
 *************************************************************************
 *
 * Load Group (including Views and child Groups) from an hdf5 handle
 *
 *************************************************************************
 */
void DataGroup::load(const hid_t& h5_id,
                     const std::string &protocol)
{
  // supported here:
  // "sidre_hdf5"
  // "conduit_hdf5"
  if(protocol == "sidre_hdf5")
  {
    Node n;
    conduit::relay::io::hdf5_read(h5_id,n);
    SLIC_ASSERT(n.has_path("sidre"));
    importFrom(n["sidre"]);
  }
  else if( protocol == "conduit_hdf5")
  {
    SLIC_ERROR("Protocol " << protocol << " not yet supported for file load.");
    Node n;
    conduit::relay::io::hdf5_read(h5_id, n);
    importConduitTree(n);
  }
  else
  {
    SLIC_ERROR("Invalid protocol " << protocol << " for file load.");
  }
}

/*
 *************************************************************************
 *
 * Load External Data from a file
 * 
 *************************************************************************
 */
void DataGroup::loadExternalData(const std::string& path)
{
  Node n;
  createExternalLayout(n);
  // CYRUS'-NOTE, not sure ":" will work with multiple trees per
  // output file
  conduit::relay::io::hdf5_read( path + ":sidre/external", n);
}

/*
 *************************************************************************
 *
 * Load External Data from an hdf5 file
 *
 * Note: this ASSUMES uses the "sidre_hdf5" protocol
 *************************************************************************
 */
void DataGroup::loadExternalData(const hid_t& h5_id)
{
  Node n;
  createExternalLayout(n);
  conduit::relay::io::hdf5_read(h5_id, "sidre/external", n);
}


////////////////////////////////////////////////////////////////////////
//
// Private methods below
//
////////////////////////////////////////////////////////////////////////

/*
 *************************************************************************
 *
 * PRIVATE ctor makes Group with given name and make it a child of parent.
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
 * PRIVATE ctor makes Group with given name and make it a child of
 * root Group in datastore.
 *
 *************************************************************************
 */
DataGroup::DataGroup(const std::string& name,
                     DataStore * datastore)
  : m_name(name),
  m_parent(datastore->getRoot()),
  m_datastore(datastore)
{}

/*
 *************************************************************************
 *
 * PRIVATE dtor destroys Group and all its contents.
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
 * PRIVATE method to attach given View to Group.
 *
 *************************************************************************
 */
DataView * DataGroup::attachView(DataView * view)
{
  if ( view == ATK_NULLPTR || hasChildView(view->getName()) )
  {
    return ATK_NULLPTR;
  }
  else
  {
    SLIC_ASSERT(view->m_owning_group == ATK_NULLPTR);
    view->m_owning_group = this;
    m_view_coll.insertItem(view, view->getName());
    return view;
  }
}

/*
 *************************************************************************
 *
 * PRIVATE method to detach View with given name from Group.
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
 * PRIVATE method to detach View with given index from Group.
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
 * PRIVATE method to destroy View in this Group and its data.
 *
 *************************************************************************
 */
void DataGroup::destroyViewAndData( DataView * view )
{
  if ( view != ATK_NULLPTR )
  {
    DataGroup * group = view->getOwningGroup();
    group->detachView( view->getName() );
    DataBuffer * const buffer = view->detachBuffer();
    if ( buffer != ATK_NULLPTR && buffer->getNumViews() == 0 )
    {
      getDataStore()->destroyBuffer(buffer);
    }
    delete view;
  }
}

/*
 *************************************************************************
 *
 * PRIVATE method to make given Group a child of this Group.
 *
 *************************************************************************
 */
DataGroup * DataGroup::attachGroup(DataGroup * group)
{
  if ( group == ATK_NULLPTR || hasChildGroup(group->getName()) )
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
 * PRIVATE method to detach child Group with given name from Group.
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
 * PRIVATE method to detach child Group with given index from Group.
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
 * Serialize tree identified by a Group into a conduit node.  Include
 * any Buffers attached to Views in that tree.
 *
 *
 * Note: this is for the "sidre_hdf5" protocol
 *
 *************************************************************************
 */
void DataGroup::exportTo(conduit::Node & result) const
{ 
  result.set(DataType::object());
  // TODO - This implementation will change in the future.  We want to write
  // out some separate set of conduit nodes:
  // #1 A set of nodes representing the Group and Views (hierarchy), with
  // the data descriptions ( schemas ).
  // #2 A set of nodes for our data ( Buffers, external data, etc ).
  // On a load, we want to be able to create our DataStore tree first,
  // then call allocate ourself, then have conduit load the data directly
  // into our allocated memory areas.  Conduit can do this, as long as the
  // conduit node set is compatible with what's in the file.
  std::set<IndexType> buffer_indices;

  // Tell Group to add itself and all sub-Groups and Views to node.
  // Any Buffers referenced by those Views will be tracked in the
  // buffer_indices
  exportTo(result, buffer_indices);

  if (!buffer_indices.empty())
  {
    // Now, add all the referenced buffers to the node.
    Node & bnode = result["buffers"];
    for (std::set<IndexType>::iterator s_it = buffer_indices.begin() ;
         s_it != buffer_indices.end() ; ++s_it)
    {
      // Use a dictionary layout here instead of conduit list.
      // Conduit IO HDF5 doesn't support conduit list objects.
      std::ostringstream oss;
      oss << "buffer_id_" << *s_it;
      Node& n_buffer = bnode.fetch( oss.str() );
      getDataStore()->getBuffer( *s_it )->exportTo(n_buffer);
    }
  }
}

/*
 *************************************************************************
 *
 * PRIVATE method to copy from Group to given Conduit node using
 * given set of ids to maintain correct association of data Buffers
 * to data Views.
 *
 *
 * Note: this is for the "sidre_hdf5" protocol
 *
 *************************************************************************
 */
void DataGroup::exportTo(conduit::Node& result,
                         std::set<IndexType>& buffer_indices) const
{
  if (getNumViews() > 0)
  {
    Node & vnode = result["views"];
    IndexType vidx = getFirstValidViewIndex();
    while ( indexIsValid(vidx) )
    {
      const DataView * view = getView(vidx);
      Node& n_view = vnode.fetch(view->getName());
      view->exportTo( n_view, buffer_indices );
      vidx = getNextValidViewIndex(vidx);
    }
  }

  if (getNumGroups() > 0)
  {
    Node & gnode = result["groups"];
    IndexType gidx = getFirstValidGroupIndex();
    while ( indexIsValid(gidx) )
    {
      const DataGroup * group =  getGroup(gidx);
      Node& n_group = gnode.fetch(group->getName());
      group->exportTo(n_group, buffer_indices);

      gidx = getNextValidGroupIndex(gidx);
    }
  }

  result.set(DataType::object());
}

/*
 *************************************************************************
 *
 * Imports tree from a conduit node into this DataGroup.  Includes
 * any Buffers attached to Views in that tree.
 *
 *
 * Note: this is for the "sidre_hdf5" protocol
 *
 *************************************************************************
 */

void DataGroup::importFrom(conduit::Node & node)
{
  // TODO - May want to put in a little meta-data into these files like a 'version'
  // or tag identifying the data.  We don't want someone giving us a file that
  // doesn't have our full multiView->buffer connectivity in there.

  destroyGroups();
  destroyViews();

  // First - Import Buffers into the DataStore.
  std::map<IndexType, IndexType> buffer_indices_map;

  if (node.has_path("buffers") )
  {
    conduit::NodeIterator buffs_itr = node["buffers"].children();
    while (buffs_itr.has_next())
    {
      Node& n_buffer = buffs_itr.next();
      IndexType old_buffer_id = n_buffer["id"].as_int32();

      DataBuffer * buffer = getDataStore()->createBuffer();

      // track change of old Buffer id to new Buffer id
      buffer_indices_map[ old_buffer_id ] = buffer->getIndex();

      // populate the new Buffer's state
      buffer->importFrom(n_buffer);
    }
  }

  // Next - import tree of Groups, sub-Groups, Views into the DataStore.
  // Use the mapping of old to new Buffer ids to connect the Views to the
  // right Buffers.
  importFrom(node, buffer_indices_map);

}


/*
 *************************************************************************
 *
 * PRIVATE method to copy from given Conduit node to this Group using
 * given map of ids to indicate association of Buffer ids in node to
 * those in datastore.
 *
 *
 * Note: this is for the "sidre_hdf5" protocol
 *
 *************************************************************************
 */
void DataGroup::importFrom(conduit::Node& node,
                           const std::map<IndexType, IndexType>& buffer_id_map)
{
  if ( node.has_path("views") )
  {
    // create the Views
    conduit::NodeIterator views_itr = node["views"].children();
    while (views_itr.has_next())
    {
      Node& n_view = views_itr.next();
      std::string view_name = views_itr.name();

      DataView * view = createView( view_name );
      view->importFrom(n_view, buffer_id_map);
    }
  }
  if ( node.has_path("groups") )
  {
    // create the child Groups
    conduit::NodeIterator groups_itr = node["groups"].children();
    while (groups_itr.has_next())
    {
      Node& n_group = groups_itr.next();
      std::string group_name = groups_itr.name();
      DataGroup * group = createGroup(group_name);
      group->importFrom(n_group, buffer_id_map);
    }
  }
}


/*
 *************************************************************************
 *
 * Imports tree from a conduit node into this DataGroup. 
 * This takes a generic conduit tree, not one with sidre conventions.
 *
 *************************************************************************
 */

void DataGroup::importConduitTree(conduit::Node &node)
{
  destroyGroups();
  destroyViews();

 //
  DataType node_dtype = node.dtype();
  if(node_dtype.is_object())
  {
      conduit::NodeIterator itr = node.children();
      while (itr.has_next())
      {
          Node&       cld_node  = itr.next();
          std::string cld_name  = itr.name();
          DataType    cld_dtype = cld_node.dtype();

          if(cld_dtype.is_object())
          {
              // create group
              DataGroup *grp = createGroup(cld_name);
              grp->importConduitTree(cld_node);
          }
          else if(cld_dtype.is_string())
          {
              //create string view
              createViewString(cld_name,cld_node.as_string());
          }
          else if(cld_dtype.is_number())
          {
             if(cld_dtype.number_of_elements() == 1)
             {
                  // create scalar view
                 DataView * view = createView(cld_name);
                 view->setScalar(cld_node);
             }
             else
             {
                 // create view with buffer
                 DataBuffer * buff = getDataStore()->createBuffer();
                 
                 conduit::index_t num_ele   = cld_dtype.number_of_elements();
                 conduit::index_t ele_bytes = DataType::default_bytes(cld_dtype.id());
                 
                 buff->allocate((TypeID)cld_dtype.id(),
                                num_ele);
                 // copy the data in a way that matches
                 // to compact representation of the buffer
                 conduit::uint8 * data_ptr = (conduit::uint8*) buff->getVoidPtr();
                 for(conduit::index_t i=0;i<num_ele;i++)
                 {
                     memcpy(data_ptr,
                            cld_node.element_ptr(i),
                            ele_bytes);
                     data_ptr+=ele_bytes;
                 }
                 
                 
                 DataView * view = createView(cld_name);
                 view->attachBuffer(buff);
                 // it is important to not use the data type directly
                 // it could contain offsets that are no longer 
                 // valid our new buffer
                 view->apply((TypeID)cld_dtype.id(),
                             cld_dtype.number_of_elements());
             }
          }
      }
  }
  else
  {
      SLIC_ERROR( "DataGroup cannot import non-object Conduit Node");
  }
  
}

/*
 *************************************************************************
 *
 * PRIVATE method to walk down a path to the next-to-last entry.
 *
 * If an error is encountered, this private function will return ATK_NULLPTR
 *
 *************************************************************************
 */
DataGroup * DataGroup::walkPath( std::string& path,
                                 bool create_groups_in_path )
{
  DataGroup * group_ptr = this;

  std::string::size_type pos = detail::find_exclusive( path, s_path_delimiter);
  if (pos != std::string::npos)
  {
    std::vector<std::string> tokens =
      detail::split(path, s_path_delimiter, pos);
    std::vector<std::string>::iterator stop = tokens.end() - 1;

    // Navigate path down to desired Group
    for (std::vector<std::string>::const_iterator iter = tokens.begin() ;
         iter < stop ; ++iter)
    {
      SLIC_ASSERT( iter->size() > 0 );

      if ( group_ptr->hasChildGroup(*iter) )
      {
        group_ptr = group_ptr->getGroup(*iter);
      }
      else if (create_groups_in_path)
      {
        group_ptr = group_ptr->createGroup(*iter);

        if ( group_ptr == ATK_NULLPTR )
        {
          iter = stop;
        }
      }
      else
      {
        iter = stop;
	group_ptr = ATK_NULLPTR;
      }
    }
    path = tokens.back();
  }

  return group_ptr;
}

/*
 *************************************************************************
 *
 * PRIVATE const method to walk down a path to the next-to-last entry.
 *
 * If an error is encountered, this private function will return ATK_NULLPTR
 *
 *************************************************************************
 */
const DataGroup * DataGroup::walkPath( std::string& path ) const
{
  const DataGroup * group_ptr = this;

  std::string::size_type pos = detail::find_exclusive( path, s_path_delimiter);
  if (pos != std::string::npos)
  {
    std::vector<std::string> tokens =
      detail::split(path, s_path_delimiter, pos);
    std::vector<std::string>::iterator stop = tokens.end() - 1;

    // Navigate path down to desired Group
    for (std::vector<std::string>::const_iterator iter = tokens.begin() ;
         iter < stop ; ++iter)
    {
      SLIC_ASSERT( iter->size() > 0 );

      if ( group_ptr->hasChildGroup(*iter) )
      {
        group_ptr = group_ptr->getGroup(*iter);
      }
      else
      {
	group_ptr = ATK_NULLPTR;
	iter = stop;
      }
    }
    path = tokens.back();
  }

  return group_ptr;
}



} /* end namespace sidre */
} /* end namespace asctoolkit */

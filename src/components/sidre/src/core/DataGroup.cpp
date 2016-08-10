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

  SLIC_CHECK_MSG( !path.empty() && group->hasChildGroup(intpath),
                  "Group " << getName() <<
                  " has no child Group with name '" << intpath << "'");

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
                  " has no child Group with name '" << path << "'");

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

  DataGroup * new_group = new(std::nothrow) DataGroup(intpath, this);
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
  //  n.reset();

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
  // n.reset();

  // Dump the group's views
  IndexType vidx = getFirstValidViewIndex();
  while ( indexIsValid(vidx) )
  {
    const DataView * view = getView(vidx);

    // Check that the view's name is not also a child group name
    SLIC_CHECK_MSG( !hasChildGroup(view->getName())
                    , view->getName() << " is the name of a groups and a view");

    view->createExternalLayout( n );
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
 * PRIVATE method to copy from Group to given Conduit node using
 * given set of ids to maintain correct association of data Buffers
 * to data Views.
 *
 *************************************************************************
 */

void DataGroup::exportTo(conduit::Node& data_holder,
                         std::set<IndexType>& buffer_indices) const
{
  if (getNumViews() > 0)
  {
    Node & vnode = data_holder["views"];
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
    Node & gnode = data_holder["groups"];
    IndexType gidx = getFirstValidGroupIndex();
    while ( indexIsValid(gidx) )
    {
      const DataGroup * group =  getGroup(gidx);
      Node& n_group = gnode.fetch(group->getName());
      group->exportTo(n_group, buffer_indices);

      gidx = getNextValidGroupIndex(gidx);
    }
  }

  // TODO - take this out when CON-131 resolved ( can't write out empty node ).
  if (data_holder.dtype().is_empty() )
  {
    data_holder.set_string("empty");
  }

}

/*
 *************************************************************************
 *
 * PRIVATE method to copy from given Conduit node to this Group using
 * given map of ids to indicate association of Buffer ids in node to
 * those in datastore.
 *
 *************************************************************************
 */
void DataGroup::importFrom(conduit::Node& data_holder,
                           const std::map<IndexType, IndexType>& buffer_id_map)
{
// If the Group is empty, conduit will complain if you call 'has_path'.

  // Added CON-132 ticket asking if has_path can just return false if node
  // is empty or not an object type.
  if ( data_holder.dtype().is_object() && data_holder.has_path("views") )
  {
    // create the Views
    conduit::NodeIterator views_itr = data_holder["views"].children();
    while (views_itr.has_next())
    {
      Node& n_view = views_itr.next();
      std::string view_name = views_itr.path();

      DataView * view = createView( view_name );
      view->importFrom(n_view, buffer_id_map);
    }
  }
  if ( data_holder.dtype().is_object() && data_holder.has_path("groups") )
  {
    // create the child Groups
    conduit::NodeIterator groups_itr = data_holder["groups"].children();
    while (groups_itr.has_next())
    {
      Node& n_group = groups_itr.next();
      std::string group_name = groups_itr.path();
      DataGroup * group = createGroup(group_name);
      group->importFrom(n_group, buffer_id_map);
    }
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

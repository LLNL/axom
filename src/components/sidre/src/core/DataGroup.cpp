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
// View access methods.
//
////////////////////////////////////////////////////////////////////////

/*
 *************************************************************************
 *
 * Return pointer to non-const view with given name or path if it exists.
 *
 *************************************************************************
 */
DataView * DataGroup::getView( const std::string& name )
{
  std::string path = name;
  DataGroup * group = walkPath( path, false );

  SLIC_CHECK_MSG( !path.empty() && group->hasView(name),
                  "Group " << getName() << 
                  " has no view with name '" << path << "'");

  return group->m_view_coll.getItem(path);
}

/*
 *************************************************************************
 *
 * Return pointer to const view with given name or path if it exists.
 *
 *************************************************************************
 */
const DataView * DataGroup::getView( const std::string& name ) const
{
// XXXX: Add path implementation
  SLIC_CHECK_MSG( !name.empty() && hasView(name),
                  "Group " << getName() << 
                  " has no view with name '" << name << "'");

  return m_view_coll.getItem(name);
}



////////////////////////////////////////////////////////////////////////
//
//  View creation methods.
//
////////////////////////////////////////////////////////////////////////

/*
 *************************************************************************
 *
 * Create view with given name, but no data description.
 *
 *************************************************************************
 */
DataView * DataGroup::createView( const std::string& name )
{
  std::string path = name;
  DataGroup * group = walkPath( path, true );

  if ( group == ATK_NULLPTR )
  {
    SLIC_CHECK( group != ATK_NULLPTR );
    return ATK_NULLPTR;
  }
  else if ( path.empty() || group->hasView(path) )
  {
    SLIC_CHECK( !path.empty() );
    SLIC_CHECK( !group->hasView(path) );
    return ATK_NULLPTR;
  }

  // Want the C++ new operator to return null pointer on failure instead of
  // throwing an exception.
  DataView * view = new(std::nothrow) DataView(path);
  if ( view != ATK_NULLPTR )
  {
    group->attachView(view);
  }
  return view;
}

/*
 *************************************************************************
 *
 * Create view with given name, and data description using type and
 * number of elements.
 *
 *************************************************************************
 */
DataView * DataGroup::createView( const std::string& name,
                                  TypeID type,
                                  SidreLength num_elems )
{
  if ( type == NO_TYPE_ID || num_elems < 0 )
  {
    SLIC_CHECK(type != NO_TYPE_ID);
    SLIC_CHECK_MSG(num_elems >= 0,
                   "Must define view with number of elems >=0 ");
    return ATK_NULLPTR;
  }

  DataView * view = createView(name);
  if (view != ATK_NULLPTR)
  {
    view->describe(type, num_elems);
  }
  return view;
}

/*
 *************************************************************************
 *
 * Create view with given name, and data description using type and shape.
 *
 *************************************************************************
 */
DataView * DataGroup::createView( const std::string& name,
                                  TypeID type,
                                  int ndims,
                                  SidreLength * shape )
{
  if ( type == NO_TYPE_ID || !(ndims >= 0) )
  {
    SLIC_CHECK(type != NO_TYPE_ID);
    SLIC_CHECK_MSG(ndims >= 0,
                   "Must define view with number of ndims >=0 ");
    return ATK_NULLPTR;
  }

  DataView * view = createView(name);
  if (view != ATK_NULLPTR)
  {
    view->describe(type, ndims, shape);
  }
  return view;
}

/*
 *************************************************************************
 *
 * Create view with given name, and data description using a conduit
 * Datatype class.
 *
 *************************************************************************
 */
DataView * DataGroup::createView( const std::string& name,
                                  const DataType& dtype )
{
  DataView * view = createView(name);
  if (view != ATK_NULLPTR)
  {
    view->describe(dtype);
  }

  return view;
}

/*
 *************************************************************************
 *
 * Create view with given name, but no data description.  Attach provided
 * buffer to view.
 *
 *************************************************************************
 */
DataView * DataGroup::createView( const std::string& name,
                                  DataBuffer * buff)
{
  DataView * view = createView(name);
  if ( view != ATK_NULLPTR )
  {
    view->attachBuffer( buff );
  }
  return view;
}

/*
 *************************************************************************
 *
 * Create view with given name, pointing to undescribed external data.
 *
 *************************************************************************
 */
DataView * DataGroup::createView( const std::string& name,
                                  void * external_ptr )
{
  DataView * view = createView(name);
  if ( view != ATK_NULLPTR )
  {
    view->setExternalDataPtr(external_ptr);
  }
  return view;
}

/*
 *************************************************************************
 *
 * Create view with given name, and data description using type and
 * number of elements.
 *
 * In addition, create an associated buffer, allocate the data, and attach
 * view to new buffer.
 *
 *************************************************************************
 */
DataView * DataGroup::createViewAndAllocate( const std::string& name,
                                             TypeID type,
                                             SidreLength num_elems )
{
  // createView will verify args.
  DataView * view = createView(name, type, num_elems);
  if ( view != ATK_NULLPTR )
  {
    view->allocate();
  }
  return view;
}

/*
 *************************************************************************
 *
 * Create view with given name, and data description using type,
 * number of dimnesions, and number of elements per dimension.
 *
 * In addition, create an associated buffer, allocate the data, and attach
 * view to new buffer.
 *
 *************************************************************************
 */
DataView * DataGroup::createViewAndAllocate( const std::string& name,
                                             TypeID type,
                                             int ndims,
                                             SidreLength * shape )
{
  // createView will verify args.
  DataView * view = createView(name, type, ndims, shape);
  if ( view != ATK_NULLPTR )
  {
    view->allocate();
  }
  return view;
}

/*
 *************************************************************************
 *
 * Create view with given name, and data description using a conduit
 * Datatype class.
 *
 * In addition, create an associated buffer, allocate the data, and attach
 * view to new buffer.
 *
 *************************************************************************
 */
DataView * DataGroup::createViewAndAllocate( const std::string& name,
                                             const DataType& dtype)
{
  // createView will verify args.
  DataView * view = createView(name, dtype);
  if ( view != ATK_NULLPTR )
  {
    view->allocate();
  }
  return view;
}

/*
 *************************************************************************
 *
 * Create view with given name for a string.
 *
 *************************************************************************
 */
DataView * DataGroup::createViewString( const std::string& name,
                                        const std::string& value)
{
  DataView * view = createView(name);
  if (view != ATK_NULLPTR)
  {
    view->setString(value);
  }

  return view;
}

////////////////////////////////////////////////////////////////////////
//
//  Methods for deleting, detaching, copying, or moving DataView object
//  in a DataGroup
//
////////////////////////////////////////////////////////////////////////


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
 * Detach view from group and destroy view.
 *
 * Destroy DataBuffer if buffer has no other view attached to it.
 *
 *************************************************************************
 */
void DataGroup::destroyViewAndData( DataView * view )
{
  if ( view != ATK_NULLPTR)
  {
    detachView( view->getName() );
    DataBuffer * const buffer = view->detachBuffer();
    if ( buffer != ATK_NULLPTR )
    {
      if (buffer->getNumViews() == 0)
      {
        getDataStore()->destroyBuffer(buffer);
      }
    }
    delete view;
  }
}

/*
 *************************************************************************
 *
 * Detach view with given name from group and destroy view.
 *
 * Destroy DataBuffer if buffer has no other view attached to it.
 *
 *************************************************************************
 */
void DataGroup::destroyViewAndData( const std::string& name )
{
  destroyViewAndData(getView(name));
}

/*
 *************************************************************************
 *
 * Detach view with given index from group and destroy view.
 *
 * Destroy DataBuffer if buffer has no other view attached to it.
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
    destroyViewAndData(vidx);
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
  if ( view == ATK_NULLPTR )
  {
    SLIC_CHECK( view != ATK_NULLPTR );
    return ATK_NULLPTR;
  }

  DataGroup * curr_group = view->getOwningGroup();
  if (curr_group == this)
  {
    // this group already owns the view
    return view;
  }
  else if (hasView(view->getName()))
  {
    SLIC_CHECK_MSG(!hasView(view->getName()),
                   "Group '" << getName() <<
                   "' already has a view named'" << view->getName() <<
                   "' so view move operation cannot happen");
    return ATK_NULLPTR;
  }

  curr_group->detachView(view);
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
DataView * DataGroup::copyView(DataView * view)
{
  if ( view == ATK_NULLPTR || hasView(view->getName()) )
  {
    SLIC_CHECK( view != ATK_NULLPTR );
    SLIC_CHECK_MSG(!hasView(view->getName()),
                   "Group '" << getName() <<
                   "' already has a view named'" << view->getName() <<
                   "' so view copy operation cannot happen");

    return ATK_NULLPTR;
  }

  DataView * copy = createView(view->getName());
  view->copyView(copy);
  return copy;
}



////////////////////////////////////////////////////////////////////////
//
// Child group access methods.
//
////////////////////////////////////////////////////////////////////////

/*
 *************************************************************************
 *
 * Return pointer to non-const child group with given name or path 
 * if it exists.
 *
 *************************************************************************
 */
DataGroup * DataGroup::getGroup( const std::string& name )
{
  std::string path = name;
  DataGroup * group = walkPath( path, false );

  SLIC_CHECK_MSG( !path.empty() && group->hasGroup(name),
                  "Group " << getName() << 
                  " has no child group with name '" << path << "'");

  return group->m_group_coll.getItem(path);
}

/*
 *************************************************************************
 *
 * Return pointer to const child group with given name or path if it exists.
 *
 *************************************************************************
 */
const DataGroup * DataGroup::getGroup( const std::string& name ) const
{
// XXXX: Add path implementation
  SLIC_CHECK_MSG( !name.empty() && hasGroup(name),
                  "Group " << getName() << 
                  " has no child group with name '" << name << "'");

  return m_group_coll.getItem(name);
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

  std::string path = name;
  DataGroup * group = walkPath( path, true );

  if ( group == ATK_NULLPTR )
  {
    SLIC_CHECK( group != ATK_NULLPTR );
    return ATK_NULLPTR;
  }
  else if ( path.empty() || group->hasGroup(name) )
  {
    SLIC_CHECK( !path.empty() );
    SLIC_CHECK_MSG( !group->hasGroup(name), "name == " << name );

    return ATK_NULLPTR;
  }

  DataGroup * new_group = new(std::nothrow) DataGroup( name, this);
  if ( new_group == ATK_NULLPTR )
  {
    return ATK_NULLPTR;
  }
  return group->attachGroup(new_group);
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
  if ( group == ATK_NULLPTR || hasGroup(group->getName()))
  {
    SLIC_CHECK( group != ATK_NULLPTR );
    SLIC_CHECK_MSG(!hasGroup(group->getName()),
                   "Group '" << getName() << 
                   "' already has a child group named'" << group->getName() <<
                   "' so group move operation cannot happen");

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
 * Create a copy of given group and make it a child of this group.
 *
 * Copying a group does not perform a deep copy of any of its buffers.
 *
 *************************************************************************
 */
DataGroup * DataGroup::copyGroup(DataGroup * group)
{
  if ( group == ATK_NULLPTR || hasGroup(group->getName()) )
  {
    SLIC_CHECK( group != ATK_NULLPTR );
    SLIC_CHECK_MSG(!hasGroup(group->getName()),
                   "Group '" << getName() << 
                   "' already has a child group named'" << group->getName() <<
                   "' so group copy operation cannot happen");

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
 * Test this DataGroup for equavalence to another DataGroup.
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

  // Test equivalence of DataViews
  if (is_equiv)
  {
    IndexType vidx = getFirstValidViewIndex();
    IndexType other_vidx = other->getFirstValidViewIndex();
    while ( is_equiv && indexIsValid(vidx) && indexIsValid(other_vidx) )
    {
      const DataView * view = getView(vidx);
      const DataView * other_view = other->getView(other_vidx);
      is_equiv = view->isEquivalentTo(other_view);
      vidx = getNextValidViewIndex(vidx);
      other_vidx = getNextValidViewIndex(other_vidx);
    }
  }

  // Recursively call this method to test equivalence of child DataGroups
  if (is_equiv)
  {
    IndexType gidx = getFirstValidGroupIndex();
    IndexType other_gidx = getFirstValidGroupIndex();
    while ( is_equiv && indexIsValid(gidx) && indexIsValid(other_gidx) )
    {
      const DataGroup * group =  getGroup(gidx);
      const DataGroup * other_group =  other->getGroup(other_gidx);
      is_equiv = group->isEquivalentTo(other_group);
      gidx = getNextValidGroupIndex(gidx);
      other_gidx = getNextValidGroupIndex(other_gidx);
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
{}

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
 * PRIVATE method to walk down a path to the next-to-last entry.
 *
 * If an empty path or a path with just a single entry, 'foo', is passed
 * in, then this routine will simply return the current group.
 *
 * If an error is encoutered, this private function will return ATK_NULLPTR
 *
 * 'path': The path to traverse.  Traversed entries will be stripped
 * from this parameter at algorithm completion, leaving the last entry
 * in the path remaining.
 *
 * 'create_on_demand': If true, will create any groups that are not
 * found while walking the path.
 *
 *************************************************************************
 */

// Developer notes:
// At next pass, should optimize this.
// #1 the split routine performs string copies when populating a new
//   vector of the tokens.  That could be rewritten to instead populate a
//   vector of delimiter positions and avoid string copying.
// #2 if the provided path is invalid, our code will halt.  This should be improved
//   to handle it.  (Figure out what an appropriate return value is, and how
//   caller public function should handle it).
// #3 Since this function can create groups on demand, it's not a const function so
//   I left it out of our const versions of getGroup, getView.  Need to just split
//   this function into two, one const and one non-const.
// -- AB

DataGroup * DataGroup::walkPath( std::string& path, bool create_on_demand )
{
  // TODO - write tests that pass in some error conditions and verify code for:
  // path = ""
  // path = "foo"
  // these should result in just getting back the same group you just called this
  // routine from and issue a warning, but not cause a code crash.

  DataGroup * group_ptr = this;

  std::string::size_type pos = detail::find_exclusive( path, s_path_delimiter);
  if (pos != std::string::npos)
  {
    std::vector<std::string> tokens =
      detail::split(path, s_path_delimiter, pos);
    std::vector<std::string>::iterator stop = tokens.end() - 1;

    // Navigate path down to desired group.
    for (std::vector<std::string>::const_iterator iter = tokens.begin() ;
         iter < stop ; ++iter)
    {
      SLIC_ASSERT( iter->size() > 0 );

      if ( group_ptr->hasGroup(*iter) )
      {
        group_ptr = group_ptr->getGroup(*iter);
      }
      else if (create_on_demand)
      {
        group_ptr = group_ptr->createGroup(*iter);

        if ( group_ptr == ATK_NULLPTR )
        {
          iter = stop;
        }
      }
      else
      {
//        SLIC_CHECK_MSG(false, "Path is invalid, group " << group_ptr->getName() << " does not have group with name " << *iter);
        SLIC_ERROR(
          "Path is invalid, group " << group_ptr->getName() << " does not have group with name " <<
          *iter);
      }
    }
    path = tokens.back();
  }

  return group_ptr;
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
    SLIC_ASSERT(view->m_owning_group == ATK_NULLPTR);
    view->m_owning_group = this;
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
 * PRIVATE method to copy from group to given Conduit node using
 * given set of ids to maintain correct association of data buffers
 * to data views.
 *
 *************************************************************************
 */

void DataGroup::exportTo(conduit::Node& data_holder,
                         std::set<IndexType>& buffer_indices) const
{
  IndexType vidx = getFirstValidViewIndex();
  while ( indexIsValid(vidx) )
  {
    const DataView * view = getView(vidx);
    Node& n_view = data_holder["views"].fetch(view->getName());
    view->exportTo( n_view, buffer_indices );
    vidx = getNextValidViewIndex(vidx);
  }

  IndexType gidx = getFirstValidGroupIndex();
  while ( indexIsValid(gidx) )
  {
    const DataGroup * group =  getGroup(gidx);
    Node& n_group = data_holder["groups"].fetch(group->getName());
    group->exportTo(n_group, buffer_indices);

    gidx = getNextValidGroupIndex(gidx);
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
 * PRIVATE method to copy from given Conduit node to this group using
 * given map of ids to indicate association of buffer ids in node to
 * those in datastore.
 *
 *************************************************************************
 */
void DataGroup::importFrom(conduit::Node& data_holder,
                           const std::map<IndexType, IndexType>& buffer_id_map)
{
// If the group is empty, conduit will complain if you call 'has_path'.

  // Added CON-132 ticket asking if has_path can just return false if node is empty or not an object type.
  if ( data_holder.dtype().is_object() && data_holder.has_path("views") )
  {
    // create the views
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
    // create the child groups
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



} /* end namespace sidre */
} /* end namespace asctoolkit */

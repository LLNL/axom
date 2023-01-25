// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Associated header file
#include "Group.hpp"

#include "conduit_relay.hpp"

#ifdef AXOM_USE_HDF5
  #include "conduit_relay_io_hdf5.hpp"
#endif

#include "axom/core/Macros.hpp"
#include "axom/core/Path.hpp"

// Sidre headers
#include "ListCollection.hpp"
#include "MapCollection.hpp"
#include "Buffer.hpp"
#include "DataStore.hpp"

namespace axom
{
namespace sidre
{
// Helper macro for defining a prepend string for sidre::Group log messages
// We are using it to add the pathName() of the group
#ifndef SIDRE_GROUP_LOG_PREPEND
  #define SIDRE_GROUP_LOG_PREPEND             \
    "[Group: '" << this->getPathName() << "'" \
                << (this->isRoot() ? " (root)" : "") << "] "
#endif

// Initialization of static path delimiter character for methods that
// support path syntax.
const char Group::s_path_delimiter = '/';

///////////////////////////////////////////////////////////////////////////////
//
// Private utility functions to cast ItemCollections to (named) MapCollections.
//
///////////////////////////////////////////////////////////////////////////////

MapCollection<View>* Group::getNamedViews()
{
  SLIC_ASSERT_MSG(this->isUsingMap(),
                  "Invalid cast: The views in this group do not have names");

  return static_cast<MapCollection<View>*>(m_view_coll);
}

const MapCollection<View>* Group::getNamedViews() const
{
  SLIC_ASSERT_MSG(this->isUsingMap(),
                  "Invalid cast: The views in this group do not have names");

  return static_cast<const MapCollection<View>*>(m_view_coll);
}

MapCollection<Group>* Group::getNamedGroups()
{
  SLIC_ASSERT_MSG(this->isUsingMap(),
                  "Invalid cast: The groups in this group do not have names");

  return static_cast<MapCollection<Group>*>(m_group_coll);
}

const MapCollection<Group>* Group::getNamedGroups() const
{
  SLIC_ASSERT_MSG(this->isUsingMap(),
                  "Invalid cast: The groups in this group do not have names");

  return static_cast<const MapCollection<Group>*>(m_group_coll);
}

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
std::string Group::getPath() const
{
  const Group* root = getDataStore()->getRoot();
  const Group* curr = getParent();
  std::string thePath = curr->getName();
  curr = curr->getParent();

  while(curr != root)
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
bool Group::hasView(const std::string& path) const
{
  std::string intpath(path);
  const Group* group = walkPath(intpath);

  if(group == nullptr)
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
View* Group::getView(const std::string& path)
{
  std::string intpath(path);
  bool create_groups_in_path = false;
  Group* group = walkPath(intpath, create_groups_in_path);

  if(group == nullptr)
  {
    SLIC_CHECK_MSG(
      group != nullptr,
      SIDRE_GROUP_LOG_PREPEND << "Non-existent group in path " << path);
    return nullptr;
  }

  SLIC_CHECK_MSG(
    !intpath.empty() && group->hasChildView(intpath),
    SIDRE_GROUP_LOG_PREPEND << "No View with name '" << intpath << "'");

  return group->getNamedViews()->getItem(intpath);
}

/*
 *************************************************************************
 *
 * Return pointer to const View with given name or path if it exists.
 *
 *************************************************************************
 */
const View* Group::getView(const std::string& path) const
{
  std::string intpath(path);
  const Group* group = walkPath(intpath);

  if(group == nullptr)
  {
    SLIC_CHECK_MSG(
      group != nullptr,
      SIDRE_GROUP_LOG_PREPEND << "Non-existent group in path " << path);
    return nullptr;
  }

  SLIC_CHECK_MSG(
    !intpath.empty() && group->hasChildView(intpath),
    SIDRE_GROUP_LOG_PREPEND << "No View with name '" << intpath << "'");

  return group->getNamedViews()->getItem(intpath);
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
View* Group::createView(const std::string& path)
{
  std::string intpath(path);

  Group* group;
  if(intpath.empty())
  {
    SLIC_CHECK_MSG(m_is_list,
                   SIDRE_GROUP_LOG_PREPEND
                     << "Could not create View with empty string "
                     << "for the path.");
    if(m_is_list)
    {
      group = this;
    }
    else
    {
      return nullptr;
    }
  }
  else
  {
    bool create_groups_in_path = true;
    group = walkPath(intpath, create_groups_in_path);

    if(group == nullptr)
    {
      if(m_is_list)
      {
        SLIC_CHECK_MSG(group != nullptr,
                       SIDRE_GROUP_LOG_PREPEND
                         << "Could not find or create View '" << path << "'."
                         << "for a Group using the list format.");
      }
      else
      {
        SLIC_CHECK_MSG(group != nullptr,
                       SIDRE_GROUP_LOG_PREPEND
                         << "Could not find or create path '" << path << "'."
                         << "There is already a View with that name.");
      }
      return nullptr;
    }
    else if(intpath.empty() || group->hasChildView(intpath) ||
            group->hasChildGroup(intpath))
    {
      SLIC_CHECK_MSG(
        !intpath.empty(),
        SIDRE_GROUP_LOG_PREPEND << "Cannot create a View with an empty path.");
      SLIC_CHECK_MSG(!group->hasChildView(intpath),
                     SIDRE_GROUP_LOG_PREPEND
                       << "Cannot create View with name '" << intpath << "'. "
                       << "There is already a View with that name.");
      SLIC_CHECK_MSG(!group->hasChildGroup(intpath),
                     SIDRE_GROUP_LOG_PREPEND
                       << "Cannot create View with name '" << intpath << "'. "
                       << "There is already has a Group with that name.");
      return nullptr;
    }
  }

  View* view = new(std::nothrow) View(intpath);
  if(view != nullptr)
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
View* Group::createView(const std::string& path, TypeID type, IndexType num_elems)
{
  if(type == NO_TYPE_ID || num_elems < 0)
  {
    SLIC_CHECK_MSG(type != NO_TYPE_ID,
                   SIDRE_GROUP_LOG_PREPEND << "Cannot create View with name '"
                                           << path << "'. "
                                           << " Invalid type << " << type);
    SLIC_CHECK_MSG(num_elems >= 0,
                   SIDRE_GROUP_LOG_PREPEND
                     << "Cannot create View with name '" << path << "'. "
                     << "Number of elements cannot be less than zero.");
    return nullptr;
  }

  View* view = createView(path);
  if(view != nullptr)
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
View* Group::createViewWithShape(const std::string& path,
                                 TypeID type,
                                 int ndims,
                                 const IndexType* shape)
{
  if(type == NO_TYPE_ID || ndims < 0 || shape == nullptr)
  {
    SLIC_CHECK_MSG(type != NO_TYPE_ID,
                   SIDRE_GROUP_LOG_PREPEND
                     << "Problem creating View with name '" << path << "'. "
                     << " Invalid type: " << type);
    SLIC_CHECK_MSG(ndims >= 0,
                   SIDRE_GROUP_LOG_PREPEND
                     << "Problem creating View with name '" << path << "'. "
                     << "ndims must be greater than 0.");
    SLIC_CHECK_MSG(shape != nullptr,
                   SIDRE_GROUP_LOG_PREPEND
                     << "Problem creating View with name '" << path << "'. "
                     << "shape pointer was null.");
    return nullptr;
  }

  View* view = createView(path);
  if(view != nullptr)
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
View* Group::createView(const std::string& path, const DataType& dtype)
{
  View* view = createView(path);
  if(view != nullptr)
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
View* Group::createView(const std::string& path, Buffer* buff)
{
  View* view = createView(path);
  if(view != nullptr)
  {
    view->attachBuffer(buff);
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
View* Group::createView(const std::string& path,
                        TypeID type,
                        IndexType num_elems,
                        Buffer* buff)
{
  View* view = createView(path, type, num_elems);
  if(view != nullptr)
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
View* Group::createViewWithShape(const std::string& path,
                                 TypeID type,
                                 int ndims,
                                 const IndexType* shape,
                                 Buffer* buff)
{
  View* view = createViewWithShape(path, type, ndims, shape);
  if(view != nullptr)
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
View* Group::createView(const std::string& path, const DataType& dtype, Buffer* buff)
{
  View* view = createView(path, dtype);
  if(view != nullptr)
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
View* Group::createView(const std::string& path, void* external_ptr)
{
  View* view = createView(path);
  if(view != nullptr)
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
View* Group::createView(const std::string& path,
                        TypeID type,
                        IndexType num_elems,
                        void* external_ptr)
{
  View* view = createView(path, type, num_elems);
  if(view != nullptr)
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
View* Group::createViewWithShape(const std::string& path,
                                 TypeID type,
                                 int ndims,
                                 const IndexType* shape,
                                 void* external_ptr)
{
  View* view = createViewWithShape(path, type, ndims, shape);
  if(view != nullptr)
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
View* Group::createView(const std::string& path,
                        const DataType& dtype,
                        void* external_ptr)
{
  View* view = createView(path, dtype);
  if(view != nullptr)
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
View* Group::createViewAndAllocate(const std::string& path,
                                   TypeID type,
                                   IndexType num_elems,
                                   int allocID)
{
  allocID = getValidAllocatorID(allocID);

  View* view = createView(path, type, num_elems);
  if(view != nullptr)
  {
    view->allocate(allocID);
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
View* Group::createViewWithShapeAndAllocate(const std::string& path,
                                            TypeID type,
                                            int ndims,
                                            const IndexType* shape,
                                            int allocID)
{
  allocID = getValidAllocatorID(allocID);

  View* view = createViewWithShape(path, type, ndims, shape);
  if(view != nullptr)
  {
    view->allocate(allocID);
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
View* Group::createViewAndAllocate(const std::string& path,
                                   const DataType& dtype,
                                   int allocID)
{
  allocID = getValidAllocatorID(allocID);

  View* view = createView(path, dtype);
  if(view != nullptr)
  {
    view->allocate(allocID);
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
View* Group::createViewString(const std::string& path, const std::string& value)
{
  View* view = createView(path);
  if(view != nullptr)
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
void Group::destroyView(const std::string& path)
{
  std::string intpath(path);
  bool create_groups_in_path = false;
  Group* group = walkPath(intpath, create_groups_in_path);

  if(group != nullptr)
  {
    View* view = group->detachView(intpath);
    if(view != nullptr)
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
void Group::destroyView(IndexType idx)
{
  View* view = detachView(idx);
  if(view != nullptr)
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
void Group::destroyViews()
{
  IndexType vidx = getFirstValidViewIndex();
  while(indexIsValid(vidx))
  {
    View* view = detachView(vidx);
    if(view != nullptr)
    {
      delete view;
    }

    vidx = getNextValidViewIndex(vidx);
  }

  m_view_coll->removeAllItems();
}

/*
 *************************************************************************
 *
 * Destroy View with given name or path and its data if it's the only View
 * associated with that data.
 *
 *************************************************************************
 */
void Group::destroyViewAndData(const std::string& path)
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
void Group::destroyViewAndData(IndexType idx)
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
void Group::destroyViewsAndData()
{
  IndexType vidx = getFirstValidViewIndex();
  while(indexIsValid(vidx))
  {
    destroyViewAndData(vidx);
    vidx = getNextValidViewIndex(vidx);
  }

  m_view_coll->removeAllItems();
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
View* Group::moveView(View* view)
{
  if(view == nullptr)
  {
    SLIC_CHECK_MSG(
      view != nullptr,
      SIDRE_GROUP_LOG_PREPEND << "Null pointer provided, no View to move.");
    return nullptr;
  }

  Group* curr_group = view->getOwningGroup();
  if(curr_group == this)
  {
    // this Group already owns the View
    return view;
  }
  else if(hasChildView(view->getName()))
  {
    SLIC_CHECK_MSG(!hasChildView(view->getName()),
                   SIDRE_GROUP_LOG_PREPEND
                     << "Group already has a View named '" << view->getName()
                     << "' so View move operation cannot happen");
    return nullptr;
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
View* Group::copyView(View* view)
{
  if(view == nullptr || hasChildView(view->getName()))
  {
    SLIC_CHECK_MSG(
      view != nullptr,
      SIDRE_GROUP_LOG_PREPEND << "Null pointer provided, no View to copy.");

    if(view != nullptr)
    {
      SLIC_CHECK_MSG(!hasChildView(view->getName()),
                     SIDRE_GROUP_LOG_PREPEND
                       << "Group already has a View named '" << view->getName()
                       << "' so View copy operation cannot happen");
    }

    return nullptr;
  }

  View* copy = createView(view->getName());
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
bool Group::hasGroup(const std::string& path) const
{
  std::string intpath(path);
  const Group* group = walkPath(intpath);

  if(group == nullptr)
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
Group* Group::getGroup(const std::string& path)
{
  std::string intpath(path);
  bool create_groups_in_path = false;
  Group* group = walkPath(intpath, create_groups_in_path);

  if(group == nullptr)
  {
    SLIC_CHECK_MSG(group != nullptr,
                   SIDRE_GROUP_LOG_PREPEND << "Non-existent group in path '"
                                           << path << "'.");
    return nullptr;
  }

  SLIC_CHECK_MSG(!intpath.empty() && group->hasChildGroup(intpath),
                 SIDRE_GROUP_LOG_PREPEND
                   << "Group has no descendant Group named '" << path << "'.");

  return group->getNamedGroups()->getItem(intpath);
}

/*
 *************************************************************************
 *
 * Return pointer to const child Group with given name or path if it exists.
 *
 *************************************************************************
 */
const Group* Group::getGroup(const std::string& path) const
{
  std::string intpath(path);
  const Group* group = walkPath(intpath);

  if(group == nullptr)
  {
    SLIC_CHECK_MSG(
      group != nullptr,
      SIDRE_GROUP_LOG_PREPEND << "Non-existent group in path " << path);
    return nullptr;
  }

  SLIC_CHECK_MSG(!intpath.empty() && group->hasChildGroup(intpath),
                 SIDRE_GROUP_LOG_PREPEND
                   << "Group has no descendant Group with name '" << path
                   << "'.");

  return group->getNamedGroups()->getItem(intpath);
}

////////////////////////////////////////////////////////////////////////
//
//  Methods for managing child Group objects in Group
//
////////////////////////////////////////////////////////////////////////

/*
 *************************************************************************
 *
 * Create Group with given name or path and make it a child of this Group.
 *
 *************************************************************************
 */
Group* Group::createGroup(const std::string& path, bool is_list)
{
  std::string intpath(path);
  bool create_groups_in_path = true;
  Group* group = walkPath(intpath, create_groups_in_path);

  if(group == nullptr)
  {
    if(m_is_list)
    {
      SLIC_CHECK_MSG(group != nullptr,
                     SIDRE_GROUP_LOG_PREPEND
                       << "Could not find or create Group '" << path << "'."
                       << "for a Group using the list format.");
    }
    else
    {
      SLIC_CHECK_MSG(group != nullptr,
                     SIDRE_GROUP_LOG_PREPEND
                       << "Could not find or create path '" << path
                       << "'. There is already a Group with that name");
    }
    return nullptr;
  }
  else if(intpath.empty() || group->hasChildGroup(intpath) ||
          group->hasChildView(intpath))
  {
    SLIC_CHECK_MSG(
      !intpath.empty(),
      SIDRE_GROUP_LOG_PREPEND << "Cannot create a Group with an empty path.");
    SLIC_CHECK_MSG(!group->hasChildGroup(intpath),
                   SIDRE_GROUP_LOG_PREPEND
                     << "Cannot create Group with name '" << path
                     << "'. There is already has a Group with that name.");
    SLIC_CHECK_MSG(!group->hasChildView(intpath),
                   SIDRE_GROUP_LOG_PREPEND
                     << "Cannot create Group with name '" << path
                     << "'. There is already has a View with that name.");

    return nullptr;
  }

  Group* new_group =
    new(std::nothrow) Group(intpath, group->getDataStore(), is_list);
  if(new_group == nullptr)
  {
    return nullptr;
  }

#ifdef AXOM_USE_UMPIRE
  new_group->setDefaultAllocator(group->getDefaultAllocator());
#endif
  return group->attachGroup(new_group);
}

Group* Group::createUnnamedGroup(bool is_list)
{
  SLIC_CHECK_MSG(m_is_list,
                 SIDRE_GROUP_LOG_PREPEND
                   << "Cannot create an unnamed Group when not using "
                   << "list format.");

  Group* new_group;
  if(m_is_list)
  {
    new_group = new(std::nothrow) Group("", getDataStore(), is_list);
  }
  else
  {
    new_group = nullptr;
  }

  if(new_group == nullptr)
  {
    return nullptr;
  }

#ifdef AXOM_USE_UMPIRE
  new_group->setDefaultAllocator(getDefaultAllocator());
#endif
  return attachGroup(new_group);
}

/*
 *************************************************************************
 *
 * Detach child Group with given name or path and destroy it.
 *
 *************************************************************************
 */
void Group::destroyGroup(const std::string& path)
{
  std::string intpath(path);
  bool create_groups_in_path = false;
  Group* group = walkPath(intpath, create_groups_in_path);

  if(group != nullptr)
  {
    Group* targetgroup = group->detachGroup(intpath);
    if(targetgroup != nullptr)
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
void Group::destroyGroup(IndexType idx)
{
  Group* group = detachGroup(idx);
  if(group != nullptr)
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
void Group::destroyGroups()
{
  IndexType gidx = getFirstValidGroupIndex();
  while(indexIsValid(gidx))
  {
    Group* group = this->getGroup(gidx);
    delete group;

    gidx = getNextValidGroupIndex(gidx);
  }

  m_group_coll->removeAllItems();
}

/*
 *************************************************************************
 *
 * Remove given Group from its owning Group and make it a child of this Group.
 *
 *************************************************************************
 */
Group* Group::moveGroup(Group* group)
{
  if(group == nullptr || hasChildGroup(group->getName()))
  {
    SLIC_CHECK_MSG(
      group != nullptr,
      SIDRE_GROUP_LOG_PREPEND << "Null pointer provided, no Group to move.");

    if(group != nullptr)
    {
      SLIC_CHECK_MSG(!hasChildGroup(group->getName()),
                     SIDRE_GROUP_LOG_PREPEND
                       << "Invalid move operation. Group already has "
                       << "a child named '" << group->getName() << "'.");
    }

    return nullptr;
  }

  Group* curr_group = group->getParent();
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
Group* Group::copyGroup(Group* group)
{
  if(group == nullptr || hasChildGroup(group->getName()))
  {
    SLIC_CHECK_MSG(
      group != nullptr,
      SIDRE_GROUP_LOG_PREPEND << "Null pointer provided, no Group to copy.");

    if(group != nullptr)
    {
      SLIC_CHECK_MSG(!hasChildGroup(group->getName()),
                     SIDRE_GROUP_LOG_PREPEND
                       << "Invalid copy operation. Group already has "
                       << "a child named '" << group->getName() << "'.");
    }

    return nullptr;
  }

  Group* res = createGroup(group->getName());

  // copy child Groups to new Group
  IndexType gidx = group->getFirstValidGroupIndex();
  while(indexIsValid(gidx))
  {
    res->copyGroup(group->getGroup(gidx));
    gidx = group->getNextValidGroupIndex(gidx);
  }

  // copy Views to new Group
  IndexType vidx = group->getFirstValidViewIndex();
  while(indexIsValid(vidx))
  {
    res->copyView(group->getView(vidx));
    vidx = group->getNextValidViewIndex(vidx);
  }

  return res;
}

/*
 *************************************************************************
 *
 * Copy Group native layout to given Conduit node.
 *
 *************************************************************************
 */
bool Group::createNativeLayout(Node& n, const Attribute* attr) const
{
  n.set(DataType::object());
  bool hasSavedViews = false;

  // Dump the group's views
  IndexType vidx = getFirstValidViewIndex();
  while(indexIsValid(vidx))
  {
    const View* view = getView(vidx);

    // Check that the view's name is not also a child group name
    SLIC_CHECK_MSG(m_is_list || !hasChildGroup(view->getName()),
                   SIDRE_GROUP_LOG_PREPEND
                     << "'" << view->getName()
                     << "' is the name of both a group and a view.");

    if(attr == nullptr || view->hasAttributeValue(attr))
    {
      conduit::Node& child_node = m_is_list ? n.append() : n[view->getName()];
      view->createNativeLayout(child_node);
      hasSavedViews = true;
    }
    vidx = getNextValidViewIndex(vidx);
  }

  // Recursively dump the child groups
  IndexType gidx = getFirstValidGroupIndex();
  while(indexIsValid(gidx))
  {
    const Group* group = getGroup(gidx);
    conduit::Node& child_node = m_is_list ? n.append() : n[group->getName()];
    if(group->createNativeLayout(child_node, attr))
    {
      hasSavedViews = true;
    }
    else
    {
      if(m_is_list)
      {
        n.remove(group->getName());
      }
      else
      {
        n.remove(n.number_of_children() - 1);
      }
    }
    gidx = getNextValidGroupIndex(gidx);
  }

  return hasSavedViews;
}

/*
 *************************************************************************
 *
 * Adds a conduit node for this group if it has external views,
 * or if any of its children groups has an external view
 *
 *************************************************************************
 * see ATK-736 - Improvements to createNativeLayout and createExternalLayout
 */
bool Group::createExternalLayout(Node& n, const Attribute* attr) const
{
  n.set(DataType::object());

  bool hasExternalViews = false;

  // Dump the group's views
  IndexType vidx = getFirstValidViewIndex();
  while(indexIsValid(vidx))
  {
    const View* view = getView(vidx);

    // Check that the view's name is not also a child group name
    SLIC_CHECK_MSG(m_is_list || !hasChildGroup(view->getName()),
                   SIDRE_GROUP_LOG_PREPEND
                     << "'" << view->getName()
                     << "' is the name of both a group and a view.");

    if(attr == nullptr || view->hasAttributeValue(attr))
    {
      if(view->isExternal())
      {
        if(view->isDescribed())
        {
          conduit::Node& vnode = m_is_list ? n.append() : n[view->getName()];
          view->createNativeLayout(vnode);
        }
        hasExternalViews = true;
      }
    }

    vidx = getNextValidViewIndex(vidx);
  }

  // Recursively dump the child groups
  IndexType gidx = getFirstValidGroupIndex();
  while(indexIsValid(gidx))
  {
    const Group* group = getGroup(gidx);

    conduit::Node& gnode = m_is_list ? n.append() : n[group->getName()];
    if(group->createExternalLayout(gnode, attr))
    {
      hasExternalViews = true;
    }
    else
    {
      // Remove nodes that do not have any external views
      if(m_is_list)
      {
        n.remove(n.number_of_children() - 1);
      }
      else
      {
        n.remove(group->getName());
      }
    }

    gidx = getNextValidGroupIndex(gidx);
  }

  return hasExternalViews;
}

/*
 *************************************************************************
 *
 * Print JSON description of data Group to stdout.
 *
 *************************************************************************
 */
void Group::print() const { print(std::cout); }

/*
 *************************************************************************
 *
 * Print JSON description of data Group to an ostream.
 *
 *************************************************************************
 */
void Group::print(std::ostream& os) const
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
void Group::printTree(const int nlevels, std::ostream& os) const
{
  for(int i = 0; i < nlevels; ++i)
  {
    os << "    ";
  }
  os << "Group " << this->getName() << std::endl;

  IndexType vidx = getFirstValidViewIndex();
  while(indexIsValid(vidx))
  {
    const View* view = getView(vidx);

    for(int i = 0; i < nlevels + 1; ++i)
    {
      os << "    ";
    }
    os << "View " << view->getName() << std::endl;

    vidx = getNextValidViewIndex(vidx);
  }

  IndexType gidx = getFirstValidGroupIndex();
  while(indexIsValid(gidx))
  {
    const Group* group = getGroup(gidx);

    group->printTree(nlevels + 1, os);

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
void Group::copyToConduitNode(Node& n) const
{
  n["name"] = m_name;

  IndexType vidx = getFirstValidViewIndex();
  while(indexIsValid(vidx))
  {
    const View* view = getView(vidx);
    if(isUsingMap())
    {
      Node& v = n["views"].fetch(view->getName());
      view->copyToConduitNode(v);
    }
    else
    {
      Node& v = n["views"].append();
      view->copyToConduitNode(v);
    }
    vidx = getNextValidViewIndex(vidx);
  }

  IndexType gidx = getFirstValidGroupIndex();
  while(indexIsValid(gidx))
  {
    const Group* group = getGroup(gidx);
    if(isUsingMap())
    {
      Node& g = n["groups"].fetch(group->getName());
      group->copyToConduitNode(g);
    }
    else
    {
      Node& g = n["groups"].append();
      group->copyToConduitNode(g);
    }
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
bool Group::isEquivalentTo(const Group* other, bool checkName) const
{
  // Equality of names
  bool is_equiv = true;
  if(checkName)
  {
    is_equiv = (m_name == other->m_name);
  }

  // Sizes of collections of child items must be equal
  if(is_equiv)
  {
    is_equiv = (m_view_coll->getNumItems() == other->m_view_coll->getNumItems()) &&
      (m_group_coll->getNumItems() == other->m_group_coll->getNumItems());
  }

  // Test equivalence of Views
  if(is_equiv)
  {
    IndexType vidx = getFirstValidViewIndex();
    while(is_equiv && indexIsValid(vidx))
    {
      const View* view = getView(vidx);
      const std::string& name = view->getName();

      is_equiv =
        other->hasChildView(name) && view->isEquivalentTo(other->getView(name));

      vidx = getNextValidViewIndex(vidx);
    }
  }

  // Recursively call this method to test equivalence of child Groups
  if(is_equiv)
  {
    IndexType gidx = getFirstValidGroupIndex();
    while(is_equiv && indexIsValid(gidx))
    {
      const Group* group = getGroup(gidx);
      const std::string& name = group->getName();

      is_equiv = other->hasChildGroup(name) &&
        group->isEquivalentTo(other->getGroup(name));

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

void Group::save(const std::string& path,
                 const std::string& protocol,
                 const Attribute* attr) const
{
  const DataStore* ds = getDataStore();

  if(protocol == "sidre_hdf5")
  {
    Node n;
    exportTo(n["sidre"], attr);
    ds->saveAttributeLayout(n["sidre/attribute"]);
    createExternalLayout(n["sidre/external"], attr);
    n["sidre_group_name"] = m_name;
    conduit::relay::io::save(n, path, "hdf5");
  }
  else if(protocol == "sidre_conduit_json")
  {
    Node n;
    exportTo(n["sidre"], attr);
    ds->saveAttributeLayout(n["sidre/attribute"]);
    createExternalLayout(n["sidre/external"], attr);
    n["sidre_group_name"] = m_name;
    conduit::relay::io::save(n, path, "conduit_json");
  }
  else if(protocol == "sidre_json")
  {
    Node n;
    exportTo(n["sidre"], attr);
    ds->saveAttributeLayout(n["sidre/attribute"]);
    createExternalLayout(n["sidre/external"], attr);
    n["sidre_group_name"] = m_name;
    conduit::relay::io::save(n, path, "json");
  }
  else if(protocol == "conduit_hdf5")
  {
    Node n;
    createNativeLayout(n, attr);
    n["sidre_group_name"] = m_name;
    conduit::relay::io::save(n, path, "hdf5");
  }
  else if(protocol == "conduit_bin" || protocol == "conduit_json" ||
          protocol == "json")
  {
    Node n;
    createNativeLayout(n, attr);
    n["sidre_group_name"] = m_name;
    conduit::relay::io::save(n, path, protocol);
  }
  else
  {
    SLIC_ERROR(SIDRE_GROUP_LOG_PREPEND << "Invalid protocol '" << protocol
                                       << "' for file save.");
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
void Group::load(const std::string& path,
                 const std::string& protocol,
                 bool preserve_contents)
{
  std::string new_name;
  load(path, protocol, preserve_contents, new_name);
}

void Group::load(const std::string& path,
                 const std::string& protocol,
                 bool preserve_contents,
                 std::string& name_from_file)
{
  if(protocol == "sidre_hdf5")
  {
    Node n;
    conduit::relay::io::load(path, "hdf5", n);
    SLIC_ASSERT_MSG(n.has_path("sidre"),
                    SIDRE_GROUP_LOG_PREPEND
                      << "Conduit Node " << n.path() << " does not have sidre "
                      << "data for this Group " << getPathName() << ".");
    importFrom(n["sidre"], preserve_contents);
    if(n.has_path("sidre_group_name"))
    {
      name_from_file = n["sidre_group_name"].as_string();
    }
  }
  else if(protocol == "sidre_conduit_json")
  {
    Node n;
    conduit::relay::io::load(path, "conduit_json", n);
    SLIC_ASSERT_MSG(n.has_path("sidre"),
                    SIDRE_GROUP_LOG_PREPEND
                      << "Conduit Node " << n.path() << " does not have sidre "
                      << "data for Group " << getPathName() << ".");
    importFrom(n["sidre"], preserve_contents);
    if(n.has_path("sidre_group_name"))
    {
      name_from_file = n["sidre_group_name"].as_string();
    }
  }
  else if(protocol == "sidre_json")
  {
    Node n;
    conduit::relay::io::load(path, "json", n);
    SLIC_ASSERT_MSG(n.has_path("sidre"),
                    SIDRE_GROUP_LOG_PREPEND
                      << "Conduit Node " << n.path() << " does not have sidre "
                      << "data for Group " << getPathName() << ".");
    importFrom(n["sidre"], preserve_contents);
    if(n.has_path("sidre_group_name"))
    {
      name_from_file = n["sidre_group_name"].as_string();
    }
  }
  else if(protocol == "conduit_hdf5")
  {
    Node n;
    conduit::relay::io::load(path, "hdf5", n);
    importConduitTree(n, preserve_contents);
    if(n.has_path("sidre_group_name"))
    {
      name_from_file = n["sidre_group_name"].as_string();
    }
  }
  else if(protocol == "conduit_bin" || protocol == "conduit_json" ||
          protocol == "json")
  {
    Node n;
    conduit::relay::io::load(path, protocol, n);
    importConduitTree(n, preserve_contents);
    if(n.has_path("sidre_group_name"))
    {
      name_from_file = n["sidre_group_name"].as_string();
    }
  }
  else
  {
    SLIC_ERROR(SIDRE_GROUP_LOG_PREPEND << "Invalid protocol '" << protocol
                                       << "' for file load.");
  }
}

/*
 *************************************************************************
 *
 * Load Group (including Views and child Groups) from a file
 *
 *************************************************************************
 */
Group* Group::createGroupAndLoad(std::string& group_name,
                                 const std::string& path,
                                 const std::string& protocol,
                                 bool& load_success)
{
  load_success = false;
  Group* child = createGroup(group_name);
  if(child != nullptr)
  {
    // In a forthcoming PR, load() will return a bool for success/failure
    load_success = true;
    child->load(path, protocol, false, group_name);
  }

  return child;
}

/*
 *************************************************************************
 *
 * Load External Data from a file
 *
 *************************************************************************
 */
void Group::loadExternalData(const std::string& path)
{
  Node n;
  createExternalLayout(n);

#ifdef AXOM_USE_HDF5
  // CYRUS'-NOTE, not sure ":" will work with multiple trees per
  // output file
  conduit::relay::io::hdf5_read(path + ":sidre/external", n);
#else
  AXOM_UNUSED_VAR(path);
  SLIC_WARNING(SIDRE_GROUP_LOG_PREPEND
               << "External data not loaded. "
               << "This function requires hdf5 support. "
               << " Please reconfigure with hdf5.");
#endif
}

// Functions that directly use the hdf5 API in their signature
#ifdef AXOM_USE_HDF5

/*
 *************************************************************************
 *
 * Save Group (including Views and child Groups) to a hdf5 handle
 *
 *************************************************************************
 */
void Group::save(const hid_t& h5_id,
                 const std::string& protocol,
                 const Attribute* attr) const
{
  // supported here:
  // "sidre_hdf5"
  // "conduit_hdf5"
  if(protocol == "sidre_hdf5")
  {
    Node n;
    exportTo(n["sidre"], attr);
    createExternalLayout(n["sidre/external"], attr);
    n["sidre_group_name"] = m_name;
    conduit::relay::io::hdf5_write(n, h5_id);
  }
  else if(protocol == "conduit_hdf5")
  {
    Node n;
    createNativeLayout(n, attr);
    n["sidre_group_name"] = m_name;
    conduit::relay::io::hdf5_write(n, h5_id);
  }
  else
  {
    SLIC_ERROR(SIDRE_GROUP_LOG_PREPEND << "Invalid protocol '" << protocol
                                       << "' for save with hdf5 handle.");
  }
}

/*
 *************************************************************************
 *
 * Load Group (including Views and child Groups) from an hdf5 handle
 *
 *************************************************************************
 */
void Group::load(const hid_t& h5_id,
                 const std::string& protocol,
                 bool preserve_contents)
{
  std::string name_from_file;
  load(h5_id, protocol, preserve_contents, name_from_file);
}

void Group::load(const hid_t& h5_id,
                 const std::string& protocol,
                 bool preserve_contents,
                 std::string& name_from_file)
{
  // supported here:
  // "sidre_hdf5"
  // "conduit_hdf5"
  if(protocol == "sidre_hdf5")
  {
    Node n;
    conduit::relay::io::hdf5_read(h5_id, n);
    SLIC_ASSERT_MSG(n.has_path("sidre"),
                    SIDRE_GROUP_LOG_PREPEND
                      << "Conduit Node " << n.path() << " does not have sidre "
                      << "data for Group " << getPathName() << ".");
    importFrom(n["sidre"], preserve_contents);
    if(n.has_path("sidre_group_name"))
    {
      name_from_file = n["sidre_group_name"].as_string();
    }
  }
  else if(protocol == "conduit_hdf5")
  {
    SLIC_ERROR("Protocol " << protocol << " not yet supported for file load.");
    Node n;
    conduit::relay::io::hdf5_read(h5_id, n);
    importConduitTree(n, preserve_contents);
    if(n.has_path("sidre_group_name"))
    {
      name_from_file = n["sidre_group_name"].as_string();
    }
  }
  else
  {
    SLIC_ERROR(SIDRE_GROUP_LOG_PREPEND << "Invalid protocol '" << protocol
                                       << "' for file load.");
  }
}

/*
 *************************************************************************
 *
 * Load External Data from an hdf5 file
 *
 * Note: this ASSUMES uses the "sidre_hdf5" protocol
 *************************************************************************
 */
void Group::loadExternalData(const hid_t& h5_id)
{
  Node n;
  createExternalLayout(n);
  conduit::relay::io::hdf5_read(h5_id, "sidre/external", n);
}

#endif /* AXOM_USE_HDF5 */

////////////////////////////////////////////////////////////////////////
//
// Private methods below
//
////////////////////////////////////////////////////////////////////////

/*
 *************************************************************************
 *
 * PRIVATE ctor makes Group with given name and make it a child of
 * root Group in datastore.
 *
 *************************************************************************
 */
Group::Group(const std::string& name, DataStore* datastore, bool is_list)
  : m_name(name)
  , m_index(InvalidIndex)
  , m_parent(nullptr)
  , m_datastore(datastore)
  , m_is_list(is_list)
  , m_view_coll(nullptr)
  , m_group_coll(nullptr)
#ifdef AXOM_USE_UMPIRE
  , m_default_allocator_id(axom::getDefaultAllocatorID())
#endif
{
  if(is_list)
  {
    m_view_coll = new ListCollection<View>();
    m_group_coll = new ListCollection<Group>();
  }
  else
  {
    m_view_coll = new MapCollection<View>();
    m_group_coll = new MapCollection<Group>();
  }
}

/*
 *************************************************************************
 *
 * PRIVATE dtor destroys Group and all its contents.
 *
 *************************************************************************
 */
Group::~Group()
{
  destroyViews();
  destroyGroups();
  delete m_view_coll;
  delete m_group_coll;
}

/*
 *************************************************************************
 *
 * PRIVATE method to attach given View to Group.
 *
 *************************************************************************
 */
View* Group::attachView(View* view)
{
  if(view == nullptr ||
     (!view->getName().empty() && hasChildView(view->getName())))
  {
    return nullptr;
  }
  else
  {
    SLIC_ASSERT_MSG(view->m_owning_group == nullptr,
                    SIDRE_GROUP_LOG_PREPEND
                      << "Provided View " << view->getPathName() << " is already "
                      << "attatched to Group "
                      << view->m_owning_group->getPathName() << " and can't be "
                      << "attatched to Group " << getPathName() << ".");
    view->m_owning_group = this;
    view->m_index = m_view_coll->insertItem(view, view->getName());
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
View* Group::detachView(const std::string& name)
{
  View* view = getNamedViews()->removeItem(name);
  if(view != nullptr)
  {
    view->m_owning_group = nullptr;
    view->m_index = InvalidIndex;
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
View* Group::detachView(IndexType idx)
{
  View* view = m_view_coll->removeItem(idx);
  if(view != nullptr)
  {
    view->m_owning_group = nullptr;
    view->m_index = InvalidIndex;
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
void Group::destroyViewAndData(View* view)
{
  if(view != nullptr)
  {
    Group* group = view->getOwningGroup();
    group->detachView(view->getName());
    Buffer* const buffer = view->detachBuffer();
    if(buffer != nullptr && buffer->getNumViews() == 0)
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
Group* Group::attachGroup(Group* group)
{
  if(group == nullptr ||
     (!group->getName().empty() && hasChildGroup(group->getName())))
  {
    return nullptr;
  }
  else
  {
    group->m_parent = this;
    group->m_index = m_group_coll->insertItem(group, group->getName());
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
Group* Group::detachGroup(const std::string& name)
{
  Group* group = getNamedGroups()->removeItem(name);
  if(group != nullptr)
  {
    group->m_parent = nullptr;
    group->m_index = InvalidIndex;
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
Group* Group::detachGroup(IndexType idx)
{
  Group* group = m_group_coll->removeItem(idx);
  if(group != nullptr)
  {
    group->m_parent = nullptr;
    group->m_index = InvalidIndex;
  }

  return group;
}

/*
 *************************************************************************
 *
 * Serialize tree identified by a Group into a conduit node.  Include
 * any Buffers attached to Views in that tree.
 *
 * Note: This is for the "sidre_{zzz}" protocols.
 *
 *************************************************************************
 */
bool Group::exportTo(conduit::Node& result, const Attribute* attr) const
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
  bool hasSavedViews = exportTo(result, attr, buffer_indices);

  if(!buffer_indices.empty())
  {
    // Now, add all the referenced buffers to the node.
    Node& bnode = result["buffers"];
    for(std::set<IndexType>::iterator s_it = buffer_indices.begin();
        s_it != buffer_indices.end();
        ++s_it)
    {
      // Use a dictionary layout here instead of conduit list.
      // Conduit IO HDF5 doesn't support conduit list objects.
      std::ostringstream oss;
      oss << "buffer_id_" << *s_it;
      Node& n_buffer = bnode.fetch(oss.str());
      getDataStore()->getBuffer(*s_it)->exportTo(n_buffer);
    }
  }

  return hasSavedViews;
}

/*
 *************************************************************************
 *
 * PRIVATE method to copy from Group to given Conduit node using
 * given set of ids to maintain correct association of data Buffers
 * to data Views.
 *
 * Note: This is for the "sidre_{zzz}" protocols.
 *
 *************************************************************************
 */
bool Group::exportTo(conduit::Node& result,
                     const Attribute* attr,
                     std::set<IndexType>& buffer_indices) const
{
  result.set(DataType::object());
  bool hasSavedViews = false;

  if(getNumViews() > 0)
  {
    Node& vnode = result["views"];
    IndexType vidx = getFirstValidViewIndex();
    while(indexIsValid(vidx))
    {
      const View* view = getView(vidx);
      if(attr == nullptr || view->hasAttributeValue(attr))
      {
        Node& n_view = m_is_list ? vnode.append() : vnode.fetch(view->getName());
        view->exportTo(n_view, buffer_indices);
        hasSavedViews = true;
      }
      vidx = getNextValidViewIndex(vidx);
    }
    if(!hasSavedViews)
    {
      result.remove("views");
    }
  }

  bool hasSavedGroups = false;
  if(getNumGroups() > 0)
  {
    Node& gnode = result["groups"];
    IndexType gidx = getFirstValidGroupIndex();
    while(indexIsValid(gidx))
    {
      const Group* group = getGroup(gidx);
      Node& n_group = m_is_list ? gnode.append() : gnode.fetch(group->getName());
      bool hsv = group->exportTo(n_group, attr, buffer_indices);
      hasSavedViews = hasSavedViews || hsv;
      hasSavedGroups = true;

      gidx = getNextValidGroupIndex(gidx);
    }
    if(!hasSavedGroups)
    {
      result.remove("groups");
    }
  }

  return hasSavedViews;
}

/*
 *************************************************************************
 *
 * Imports tree from a conduit node into this Group.  Includes
 * any Buffers attached to Views in that tree.
 *
 *
 * Note: This is for the "sidre_{zzz}" protocols.
 *
 *************************************************************************
 */

void Group::importFrom(conduit::Node& node, bool preserve_contents)
{
  // TODO - May want to put in a little meta-data into these files like a
  // 'version'
  // or tag identifying the data.  We don't want someone giving us a file that
  // doesn't have our full multiView->buffer connectivity in there.

  if(!preserve_contents)
  {
    destroyGroups();
    destroyViews();
  }

  getDataStore()->loadAttributeLayout(node);

  // First - Import Buffers into the DataStore.
  std::map<IndexType, IndexType> buffer_indices_map;

  if(node.has_path("buffers"))
  {
    conduit::NodeIterator buffs_itr = node["buffers"].children();
    while(buffs_itr.has_next())
    {
      Node& n_buffer = buffs_itr.next();
      IndexType old_buffer_id = n_buffer["id"].to_int64();

      Buffer* buffer = getDataStore()->createBuffer();

      // track change of old Buffer id to new Buffer id
      buffer_indices_map[old_buffer_id] = buffer->getIndex();

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
 * Note: This is for the "sidre_{zzz}" protocols.
 *
 *************************************************************************
 */
void Group::importFrom(conduit::Node& node,
                       const std::map<IndexType, IndexType>& buffer_id_map)
{
  if(node.has_path("views"))
  {
    // create the Views
    conduit::NodeIterator views_itr = node["views"].children();
    while(views_itr.has_next())
    {
      Node& n_view = views_itr.next();
      std::string view_name = m_is_list ? "" : views_itr.name();

      View* view = createView(view_name);
      view->importFrom(n_view, buffer_id_map);
    }
  }
  if(node.has_path("groups"))
  {
    // create the child Groups
    conduit::NodeIterator groups_itr = node["groups"].children();
    while(groups_itr.has_next())
    {
      Group* group;
      Node& n_group = groups_itr.next();
      bool create_list = false;
      if(n_group.has_child("views"))
      {
        if(n_group["views"].dtype().is_list())
        {
          create_list = true;
        }
      }
      if(!create_list && n_group.has_child("groups"))
      {
        if(n_group["groups"].dtype().is_list())
        {
          create_list = true;
        }
      }
      std::string group_name = groups_itr.name();
      if(m_is_list)
      {
        group = createUnnamedGroup(create_list);
      }
      else
      {
        std::string group_name = groups_itr.name();
        group = createGroup(group_name, create_list);
      }
      group->importFrom(n_group, buffer_id_map);
    }
  }
}

/*
 *************************************************************************
 *
 * Imports tree from a conduit node into this Group.
 * This takes a generic conduit tree, not one with sidre conventions.
 *
 *************************************************************************
 */

bool Group::importConduitTree(const conduit::Node& node, bool preserve_contents)
{
  bool success = true;
  if(!preserve_contents)
  {
    destroyGroups();
    destroyViews();
  }

  //
  DataType node_dtype = node.dtype();
  if(node_dtype.is_object() || node_dtype.is_list())
  {
    conduit::NodeConstIterator itr = node.children();
    while(itr.has_next())
    {
      const Node& cld_node = itr.next();
      std::string cld_name = m_is_list ? "" : itr.name();
      DataType cld_dtype = cld_node.dtype();

      if(cld_dtype.is_object() || cld_dtype.is_list())
      {
        // create group
        Group* grp = m_is_list ? createUnnamedGroup(cld_dtype.is_list())
                               : createGroup(cld_name, cld_dtype.is_list());
        success = grp->importConduitTree(cld_node, preserve_contents);
      }
      else if(cld_dtype.is_empty())
      {
        //create empty view
        createView(cld_name);
      }
      else if(cld_dtype.is_string())
      {
        if(cld_name != "sidre_group_name")
        {
          //create string view
          createViewString(cld_name, cld_node.as_string());
        }
      }
      else if(cld_dtype.is_number())
      {
        if(cld_dtype.number_of_elements() == 1)
        {
          // create scalar view
          View* view = createView(cld_name);
          view->setScalar(cld_node);
        }
        else
        {
          View* view = createView(cld_name);
          view->importArrayNode(cld_node);
        }
      }
      else
      {
        // All Nodes should have one of the above datatypes, so if
        // we get here something is wrong.
        SLIC_ERROR(SIDRE_GROUP_LOG_PREPEND
                   << "Conduit child Node " << cld_name
                   << " does not have a recognized datatype."
                   << " Cannot import into Group " << getPathName());
      }
    }
  }
  else
  {
    SLIC_ERROR(SIDRE_GROUP_LOG_PREPEND
               << "Group cannot import non-object Conduit Node");
  }

  return success;
}

bool Group::importConduitTreeExternal(conduit::Node& node, bool preserve_contents)
{
  bool success = true;
  if(!preserve_contents)
  {
    destroyGroups();
    destroyViews();
  }

  DataType node_dtype = node.dtype();
  if(node_dtype.is_object() || node_dtype.is_list())
  {
    conduit::NodeIterator itr = node.children();
    while(itr.has_next())
    {
      Node& cld_node = itr.next();
      std::string cld_name = m_is_list ? "" : itr.name();
      DataType cld_dtype = cld_node.dtype();

      if(cld_dtype.is_object() || cld_dtype.is_list())
      {
        // create group
        Group* grp = m_is_list ? createUnnamedGroup(cld_dtype.is_list())
                               : createGroup(cld_name, cld_dtype.is_list());
        success = grp->importConduitTreeExternal(cld_node, preserve_contents);
      }
      else if(cld_dtype.is_empty())
      {
        //create empty view
        createView(cld_name);
      }
      else if(cld_dtype.is_string())
      {
        if(cld_name != "sidre_group_name")
        {
          //create string view
          createViewString(cld_name, cld_node.as_string());
        }
      }
      else if(cld_dtype.is_number())
      {
        if(cld_dtype.number_of_elements() == 1)
        {
          // create scalar view
          View* view = createView(cld_name);
          view->setScalar(cld_node);
        }
        else
        {
          void* conduit_ptr = cld_node.data_ptr();
          View* view = createView(cld_name);
          view->setExternalDataPtr(conduit_ptr);
          view->apply(cld_dtype);
        }
      }
      else
      {
        // All Nodes should have one of the above datatypes, so if
        // we get here something is wrong.
        SLIC_ERROR(SIDRE_GROUP_LOG_PREPEND
                   << "Conduit child Node " << cld_name
                   << " does not have a recognized datatype."
                   << " Cannot import into Group " << getPathName());
      }
    }
  }
  else
  {
    SLIC_ERROR(SIDRE_GROUP_LOG_PREPEND
               << "Group cannot import non-object Conduit Node");
  }

  return success;
}

/*
 *************************************************************************
 *
 * PRIVATE method to walk down a path to the next-to-last entry.
 *
 * If an error is encountered, this private function will return nullptr
 *
 *************************************************************************
 */
Group* Group::walkPath(std::string& path, bool create_groups_in_path)
{
  Group* group_ptr = this;

  // Split path into parts
  std::vector<std::string> path_parts =
    axom::Path(path, s_path_delimiter).parts();

  if(path_parts.size() > 0)
  {
    if(m_is_list && path_parts.size() > 1)
    {
      // A size > 1 indicates that a delimited path string was provided.
      // A path string is invalid when this Group uses the list format.
      SLIC_WARNING(SIDRE_GROUP_LOG_PREPEND
                   << "A delimited path string '" << path
                   << "' cannot be used as the name of an object "
                   << "to be created by a Group that uses the list format. "
                   << "A null pointer will be returned.");
      group_ptr = nullptr;
    }
    else
    {
      // Find stopping point (right before last part of path)
      std::vector<std::string>::const_iterator stop = path_parts.end() - 1;

      // Navigate path down to desired Group
      for(std::vector<std::string>::const_iterator iter = path_parts.begin();
          iter < stop;
          ++iter)
      {
        if(group_ptr->hasChildGroup(*iter))
        {
          group_ptr = group_ptr->getGroup(*iter);
        }
        else if(create_groups_in_path)
        {
          group_ptr = group_ptr->createGroup(*iter);

          if(group_ptr == nullptr)
          {
            iter = stop;
          }
        }
        else
        {
          iter = stop;
          group_ptr = nullptr;
        }
      }
      path = path_parts.back();
    }
  }

  return group_ptr;
}

/*
 *************************************************************************
 *
 * PRIVATE const method to walk down a path to the next-to-last entry.
 *
 * If an error is encountered, this private function will return nullptr
 *
 *************************************************************************
 */
const Group* Group::walkPath(std::string& path) const
{
  const Group* group_ptr = this;

  // Split path into parts
  std::vector<std::string> path_parts =
    axom::Path(path, s_path_delimiter).parts();

  if(path_parts.size() > 0)
  {
    // Find stopping point (right before last part of path)
    std::vector<std::string>::const_iterator stop = path_parts.end() - 1;

    // Navigate path down to desired Group
    for(std::vector<std::string>::const_iterator iter = path_parts.begin();
        iter < stop;
        ++iter)
    {
      if(group_ptr->hasChildGroup(*iter))
      {
        group_ptr = group_ptr->getGroup(*iter);
      }
      else
      {
        group_ptr = nullptr;
        iter = stop;
      }
    }
    path = path_parts.back();
  }

  return group_ptr;
}

/*
 *************************************************************************
 *
 * Return number of child Groups in a Group object.
 *
 *************************************************************************
 */
IndexType Group::getNumGroups() const { return m_group_coll->getNumItems(); }

/*
 *************************************************************************
 *
 * Return number of Views owned by a Group object.
 *
 *************************************************************************
 */
IndexType Group::getNumViews() const { return m_view_coll->getNumItems(); }

/*
 *************************************************************************
 *
 * Return true if this Group owns a View with given name (not path);
 * else false.
 *
 *************************************************************************
 */
bool Group::hasChildView(const std::string& name) const
{
  return isUsingMap() ? getNamedViews()->hasItem(name) : false;
}

/*
 *************************************************************************
 *
 * Return true if this Group owns a View with given index; else false.
 *
 *************************************************************************
 */
bool Group::hasView(IndexType idx) const { return m_view_coll->hasItem(idx); }

/*
 *************************************************************************
 *
 * Return index of View with given name owned by this Group object.
 *
 * If no such View exists, return sidre::InvalidIndex;
 *
 *************************************************************************
 */
IndexType Group::getViewIndex(const std::string& name) const
{
  SLIC_CHECK_MSG(
    hasChildView(name),
    SIDRE_GROUP_LOG_PREPEND << "Group has no View with name '" << name << "'");

  return getNamedViews()->getItemIndex(name);
}

/*
 *************************************************************************
 *
 * Return name of View with given index owned by Group object.
 *
 * If no such View exists, return sidre::InvalidName.
 *
 *************************************************************************
 */
const std::string& Group::getViewName(IndexType idx) const
{
  SLIC_CHECK_MSG(
    hasView(idx),
    SIDRE_GROUP_LOG_PREPEND << "Group has no View with index " << idx);

  return getNamedViews()->getItemName(idx);
}

/*
 *************************************************************************
 *
 * Return pointer to non-const View with given index.
 *
 * If no such View exists, nullptr is returned.
 *
 *************************************************************************
 */
View* Group::getView(IndexType idx)
{
  SLIC_CHECK_MSG(
    hasView(idx),
    SIDRE_GROUP_LOG_PREPEND << "Group has no View with index " << idx);

  return m_view_coll->getItem(idx);
}

/*
 *************************************************************************
 *
 * Return pointer to const View with given index.
 *
 * If no such View exists, nullptr is returned.
 *
 *************************************************************************
 */
const View* Group::getView(IndexType idx) const
{
  SLIC_CHECK_MSG(
    hasView(idx),
    SIDRE_GROUP_LOG_PREPEND << "Group has no View with index " << idx);

  return m_view_coll->getItem(idx);
}

/*
 *************************************************************************
 *
 * Return first valid View index in Group object
 *        (i.e., smallest index over all Views).
 *
 * sidre::InvalidIndex is returned if Group has no Views.
 *
 *************************************************************************
 */
IndexType Group::getFirstValidViewIndex() const
{
  return m_view_coll->getFirstValidIndex();
}

/*
 *************************************************************************
 *
 * Return next valid View index in Group object after given index
 *        (i.e., smallest index over all View indices larger than given one).
 *
 * sidre::InvalidIndex is returned if there is no valid index greater
 * than given one.
 *
 *************************************************************************
 */
IndexType Group::getNextValidViewIndex(IndexType idx) const
{
  return m_view_coll->getNextValidIndex(idx);
}

/*!
 * \brief Returns an adaptor to support iterating the collection of views
 */
typename Group::ViewCollection::iterator_adaptor Group::views()
{
  return m_view_coll->getIteratorAdaptor();
}

/*!
 * \brief Returns a const adaptor to support iterating the collection of views
 */
typename Group::ViewCollection::const_iterator_adaptor Group::views() const
{
  return m_view_coll->getIteratorAdaptor();
}

/*!
 * \brief Returns an adaptor to support iterating the collection of groups
 */
typename Group::GroupCollection::iterator_adaptor Group::groups()
{
  return m_group_coll->getIteratorAdaptor();
}

/*!
 * \brief Returns a const adaptor to support iterating the collection of groups
 */
typename Group::GroupCollection::const_iterator_adaptor Group::groups() const
{
  return m_group_coll->getIteratorAdaptor();
}

/*
 *************************************************************************
 *
 * Return true if this Group has a child Group with given
 * name; else false.
 *
 *************************************************************************
 */
bool Group::hasChildGroup(const std::string& name) const
{
  return isUsingMap() ? getNamedGroups()->hasItem(name) : false;
}

/*
 *************************************************************************
 *
 * Return true if Group has an immediate child Group
 * with given index; else false.
 *
 *************************************************************************
 */
bool Group::hasGroup(IndexType idx) const { return m_group_coll->hasItem(idx); }

/*
 *************************************************************************
 * Return the index of immediate child Group with given name.
 *
 * If no such child Group exists, return sidre::InvalidIndex;
 *
 *************************************************************************
 */
IndexType Group::getGroupIndex(const std::string& name) const
{
  SLIC_CHECK_MSG(hasChildGroup(name),
                 SIDRE_GROUP_LOG_PREPEND
                   << "Group has no child Group with name '" << name << "'");

  return getNamedGroups()->getItemIndex(name);
}

/*
 *************************************************************************
 *
 * Return the name of immediate child Group with given index.
 *
 * If no such child Group exists, return sidre::InvalidName.
 *
 *************************************************************************
 */
const std::string& Group::getGroupName(IndexType idx) const
{
  SLIC_CHECK_MSG(
    hasGroup(idx),
    SIDRE_GROUP_LOG_PREPEND << "Group has no child Group with index " << idx);

  return getNamedGroups()->getItemName(idx);
}

/*
 *************************************************************************
 *
 * Return pointer to non-const immediate child Group with given index.
 *
 * If no such Group exists, nullptr is returned.
 *
 *************************************************************************
 */
Group* Group::getGroup(IndexType idx)
{
  SLIC_CHECK_MSG(
    hasGroup(idx),
    SIDRE_GROUP_LOG_PREPEND << "Group has no child Group with index " << idx);

  return m_group_coll->getItem(idx);
}

/*
 *************************************************************************
 *
 * Return pointer to const immediate child Group with given index.
 *
 * If no such Group exists, nullptr is returned.
 *
 *************************************************************************
 */
const Group* Group::getGroup(IndexType idx) const
{
  SLIC_CHECK_MSG(
    hasGroup(idx),
    SIDRE_GROUP_LOG_PREPEND << "Group has no child Group with index " << idx);

  return m_group_coll->getItem(idx);
}

/*
 *************************************************************************
 *
 * Return first valid child Group index (i.e., smallest
 *        index over all child Groups).
 *
 * sidre::InvalidIndex is returned if Group has no child Groups.
 *
 *************************************************************************
 */
IndexType Group::getFirstValidGroupIndex() const
{
  return m_group_coll->getFirstValidIndex();
}

/*
 *************************************************************************
 *
 * Return next valid child Group index after given index
 *        (i.e., smallest index over all child Group indices larger
 *        than given one).
 *
 * sidre::InvalidIndex is returned if there is no valid index greater
 * than given one.
 *
 *************************************************************************
 */
IndexType Group::getNextValidGroupIndex(IndexType idx) const
{
  return m_group_coll->getNextValidIndex(idx);
}

/*
 *************************************************************************
 *
 * Rename this Group with a new string name.
 *
 *************************************************************************
 */
bool Group::rename(const std::string& new_name)
{
  bool do_rename = true;
  if(new_name != m_name)
  {
    if(new_name.empty())
    {
      SLIC_WARNING(SIDRE_GROUP_LOG_PREPEND
                   << "Cannot rename Group to an empty string.");
      do_rename = false;
    }
    else if(new_name.find(s_path_delimiter) != std::string::npos)
    {
      SLIC_WARNING(SIDRE_GROUP_LOG_PREPEND
                   << "Cannot rename Group to path name '" << new_name << "'. "
                   << "Only strings without path delimiters can "
                   << "be passed into the rename method.");
      do_rename = false;
    }

    if(do_rename)
    {
      const Group* root = getDataStore()->getRoot();
      Group* parent = getParent();

      //If this is the root group, we don't need
      //to do anything to change the parent's handle to this group.
      if(this != root && parent != nullptr)
      {
        if(parent->hasGroup(new_name) || parent->hasView(new_name))
        {
          SLIC_WARNING(SIDRE_GROUP_LOG_PREPEND
                       << "Parent group " << parent->getPathName()
                       << " already has a child group named '" << new_name
                       << "'. Group " << getPathName()
                       << " will not be renamed.");
          do_rename = false;
        }
        else
        {
          Group* detached_group = parent->detachGroup(m_name);
          SLIC_CHECK_MSG(detached_group == this,
                         SIDRE_GROUP_LOG_PREPEND
                           << "Group detatched from parent '"
                           << detached_group->getPathName() << "' is not "
                           << "this Group '" << getPathName() << "'.");

          m_name = new_name;

          Group* attached_group = parent->attachGroup(detached_group);
          AXOM_UNUSED_VAR(attached_group);
          SLIC_CHECK_MSG(attached_group == this,
                         SIDRE_GROUP_LOG_PREPEND
                           << "Group attached to parent '"
                           << attached_group->getPathName() << "' is not "
                           << "this Group '" << getPathName() << "'.");
        }
      }
      else
      {
        m_name = new_name;
      }
    }
  }

  return do_rename;
}

/*
 *************************************************************************
 *
 * PRIVATE method to return a valid umpire::Allocator ID.
 *
 *************************************************************************
 */
int Group::getValidAllocatorID(int allocID)
{
#ifdef AXOM_USE_UMPIRE
  if(allocID == INVALID_ALLOCATOR_ID)
  {
    allocID = m_default_allocator_id;
  }
#endif

  return allocID;
}

} /* end namespace sidre */
} /* end namespace axom */

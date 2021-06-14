// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file Group.hpp
 *
 * \brief   Header file containing definition of Group class.
 *
 ******************************************************************************
 */

#ifndef SIDRE_GROUP_HPP_
#define SIDRE_GROUP_HPP_

// axom headers
#include "axom/config.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/Types.hpp"
#include "axom/slic.hpp"

// Standard C++ headers
#include <memory>
#include <map>
#include <string>
#include <vector>
#include <set>
#include <cstring>

// third party lib headers
#ifdef AXOM_USE_HDF5
  #include "hdf5.h"
#endif

// Sidre headers
#include "SidreTypes.hpp"
#include "View.hpp"

// Define the default protocol for sidre I/O
#ifdef AXOM_USE_HDF5
  #define SIDRE_DEFAULT_PROTOCOL "sidre_hdf5"
#else
  #define SIDRE_DEFAULT_PROTOCOL "sidre_conduit_json"
#endif

namespace axom
{
namespace sidre
{
class Buffer;
class Group;
class DataStore;
template <typename TYPE>
class ItemCollection;

/*!
 * \class Group
 *
 * \brief Group holds a collection of Views and (child) Groups.
 *
 * The Group class has the following properties:
 *
 *    - Groups can be organized into a (tree) hierarchy by creating
 *      child Groups from the root Group owned by a DataStore object.
 *    - A Group object can only be created by another Group; the
 *      Group ctor is not visible externally. A Group is owned
 *      by the Group that creates it (its parent) and becomes a
 *      (direct) child Group of the parent. Groups in the subtree
 *      rooted at an ancestor Group are that Group's descendants.
 *    - A Group object has a unique name (string) within its parent
 *      Group.
 *    - A Group object maintains a pointer to its parent Group.
 *    - A Group object can be moved or copied to another Group.
 *    - Group objects can create View objects within them. The
 *      Group that creates a View owns it.
 *    - A View object has a unique name (string) within the Group
 *      that owns it.
 *    - A View object can be moved or copied to another Group.
 *
 * Note that Views and child Groups within a Group can be accessed
 * by name or index.
 *
 * Note that certain methods for querying, creating, retrieving, and
 * deleting Groups and Views take a string with path syntax,
 * while others take the name of a direct child of the current Group.
 * Methods that require the name of a direct child are marked with
 * "Child", for example hasChildView() and hasChildGroup().  When a path
 * string is passed to a method that accepts path syntax, the last item in
 * the path indicates the item to be created, accessed, etc.  For example,
 *
 *      View* view = group->createView("foo/bar/baz");
 *
 * is equivalent to
 *
 *      View* view =
 *        group->createGroup("foo")->createGroup("bar")->createView("baz");
 *
 * In particular, intermediate Groups "foo" and "bar" will be created in
 * this case if they don't already exist.
 *
 * Methods that access Views or Groups by index work with the direct
 * children of the current Group because an index has no meaning outside
 * of the indexed group.  None of these methods is marked with "Child".
 *
 * A Group can optionally be created to hold items in a "list format". In
 * this format, any number of child Group or View items can be created with
 * empty strings for their names, and none of the methods that access
 * child items by name or path will return a valid pointer.
 *
 * \attention when Views or Groups are created, destroyed, copied, or moved,
 * indices of other Views and Groups in associated Group objects may
 * become invalid. This is analogous to iterator invalidation for STL
 * containers when the container contents change.
 *
 */
class Group
{
public:
  //
  // Friend declarations to constrain usage via controlled access to
  // private members.
  //
  friend class DataStore;
  friend class View;

  //@{
  //!  @name Basic query and accessor methods.

  /*!
   * \brief Return the path delimiter
   */
  char getPathDelimiter() const { return s_path_delimiter; }

  /*!
   * \brief Return index of Group object within parent Group.
   */
  IndexType getIndex() const { return m_index; }

  /*!
   * \brief Return const reference to name of Group object.
   *
   * \sa getPath(), getPathName()
   */
  const std::string& getName() const { return m_name; }

  /*!
   * \brief Return path of Group object, not including its name.
   *
   * \sa getName(), getPathName()
   */
  std::string getPath() const;

  /*!
   * \brief Return full path of Group object, including its name.
   *
   * If a DataStore contains a Group tree structure a/b/c/d/e, the
   * following results are expected:
   *
   * Method Call      | Result
   * -----------------|----------
   * e->getName()     | e
   * e->getPath()     | a/b/c/d
   * e->getPathName() | a/b/c/d/e
   *
   * \sa getName(), getPath(), View::getPathName()
   */
  std::string getPathName() const
  {
    const auto path = getPath();

    if(path.length() < 1)
    {
      return getName();
    }

    return path + getPathDelimiter() + getName();
  }

  /*!
   * \brief Return pointer to non-const parent Group of a Group.
   *
   * Note that if this method is called on the root Group in a
   * DataStore, a pointer to itself is returned.
   * This allows root->getParent()->getParent() to always work similar
   * to how the filesystem's `cd /; cd ../..` works.
   */
  Group* getParent() { return m_parent; }

  /*!
   * \brief Return pointer to const parent Group of a Group.
   *
   * Note that if this method is called on the root Group in a
   * DataStore, a pointer to itself is returned.
   * This allows root->getParent()->getParent() to always work similar
   * to how the filesystem's `cd /; cd ../..` works.
   */
  const Group* getParent() const { return m_parent; }

  /*!
   * \brief Return number of child Groups in a Group object.
   */
  IndexType getNumGroups() const;

  /*!
   * \brief Return number of Views owned by a Group object.
   */
  IndexType getNumViews() const;

  /*!
   * \brief Return pointer to non-const DataStore object that owns this
   * object.
   */
  DataStore* getDataStore() { return m_datastore; }

  /*!
   * \brief Return pointer to const DataStore object that owns this
   * object.
   */
  const DataStore* getDataStore() const { return m_datastore; }

  /*!
   * \brief Return true if this Group is the DataStore's root Group.
   */
  bool isRoot() const { return m_parent == this; }

#ifdef AXOM_USE_UMPIRE

  /*!
   * \brief Return the ID of the default umpire::Allocator associated with this
   * Group.
   */
  int getDefaultAllocatorID() const { return m_default_allocator_id; }

  /*!
   * \brief Return the default umpire::Allocator associated with this Group.
   */
  umpire::Allocator getDefaultAllocator() const
  {
    umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
    return rm.getAllocator(m_default_allocator_id);
  }

  /*!
   * \brief Set the default umpire::Allocator associated with this Group.
   */
  Group* setDefaultAllocator(umpire::Allocator alloc)
  {
    m_default_allocator_id = alloc.getId();
    return this;
  }

  /*!
   * \brief Set the default umpire::Allocator associated with this Group.
   */
  Group* setDefaultAllocator(int allocId)
  {
    m_default_allocator_id = allocId;
    return this;
  }
#endif

  //@}

  //@{
  //!  @name View query methods.

  /*!
   * \brief Return true if Group includes a descendant View with
   * given name or path; else false.
   */
  bool hasView(const std::string& path) const;

  /*!
   * \brief Return true if this Group owns a View with given name (not path);
   * else false.
   */
  bool hasChildView(const std::string& name) const;

  /*!
   * \brief Return true if this Group owns a View with given index; else false.
   */
  bool hasView(IndexType idx) const;

  /*!
   * \brief Return index of View with given name owned by this Group object.
   *
   *        If no such View exists, return sidre::InvalidIndex;
   */
  IndexType getViewIndex(const std::string& name) const;

  /*!
   * \brief Return name of View with given index owned by Group object.
   *
   *        If no such View exists, return sidre::InvalidName.
   */
  const std::string& getViewName(IndexType idx) const;

  //@}

  //@{
  //!  @name View access methods.

  /*!

   * \brief Return pointer to non-const View with given name or path.
   *
   * This method requires that all groups in the path exist if a path is given.
   *
   * If no such View exists, nullptr is returned.
   */
  View* getView(const std::string& path);

  /*!
   * \brief Return pointer to const View with given name or path.
   *
   * This method requires that all Groups in the path exist if a path is given.
   *
   * If no such View exists, nullptr is returned.
   */
  const View* getView(const std::string& path) const;

  /*!
   * \brief Return pointer to non-const View with given index.
   *
   * If no such View exists, nullptr is returned.
   */
  View* getView(IndexType idx);

  /*!
   * \brief Return pointer to const View with given index.
   *
   * If no such View exists, nullptr is returned.
   */
  const View* getView(IndexType idx) const;

  //@}

  //@{
  //!  @name View iteration methods.
  //!
  //! Using these methods, a code can get the first View index and each
  //! succeeding index.  This allows View iteration using the same
  //! constructs in C++, C, and Fortran.  Example:
  //!
  //!      sidre::IndexType idx = grp->getFirstValidViewIndex();
  //!      while( sidre::indexIsValid(idx) )
  //!      {
  //!          View* view = grp->getView(idx);
  //!
  //!          /// code here using view
  //!
  //!          idx = grp -> getNextValidViewIndex(idx);
  //!      }

  /*!
   * \brief Return first valid View index in Group object
   *        (i.e., smallest index over all Views).
   *
   * sidre::InvalidIndex is returned if Group has no Views.
   *
   * \sa sidre::IndexType
   * \sa sidre::indexIsValid()
   */
  IndexType getFirstValidViewIndex() const;

  /*!
   * \brief Return next valid View index in Group object after given index
   *        (i.e., smallest index over all View indices larger than given one).
   *
   * sidre::InvalidIndex is returned if there is no valid index greater
   * than given one.
   *
   * \sa sidre::IndexType
   * \sa sidre::indexIsValid()
   */
  IndexType getNextValidViewIndex(IndexType idx) const;

  //@}

  //@{
  //!  @name Methods to create a View that has no associated data.
  //!
  //! \attention These methods do not allocate data or associate a View
  //! with data. Thus, to do anything useful with a View created by one
  //! of these methods, the View should be allocated, attached to a Buffer
  //! or attached to externally-owned data.
  //!
  //! Each of these methods is a no-op if the given View name is an
  //! empty string or the Group already has a View with given name or path.
  //!
  //! Additional conditions under which a method can be a no-op are described
  //! for each method.

  /*!
   * \brief Create an undescribed (i.e., empty) View object with given name
   * or path in this Group.
   *
   * If path is an empty string, an unnamed view can be created only if
   * this Group was created to hold items in a list format.  Otherwise
   * an empty string will result in a nullptr being returned.
   *
   * \return pointer to new View object or nullptr if one is not created.
   */
  View* createView(const std::string& path);

  /*!
   * \brief Create View object with given name or path in this Group that
   *  has a data description with data type and number of elements.
   *
   * If given data type is undefined, or given number of elements is < 0,
   * method is a no-op.
   *
   * \return pointer to new View object or nullptr if one is not created.
   */
  View* createView(const std::string& path, TypeID type, IndexType num_elems);

  /*!
   * \brief Create View object with given name or path in this Group that
   *  has a data description with data type and shape.
   *
   * If given data type is undefined, or given number of dimensions is < 0,
   * or given shape ptr is null, method is a no-op.
   *
   * \return pointer to new View object or nullptr if one is not created.
   */
  View* createView(const std::string& path,
                   TypeID type,
                   int ndims,
                   const IndexType* shape);

  /*!
   * \brief Create View object with given name or path in this Group that
   *  is described by a Conduit DataType object.
   *
   * \return pointer to new View object or nullptr if one is not created.
   */
  View* createView(const std::string& path, const DataType& dtype);

  //@}

  //@{
  //!  @name Methods to create a View with a Buffer attached.
  //!
  //! \attention The Buffer passed to each of these methods may or may not
  //! be allocated. Thus, to do anything useful with a View created by one
  //! of these methods, the Buffer must be allocated and it must be compatible
  //! with the View data description.
  //!
  //! Each of these methods is a no-op if Group already has a View or child
  //! Group with the given name or path.
  //!
  //! If this Group was created to hold items in list format, the path can
  //! be an empty string. Otherwise an empty string for the path will result
  //! in a no-op.
  //!
  //! Also, calling one of these methods with a null Buffer pointer is
  //! similar to creating a View with no data association.
  //!
  //! Additional conditions under which a method can be a no-op are described
  //! for each method.

  /*!
   * \brief Create an undescribed View object with given name or path in
   * this Group and attach given Buffer to it.
   *
   * \attention The View cannot be used to access data in Buffer until it
   * is described by calling a View::apply() method.
   *
   * This method is equivalent to:
   * group->createView(name)->attachBuffer(buff).
   *
   * \return pointer to new View object or nullptr if one is not created.
   *
   * \sa View::attachBuffer()
   */
  View* createView(const std::string& path, Buffer* buff);

  /*!
   * \brief Create View object with given name or path in this Group that
   * has a data description with data type and number of elements and
   * attach given Buffer to it.
   *
   * If given data type is undefined, or given number of elements is < 0,
   * method is a no-op.
   *
   * This method is equivalent to:
   * group->createView(name, type, num_elems)->attachBuffer(buff), or
   * group->createView(name)->attachBuffer(buff)->apply(type, num_elems).
   *
   * \return pointer to new View object or nullptr if one is not created.
   *
   * \sa View::attachBuffer()
   */
  View* createView(const std::string& path,
                   TypeID type,
                   IndexType num_elems,
                   Buffer* buff);

  /*!
   * \brief Create View object with given name or path in this Group that
   * has a data description with data type and shape and attach given
   * Buffer to it.
   *
   * If given data type is undefined, or given number of dimensions is < 0,
   * or given shape ptr is null, method is a no-op.
   *
   * This method is equivalent to:
   * group->createView(name, type, ndims, shape)->attachBuffer(buff), or
   * group->createView(name)->attachBuffer(buff)->apply(type, ndims, shape).
   *
   * \return pointer to new View object or nullptr if one is not created.
   *
   * \sa View::attachBuffer()
   */
  View* createView(const std::string& path,
                   TypeID type,
                   int ndims,
                   const IndexType* shape,
                   Buffer* buff);

  /*!
   * \brief Create View object with given name or path in this Group that
   *  is described by a Conduit DataType object and attach given Buffer to it.
   *
   * This method is equivalent to:
   * group->createView(name, dtype)->attachBuffer(buff), or
   * group->createView(name)->attachBuffer(buff)->apply(dtype).
   *
   * \return pointer to new View object or nullptr if one is not created.
   *
   * \sa View::attachBuffer()
   */
  View* createView(const std::string& path, const DataType& dtype, Buffer* buff);

  //@}

  //@{
  //!  @name Methods to create a View with externally-owned data attached.
  //!
  //! \attention To do anything useful with a View created by one of these
  //! methods, the external data must be allocated and compatible with the
  //! View description.
  //!
  //! Each of these methods is a no-op if the given View name is an
  //! empty string or the Group already has a View with given name or path.
  //!
  //! Additional conditions under which a method can be a no-op are described
  //! for each method.

  /*!
   * \brief Create View object with given name with given name or path in
   * this Group and attach external data ptr to it.
   *
   * \attention Note that the View is "opaque" (it has no knowledge of
   * the type or structure of the data) until a View::apply() method
   * is called.
   *
   * This method is equivalent to:
   * group->createView(name)->setExternalDataPtr(external_ptr).
   *
   * \return pointer to new View object or nullptr if one is not created.
   *
   * \sa View::setExternalDataPtr()
   */
  View* createView(const std::string& path, void* external_ptr);

  /*!
   * \brief Create View object with given name or path in this Group that
   * has a data description with data type and number of elements and
   * attach externally-owned data to it.
   *
   * If given data type is undefined, or given number of elements is < 0,
   * method is a no-op.
   *
   * This method is equivalent to:
   * group->createView(name, type, num_elems)->setExternalDataPtr(external_ptr),
   * or group->createView(name)->setExternalDataPtr(external_ptr)->
   *           apply(type, num_elems).
   *
   * \return pointer to new View object or nullptr if one is not created.
   *
   * \sa View::setExternalDataPtr()
   */
  View* createView(const std::string& path,
                   TypeID type,
                   IndexType num_elems,
                   void* external_ptr);

  /*!
   * \brief Create View object with given name or path in this Group that
   * has a data description with data type and shape and attach
   * externally-owned data to it.
   *
   * If given data type is undefined, or given number of dimensions is < 0,
   * or given shape ptr is null, method is a no-op.
   *
   * This method is equivalent to:
   * group->createView(name, type, ndims, shape)->
   *        setExternalDataPtr(external_ptr), or
   * group->createView(name)->setExternalDataPtr(external_ptr)->
   *        apply(type, ndims, shape).
   *
   * \return pointer to new View object or nullptr if one is not created.
   *
   * \sa View::setExternalDataPtr()
   */
  View* createView(const std::string& path,
                   TypeID type,
                   int ndims,
                   const IndexType* shape,
                   void* external_ptr);
  /*!
   * \brief Create View object with given name or path in this Group that
   * is described by a Conduit DataType object and attach externally-owned
   * data to it.
   *
   * This method is equivalent to:
   * group->createView(name, dtype)->setExternalDataPtr(external_ptr), or
   * group->createView(name)->setExternalDataPtr(external_ptr)->apply(dtype).
   *
   * \return pointer to new View object or nullptr if one is not created.
   *
   * \sa View::attachBuffer()
   */
  View* createView(const std::string& path,
                   const DataType& dtype,
                   void* external_ptr);

  //@}

  //@{
  //!  @name Methods to create a View and allocate data for it.
  //!
  //! Each of these methods is a no-op if the given View name is an
  //! empty string or the Group already has a View with given name or path.
  //!
  //! Additional conditions under which a method can be a no-op are described
  //! for each method.

  /*!
   * \brief Create View object with given name or path in this Group that
   * has a data description with data type and number of elements and
   * allocate data for it.
   *
   * If given data type is undefined, or given number of elements is < 0,
   * method is a no-op.
   *
   * This is equivalent to: createView(name)->allocate(type, num_elems), or
   * createView(name, type, num_elems)->allocate()
   *
   * \return pointer to new View object or nullptr if one is not created.
   *
   * \sa View::allocate()
   */
  View* createViewAndAllocate(const std::string& path,
                              TypeID type,
                              IndexType num_elems,
                              int allocID = INVALID_ALLOCATOR_ID);

  /*!
   * \brief Create View object with given name or path in this Group that
   * has a data description with data type and shape and allocate data for it.
   *
   * If given data type is undefined, or given number of dimensions is < 0,
   * or given shape ptr is null, method is a no-op.
   *
   * This method is equivalent to:
   * group->createView(name)->allocate(type, ndims, shape), or
   * createView(name, type, ndims, shape)->allocate().
   *
   * \return pointer to new View object or nullptr if one is not created.
   *
   * \sa View::allocate()
   */
  View* createViewAndAllocate(const std::string& path,
                              TypeID type,
                              int ndims,
                              const IndexType* shape,
                              int allocID = INVALID_ALLOCATOR_ID);

  /*!
   * \brief Create View object with given name or path in this Group that
   * is described by a Conduit DataType object and allocate data for it.
   *
   * This method is equivalent to:
   * group->createView(name)->allocate(dtype), or
   * group->createView(name, dtype)->allocate().
   *
   * If given data type object is empty, data will not be allocated.
   *
   * \return pointer to new View object or nullptr if one is not created.
   *
   * \sa View::allocate()
   */
  View* createViewAndAllocate(const std::string& path,
                              const DataType& dtype,
                              int allocID = INVALID_ALLOCATOR_ID);

  /*!
   * \brief Create View object with given name or path in this Group
   * set its data to given scalar value.
   *
   * This is equivalent to: createView(name)->setScalar(value);
   *
   * If given data type object is empty, data will not be allocated.
   *
   * \return pointer to new View object or nullptr if one is not created.
   *
   * \sa View::setScalar()
   */
  template <typename ScalarType>
  View* createViewScalar(const std::string& path, ScalarType value)
  {
    View* view = createView(path);
    if(view != nullptr)
    {
      view->setScalar(value);
    }

    return view;
  }

  /*!
   * \brief Create View object with given name or path in this Group
   * set its data to given string.
   *
   * This is equivalent to: createView(name)->setString(value);
   *
   * If given data type object is empty, data will not be allocated.
   *
   * \return pointer to new View object or nullptr if one is not created.
   *
   * \sa View::setString()
   */
  View* createViewString(const std::string& path, const std::string& value);

  //@}

  //@{
  //!  @name View destruction methods.
  //!
  //! Each of these methods is a no-op if the specified View does not exist.

  /*!
   * \brief Destroy View with given name or path owned by this Group, but leave
   * its data intect.
   */
  void destroyView(const std::string& path);

  /*!
   * \brief Destroy View with given index owned by this Group, but leave
   * its data intect.
   */
  void destroyView(IndexType idx);

  /*!
   * \brief Destroy all Views owned by this Group, but leave all their
   *        data intact.
   */
  void destroyViews();

  /*!
   * \brief Destroy View with given name or path owned by this Group and
   * deallocate
   * its data if it's the only View associated with that data.
   */
  void destroyViewAndData(const std::string& path);

  /*!
   * \brief Destroy View with given index owned by this Group and deallocate
   * its data if it's the only View associated with that data.
   */
  void destroyViewAndData(IndexType idx);

  /*!
   * \brief Destroy all Views owned by this Group and deallocate
   * data for each View when it's the only View associated with that data.
   */
  void destroyViewsAndData();

  //@}

  //@{
  //!  @name View move and copy methods.

  /*!
   * \brief Remove given View object from its owning Group and move it
   *        to this Group.
   *
   * If given View pointer is null or Group already has a View with
   * same name as given View, method is a no-op.
   *
   * \return pointer to given argument View object or nullptr if View
   * is not moved into this Group.
   */
  View* moveView(View* view);

  /*!
   * \brief Create a copy of given View object and add it to this Group.
   *
   * Note that View copying is a "shallow" copy; the data associated with
   * the View is not copied. The new View object is associated with
   * the same data as the original.
   *
   * If given Group pointer is null or Group already has a child Group with
   * same name as given Group, method is a no-op.
   *
   * \return pointer to given argument Group object or nullptr if Group
   * is not moved into this Group.
   */
  View* copyView(View* view);

  //@}

  //@{
  //!  @name Child Group query methods.

  /*!
   * \brief Return true if this Group has a descendant Group with given
   * name or path; else false.
   */
  bool hasGroup(const std::string& path) const;

  /*!
   * \brief Return true if this Group has a child Group with given
   * name; else false.
   */
  bool hasChildGroup(const std::string& name) const;

  /*!
   * \brief Return true if Group has an immediate child Group
   *        with given index; else false.
   */
  bool hasGroup(IndexType idx) const;

  /*!
   * \brief Return the index of immediate child Group with given name.
   *
   *        If no such child Group exists, return sidre::InvalidIndex;
   */
  IndexType getGroupIndex(const std::string& name) const;

  /*!
   * \brief Return the name of immediate child Group with given index.
   *
   *        If no such child Group exists, return sidre::InvalidName.
   */
  const std::string& getGroupName(IndexType idx) const;

  //@}

  //@{
  //!  @name Group access and iteration methods.

  /*!
   * \brief Return pointer to non-const child Group with given name or path.
   *
   * This method requires that all Groups in the path exist if a path is given.
   *
   * If no such Group exists, nullptr is returned.
   */
  Group* getGroup(const std::string& path);

  /*!
   * \brief Return pointer to const child Group with given name or path.
   *
   * This method requires that all Groups in the path exist if a path is given.
   *
   * If no such Group exists, nullptr is returned.
   */
  Group const* getGroup(const std::string& path) const;

  /*!
   * \brief Return pointer to non-const immediate child Group with given index.
   *
   * If no such Group exists, nullptr is returned.
   */
  Group* getGroup(IndexType idx);

  /*!
   * \brief Return pointer to const immediate child Group with given index.
   *
   * If no such Group exists, nullptr is returned.
   */
  const Group* getGroup(IndexType idx) const;

  //@}

  //@{
  //!  @name Group iteration methods.
  //!
  //! Using these methods, a code can get the first Group index and each
  //! succeeding index.  This allows Group iteration using the same
  //! constructs in C++, C, and Fortran.  Example:
  //!      for (sidre::IndexType idx = grp->getFirstValidGroupIndex();
  //!           sidre::indexIsValid(idx);
  //!           idx = grp->getNextValidGroupIndex(idx))
  //!      {
  //!          Group * group = grp->getGroup(idx);
  //!
  //!          /// code here using group
  //!      }

  /*!
   * \brief Return first valid child Group index (i.e., smallest
   *        index over all child Groups).
   *
   * sidre::InvalidIndex is returned if Group has no child Groups.
   *
   * \sa axom::sidre::indexIsValid()
   */
  IndexType getFirstValidGroupIndex() const;

  /*!
   * \brief Return next valid child Group index after given index
   *        (i.e., smallest index over all child Group indices larger
   *        than given one).
   *
   * sidre::InvalidIndex is returned if there is no valid index greater
   * than given one.
   *
   * \sa axom::sidre::indexIsValid()
   */
  IndexType getNextValidGroupIndex(IndexType idx) const;

  //@}

  //@{
  //!  @name Child Group creation and destruction methods.

  /*!
   * \brief Create a child Group within this Group with given name or path.
   *
   * If name is an empty string or Group already has a child Group with
   * given name or path, method is a no-op.
   *
   * The optional is_list argument is used to determine if the created
   * child Group will hold items in list format.
   *
   * \return pointer to created Group object or nullptr if new
   * Group is not created.
   */
  Group* createGroup(const std::string& path, bool is_list = false);

  /*
   * \brief Create a child Group within this Group with no name.
   *
   * This is intended only to be called when this Group holds items in list
   * format.  If this Group does not use list format, this method is a
   * no-op.
   *
   * \return pointer to created Group object or nullptr if new Group is
   * not created.
   */
  Group* createUnnamedGroup(bool is_list = false);

  /*!
   * \brief Destroy child Group in this Group with given name or path.
   *
   * If no such Group exists, method is a no-op.
   */
  void destroyGroup(const std::string& path);

  /*!
   * \brief Destroy child Group within this Group with given index.
   *
   * If no such Group exists, method is a no-op.
   */
  void destroyGroup(IndexType idx);

  /*!
   * \brief Destroy all child Groups in this Group.
   *
   * Note that this will recursively destroy entire Group sub-tree below
   * this Group.
   */
  void destroyGroups();

  //@}

  //@{
  //!  @name Group move and copy methods.

  /*!
   * \brief Remove given Group object from its parent Group and make it
   *        a child of this Group.
   *
   * If given Group pointer is null or Group already has a child Group with
   * same name as given Group, method is a no-op.
   *
   * \return pointer to given argument Group object or nullptr if Group
   * is not moved into this Group.
   */
  Group* moveGroup(Group* group);

  /*!
   * \brief Create a copy of Group hierarchy rooted at given Group and make it
   *        a child of this Group.
   *
   * Note that all Views in the Group hierarchy are copied as well.
   *
   * Note that Group copying is a "shallow" copy; the data associated
   * with Views in a Group are not copied. In particular, the new Group
   * hierarchy and all its Views is associated with the same data as the
   * given Group.
   *
   * If given Group pointer is null or Group already has a child Group with
   * same name as given Group, method is a no-op.
   *
   * \return pointer to given argument Group object or nullptr if Group
   * is not moved into this Group.
   */
  Group* copyGroup(Group* group);

  //@}

  //@{
  //!  @name Group print methods.

  /*!
   * \brief Print JSON description of data Group to stdout.
   *
   * Note that this will recursively print entire Group sub-tree
   * starting at this Group object.
   */
  void print() const;

  /*!
   * \brief Print JSON description of data Group to an ostream.
   *
   * Note that this will recursively print entire Group sub-tree
   * starting at this Group object.
   */
  void print(std::ostream& os) const;

  /*!
   * \brief Print given number of levels of Group sub-tree
   *        starting at this Group object to an output stream.
   */
  void printTree(const int nlevels, std::ostream& os) const;

  //@}

  /*!
   * \brief Copy description of Group hierarchy rooted at this Group to
   * given Conduit node.
   *
   * The description includes Views of this Group and all of its children
   * recursively.
   */
  void copyToConduitNode(Node& n) const;

  /*!
   * \brief Copy data Group native layout to given Conduit node.
   *
   * The native layout is a Conduit Node hierarchy that maps the Conduit Node
   * data
   * externally to the Sidre View data so that it can be filled in from the data
   * in the file (independent of file format) and can be accessed as a Conduit
   * tree.
   *
   * \return True if the Group or any of its children were added to the Node,
   * false otherwise.
   *
   */
  bool createNativeLayout(Node& n, const Attribute* attr = nullptr) const;

  /*!
   * \brief Copy data Group external layout to given Conduit node.
   *
   * The external layout is a Conduit Node hierarchy that maps the Conduit Node
   * data externally to the Sidre View data so that it can be filled in from the
   * data in the file (independent of file format) and can be accessed as a
   * Conduit tree.
   *
   * Only the Views which have external data are added to the node.
   *
   * \return True if the Group or any of its children have an external
   *  View, false otherwise.
   */
  bool createExternalLayout(Node& n, const Attribute* attr = nullptr) const;

  /*!
   * \brief Return true if this Group is equivalent to given Group; else false.
   *
   * Two Groups are equivalent if they are the root Groups of identical
   * Group hierarchy structures with the same names for all Views and
   * Groups in the hierarchy, and the Views are also equivalent.
   *
   * \param other     The group to compare with
   * \param checkName If true (default), groups must have the same name
   *                  to be equivalent.  If false, disregard the name.
   *
   * \sa View::isEquivalentTo()
   */
  bool isEquivalentTo(const Group* other, bool checkName = true) const;

  /*!
   * \brief Return true if this Group holds items in map format.
   */
  bool isUsingMap() const { return !m_is_list; }

  /*!
   * \brief Return true if this Group holds items in list format.
   */
  bool isUsingList() const { return m_is_list; }

  //@{
  /*!
 * @name    Group I/O methods
 *   These methods save and load Group trees to and from files.
 *   This includes the views and buffers used in by groups in the tree.
 *   We provide several "protocol" options:
 *
 *   protocols:
 *    sidre_hdf5 (default when Axom is configured with hdf5)
 *    sidre_conduit_json (default otherwise)
 *    sidre_json
 *
 *    conduit_hdf5
 *    conduit_bin
 *    conduit_json
 *    json
 *
 *   \note The sidre_hdf5 and conduit_hdf5 protocols are only available
 *   when Axom is configured with hdf5.
 *
 *   There are two overloaded versions for each of save, load, and
 *   loadExternalData.  The first of each takes a file path and is intended
 *   for use in a serial context and can be called directly using any
 *   of the supported protocols.  The second takes an hdf5 handle that
 *   has previously been created by the calling code.  These mainly exist
 *   to handle parallel I/O calls from the SPIO component.  They can only
 *   take the sidre_hdf5 or conduit_hdf5 protocols.
 *
 *   \note The hdf5 overloads are only available when Axom is configured
 *   with hdf5.
 */

  /*!
   * \brief Save the Group to a file.
   *
   *  Saves the tree starting at this Group and the Buffers used by the Views
   *  in this tree.
   *
   *  If attr is a null pointer, dump all Views.  Otherwise, only dump Views
   *  which have the Attribute set.
   *
   * \param path      file path
   * \param protocol  I/O protocol
   * \param attr      Save Views that have Attribute set.
   */
  void save(const std::string& path,
            const std::string& protocol = SIDRE_DEFAULT_PROTOCOL,
            const Attribute* attr = nullptr) const;

  /*!
   * \brief Load a Group hierarchy from a file into this Group
   *
   * This method instantiates the Group hierarchy and its Views stored
   * in the file under this Group.  The name of this Group is not
   * changed.
   *
   * If preserve_contents is true, then the names of the children held by the
   * Node cannot be the same as the names of the children already held by this
   * Group.  If there is a naming conflict, an error will occur.
   *
   * \param path     file path
   * \param protocol I/O protocol
   * \param preserve_contents  If true, any child Groups and Views held by
   *                           this Group remain in place.  If false, all
   *                           child Groups and Views are destroyed before
   *                           loading data from the file.
   */
  void load(const std::string& path,
            const std::string& protocol = SIDRE_DEFAULT_PROTOCOL,
            bool preserve_contents = false);

  /*!
   * \brief Load a Group hierarchy from a file into this Group, reporting
   *        the Group name stored in the file
   *
   * This method instantiates the Group hierarchy and its Views stored
   * in the file under this Group.  The name of this Group is not
   * changed.  The name of the group stored in the file is returned in
   * the output parameter name_from_file.  This can be used to rename the
   * group in a subsequent call.
   *
   * If preserve_contents is true, then the names of the children held by the
   * Node cannot be the same as the names of the children already held by this
   * Group.  If there is a naming conflict, an error will occur.
   *
   * \param [in]  path     file path to load
   * \param [in]  protocol I/O protocol to use
   * \param [in]  preserve_contents  If true, any child Groups and Views
   *                           held by this Group remain in place.
   *                           If false, all child Groups and Views are
   *                           destroyed before loading data from the file.
   * \param [out] name_from_file    Group name stored in the file
   */
  void load(const std::string& path,
            const std::string& protocol,
            bool preserve_contents,
            std::string& name_from_file);

  /*!
   * \brief Create a child Group and load a Group hierarchy from file
   *        into the new Group.
   *
   * This is a convenience routine for the following sequence:
   * - create a group with name or path group_name
   * - load a Group hierarchy from a file into the newly-created Group
   * - return the newly created Group, or nullptr if creation failed
   * - out-parameters return:
   *   - the group name from the file
   *   - a flag indicating success reading the file
   *
   * As with the createGroup() method, if group_name is empty or there
   * already exists a child Group with that name or path, the child Group
   * will not be created and this method will return nullptr.
   *
   * As with the load() method, after calling createGroupAndLoad() a host
   * code may choose to rename the newly-created Group with the string
   * returned in group_name.
   *
   * \param [in,out] group_name    In: name for the new group.
   *                               Out: the group name stored in the file.
   * \param [in]     path          file path
   * \param [in]     protocol      I/O protocol
   * \param [out]    load_success  Report success of the load operation
   *
   * \return pointer to created Group object or nullptr if new
   *         Group is not created.
   */
  Group* createGroupAndLoad(std::string& group_name,
                            const std::string& path,
                            const std::string& protocol,
                            bool& load_success);

  /*!
   * \brief Load data into the Group's external views from a file.
   *
   * No protocol argument is needed, as this only is used with the sidre_hdf5
   * protocol.
   *
   * \param path      file path
   */
  void loadExternalData(const std::string& path);

#ifdef AXOM_USE_HDF5

  /*!
   * \brief Save the Group to an hdf5 handle.
   *
   *  If attr is nullptr, dump all Views.  Otherwise, only dump Views
   *  which have the Attribute set.
   *
   * \param h5_id      hdf5 handle
   * \param protocol   I/O protocol sidre_hdf5 or conduit_hdf5
   * \param attr       Save Views that have Attribute set.
   */
  void save(const hid_t& h5_id,
            const std::string& protocol = SIDRE_DEFAULT_PROTOCOL,
            const Attribute* attr = nullptr) const;

  /*!
   * \brief Load the Group from an hdf5 handle.
   *
   * If preserve_contents is true, then the names of the children held by the
   * Node cannot be the same as the names of the children already held by this
   * Group.  If there is a naming conflict, an error will occur.
   *
   * \param h5_id      hdf5 handle
   * \param protocol   I/O protocol sidre_hdf5 or conduit_hdf5
   * \param preserve_contents  If true, any child Groups and Views held by
   *                           this Group remain in place.  If false, all
   *                           child Groups and Views are destroyed before
   *                           loading data from the file.
   */
  void load(const hid_t& h5_id,
            const std::string& protocol = SIDRE_DEFAULT_PROTOCOL,
            bool preserve_contents = false);

  /*!
   * \brief Load the Group from an hdf5 handle.
   *
   * If preserve_contents is true, then the names of the children held by the
   * Node cannot be the same as the names of the children already held by this
   * Group.  If there is a naming conflict, an error will occur.
   *
   * \param [in]  h5_id      hdf5 handle
   * \param [in]  protocol   I/O protocol sidre_hdf5 or conduit_hdf5
   * \param [in]  preserve_contents  If true, any child Groups and Views held by
   *                           this Group remain in place.  If false, all
   *                           child Groups and Views are destroyed before
   *                           loading data from the file.
   * \param [out] name_from_file    Group name stored in the file
   */
  void load(const hid_t& h5_id,
            const std::string& protocol,
            bool preserve_contents,
            std::string& name_from_file);

  /*!
   * \brief Load data into the Group's external views from a hdf5 handle.
   *
   * No protocol argument is needed, as this only is used with the sidre_hdf5
   * protocol.
   *
   * \param h5_id      hdf5 handle
   */
  void loadExternalData(const hid_t& h5_id);

#endif /* AXOM_USE_HDF5 */

  //@}

  /*!
   * \brief Change the name of this Group.
   *
   * The name of this group is changed to the new name.  If this group
   * has a parent group, the name for this group held by the parent is
   * also changed.
   *
   * Warnings will occur and the name will not be changed under these
   * conditions:  If the new name is an empty string, if the new name
   * contains a path delimiter (usually '/'), or if the new name is
   * identical to a name that is already held by the parent for another
   * Group or View object.
   *
   * It is possible to rename the root Group, but a code cannot
   * subsequently rename root Group back to its original empty string
   * name.
   *
   * \param new_name    The new name for this group.
   *
   * \return            Success or failure of rename.
   */
  bool rename(const std::string& new_name);

  /*!
   * \brief Import data from a conduit Node into a Group
   *
   * This imports the hierarchy from the Node into a Sidre Group with the
   * same tree structure.
   *
   * This does not support conduit's list datatype.  If the Node contains a
   * list any where in its tree, the list and any child Nodes descending from
   * the list will not be imported.  A warning will occur and an unsuccessful
   * return value will be returned.
   *
   * If preserve_contents is true, then the names of the children held by the
   * Node cannot be the same as the names of the children already held by this
   * Group.  If there is a naming conflict, an error will occur.
   *
   * \param node               A conduit Node containing hierarchical data.
   * \param preserve_contents  If true, any child Groups and Views held by
   *                           this Group remain in place.  If false, all
   *                           child Groups and Views are destroyed before
   *                           importing data from the Node.
   *
   * \return                   true for success, false if the full conduit
   *                           tree is not succesfully imported.
   */
  bool importConduitTree(const conduit::Node& node,
                         bool preserve_contents = false);

  /*!
   * \brief Import data from a conduit Node into a Group without copying arrays
   *
   * This differs from the importConduitTree in that it does not copy any
   * data held by the Node as an array.  Instead it imports the existing
   * pointer to the array as an external pointer.
   *
   * This imports the hierarchy from the Node into a Sidre Group with the
   * same tree structure.
   *
   * This does not support conduit's list datatype.  If the Node contains a
   * list any where in its tree, the list and any child Nodes descending from
   * the list will not be imported.  A warning will occur and an unsuccessful
   * return value will be returned.
   *
   * If preserve_contents is true, then the names of the children held by the
   * Node cannot be the same as the names of the children already held by this
   * Group.  If there is a naming conflict, an error will occur.
   *
   * \param node               A conduit Node containing hierarchical data.
   * \param preserve_contents  If true, any child Groups and Views held by
   *                           this Group remain in place.  If false, all
   *                           child Groups and Views are destroyed before
   *                           importing data from the Node.
   *
   * \return                   true for success, false if the full conduit
   *                           tree is not succesfully imported.
   */
  bool importConduitTreeExternal(conduit::Node& node,
                                 bool preserve_contents = false);

private:
  DISABLE_DEFAULT_CTOR(Group);
  DISABLE_COPY_AND_ASSIGNMENT(Group);
  DISABLE_MOVE_AND_ASSIGNMENT(Group);

  //@{
  //!  @name Private Group ctors and dtors
  //!        (callable only by DataStore and Group methods).

  /*!
   *  \brief Private ctor that creates a Group with given name
   *         in the given DataStore.
   *
   *  attachGroup must be called on a newly created Group to insert it
   *  into the hierarchy. The root group is an exception to this rule.
   *
   *  The boolean argument is_list, if true, allows the Group to hold its
   *  child items in list format, which allows those items to have empty
   *  strings for names.  If not in list format, all items must have unique
   *  non-empty strings for names.
   */
  Group(const std::string& name, DataStore* datastore, bool is_list);

  /*!
   * \brief Destructor destroys all Views and child Groups.
   */
  ~Group();

  //@}

  //@{
  //!  @name View attach and detach methods.

  /*!
   * \brief Attach View object to this Group.
   */
  View* attachView(View* view);

  /*!
   * \brief Detach View object from this Group.
   */
  View* detachView(const View* view) { return detachView(view->getName()); }

  /*!
   * \brief Detach View with given name from this Group.
   */
  View* detachView(const std::string& name);

  /*!
   * \brief Detach View with given index from this Group.
   */
  View* detachView(IndexType idx);

  //@}

  //@{
  //!  @name Group attach and detach methods.

  /*!
   * \brief Attach Group object to this Group.
   */
  Group* attachGroup(Group* view);

  /*!
   * \brief Detach Child Group with given name from this Group.
   */
  Group* detachGroup(const std::string& name);

  /*!
   * \brief Detach Child Group with given index from this Group.
   */
  Group* detachGroup(IndexType idx);

  //@}

  //@{
  //!  @name Private Group View manipulation methods.

  /*!
   * \brief Destroy View and its data if its data is not shared with any
   * other View.
   *
   * Data will not be destroyed as long as a View still exists that
   * references it.
   *
   * \attention this method assumes View is owned by this Group.
   */
  void destroyViewAndData(View* view);

  //@}

  //@{
  //!  @name Private Group methods for interacting with Conduit Nodes.

  /*!
   * \brief Private method to copy Group to Conduit Node.
   *
   * Note: This is for the "sidre_{zzz}" protocols.
   *
   * \return True if the group or any of its children have saved Views,
   * false otherwise.
   */
  bool exportTo(conduit::Node& result, const Attribute* attr) const;

  /*!
   * \brief Private method to copy Group to Conduit Node.
   *
   * \param buffer_indices Used to track what Buffers are referenced
   * by the Views in this Group and Groups in the sub-tree below it.
   *
   * \return True if the group or any of its children have saved Views,
   * false otherwise.
   */
  bool exportTo(conduit::Node& data_holder,
                const Attribute* attr,
                std::set<IndexType>& buffer_indices) const;

  /*!
   * \brief Private method to build a Group hierarchy from Conduit Node.
   *
   * Note: This is for the "sidre_{zzz}" protocols.
   */
  void importFrom(conduit::Node& node, bool preserve_contents = false);

  /*!
   * \brief Private method to copy Group from Conduit Node.
   *
   * Map of Buffer indices tracks old Buffer ids in the file to the
   * new Buffer ids in the datastore.  Buffer ids are not guaranteed
   * to remain the same when a tree is restored.
   *
   */
  void importFrom(conduit::Node& node,
                  const std::map<IndexType, IndexType>& buffer_id_map);

  //@}

  /*!
   * \brief Private method that returns the Group that is the next-to-last
   * entry in a slash-delimited ("/") path string.
   *
   * The string before the last "/" character, if there is one, is the
   * next-to-last path entry. In this case, the return value is that Group
   * in the path.
   *
   * If there is no "/" in the given path, the entire string is considered
   * the next-to-last path entry. In this case, the return value is this
   * Group.
   *
   * The path argument is modified while walking the path. Its value when
   * the method returns is the last entry in the path, either the string
   * following the last "/" in the input (if there is one) or the entire
   * input path string if it contains no "/".
   */
  Group* walkPath(std::string& path, bool create_groups_in_path);

  /*!
   * \brief Const private method that returns the Group that is the
   * next-to-last entry in a delimited path string.
   *
   * The path argument is modified while walking the path. Its value when
   * the method returns is the last entry in the path, either the string
   * following the last "/" in the input (if there is one) or the entire
   * input path string if it contains no "/".
   */
  const Group* walkPath(std::string& path) const;

  /*!
   * \brief Private method. If allocatorID is a valid allocator ID then return
   *  it. Otherwise return the ID of the default allocator of the owning group.
   */
  int getValidAllocatorID(int allocatorID);

  /// Name of this Group object.
  std::string m_name;

  /// Index of this Group object within m_parent.
  IndexType m_index;

  /// Parent Group of this Group object.
  Group* m_parent;

  /// This Group object lives in the tree of this DataStore object.
  DataStore* m_datastore;

  /// This identifies whether this Group holds items in list format.
  bool m_is_list;

  /// Character used to denote a path string passed to get/create calls.
  AXOM_EXPORT static const char s_path_delimiter;

  ///////////////////////////////////////////////////////////////////
  //
  using ViewCollection = ItemCollection<View>;
  //
  using GroupCollection = ItemCollection<Group>;
  ///////////////////////////////////////////////////////////////////

  /// Collection of Views
  ViewCollection* m_view_coll;

  /// Collection of child Groups
  GroupCollection* m_group_coll;

#ifdef AXOM_USE_UMPIRE
  int m_default_allocator_id;
#endif
};

} /* end namespace sidre */
} /* end namespace axom */

#endif /* SIDRE_GROUP_HPP_ */

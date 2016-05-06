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
 * \brief   Header file containing definition of DataGroup class.
 *
 ******************************************************************************
 */

#ifndef DATAGROUP_HPP_
#define DATAGROUP_HPP_

// Standard C++ headers
#include <memory>
#include <map>
#include <string>
#include <vector>
#include <set>

#ifndef USE_UNORDERED_MAP
#define USE_UNORDERED_MAP
#endif
//#ifndef USE_DENSE_HASH_MAP
//#define USE_DENSE_HASH_MAP
//#endif

#ifndef USE_NEW_MAP_COLLECTION
#define USE_NEW_MAP_COLLECTION
#endif

#if defined(USE_UNORDERED_MAP)
//STL or Boost unordered_map, depending on
#if defined(USE_CXX11)
#include <unordered_map>
#else
#include "boost/unordered_map.hpp"
#endif
#endif

#if defined(USE_DENSE_HASH_MAP)
#include <sparsehash/dense_hash_map>
#endif

// Other CS Toolkit headers
#include "slic/slic.hpp"
#include "common/CommonTypes.hpp"

// SiDRe project headers
#include "SidreTypes.hpp"
#include "Collections.hpp"
#include "DataView.hpp"


namespace asctoolkit
{
namespace sidre
{

// using directives to make Conduit usage easier and less visible
using conduit::Node;
using conduit::Schema;

class DataBuffer;
class DataGroup;
class DataStore;

/*!
 * \class DataGroup
 *
 * \brief DataGroup holds a collection of DataViews and (child) DataGroups.
 *
 * The DataGroup class has the following properties:
 *
 *    - DataGroups can be organized into a (tree) hierachy by creating
 *      child groups from the root group owned by a DataStore object.
 *    - A DataGroup object can only be created by another DataGroup; the
 *      DataGroup ctor is not visible externally. A DataGroup is owned
 *      by the DataGroup that creates it (its parent) and becomes a child
 *      group of the parent.
 *    - A DataGroup object has a unique name (string) within its parent
 *      DataGroup.
 *    - A DataGroup object maintains a pointer to its parent DataGroup.
 *    - A DataGroup object can be moved or copied to another DataGroup.
 *    - DataGroup objects can create DataView objects within them. The
 *      DataGroup that creates a DataView owns it.
 *    - A DataView object has a unique name (string) within the DataGroup
 *      that owns it.
 *    - A DataView object can be moved or copied to another DataGroup.
 *
 * Note that DataViews and child DataGroups within a group can be accessed
 * by name or index.
 *
 * Note that certain methods for creating, accessing, etc. DataGroups and 
 * DataViews that take a string name accept either the name of a child group 
 * or view within a group object or a path syntax. When a path is given, the 
 * last item in the path indicates the item to be created, accessed, etc. So, 
 * for example,
 * 
 * \verbatim
 *
 *    DataView* view = group->createView("foo/bar/baz");
 *
 *    is equivalent to:
 *
 *    DataView* view = 
 *      group->createGroup("foo")->createGroup("bar")->createView("baz"); 
 *
 * \endverbatim
 *
 * In particular, intermediate groups "foo" and "bar" will be created in 
 * this case if they don't already exist.
 *
 * IMPORTANT: when views or groups are created, destroyed, copied, or moved,
 * indices of other views and groups in associated DataGroup objects may
 * become invalid. This is analogous to iterator invalidation for STL
 * containers when the container contents change.
 *
 */
class DataGroup
{
public:

  //
  // Friend declarations to constrain usage via controlled access to
  // private members.
  //
  friend class DataStore;


//@{
//!  @name Basic query and accessor methods.

  /*!
   * \brief Return const reference to name of group object.
   */
  const std::string& getName() const
  {
    return m_name;
  }

  /*!
   * \brief Return pointer to non-const parent group of a group.
   *
   * Note that if this method is called on the root group in a
   * DataStore, ATK_NULLPTR is returned.
   */
  DataGroup * getParent()
  {
    return m_parent;
  }

  /*!
   * \brief Return pointer to const parent group of a group.
   *
   * Note that if this method is called on the root group in a
   * DataStore, ATK_NULLPTR is returned.
   */
  const DataGroup * getParent() const
  {
    return m_parent;
  }

  /*!
   * \brief Return number of child groups in a group object.
   */
  size_t getNumGroups() const
  {
    return m_group_coll.getNumItems();
  }

  /*!
   * \brief Return number of views owned by a group object.
   */
  size_t getNumViews() const
  {
    return m_view_coll.getNumItems();
  }

  /*!
   * \brief Return pointer to non-const DataStore object that owns a 
   * object.
   */
  DataStore * getDataStore()
  {
    return m_datastore;
  }

  /*!
   * \brief Return pointer to const DataStore object that owns a 
   * object.
   */
  const DataStore * getDataStore() const
  {
    return m_datastore;
  }

//@}


//@{
//!  @name View query methods.

  /*!
   * \brief Return true if group owns a view with given name or path; 
   * else false.
   */
  bool hasView( const std::string& name ) const
  {
    return m_view_coll.hasItem(name);
  }

  /*!
   * \brief Return true if group owns a view with given index; else false.
   */
  bool hasView( IndexType idx ) const
  {
    return m_view_coll.hasItem(idx);
  }

  /*!
   * \brief Return index of view with given name owned by group object.
   *
   *        If no such view exists, return sidre::InvalidIndex;
   */
  IndexType getViewIndex(const std::string& name) const
  {
    SLIC_CHECK_MSG(hasView(name),
                   "Group " << this->getName() << 
                   " has no view with name '" << name << "'");

    return m_view_coll.getItemIndex(name);
  }

  /*!
   * \brief Return name of view with given index owned by group object.
   *
   *        If no such view exists, return sidre::InvalidName.
   */
  const std::string& getViewName(IndexType idx) const
  {
    SLIC_CHECK_MSG(hasView(idx),
                   "Group " << this->getName() << 
                   " has no view with index " << idx);

    return m_view_coll.getItemName(idx);
  }

//@}


//@{
//!  @name View access and iteration methods.

  /*!

   * \brief Return pointer to non-const view with given name or path.
   *
   * Thie method requires that all groups in the path exist if a path is given.
   *
   * If no such view exists, ATK_NULLPTR is returned.
   */
  DataView * getView( const std::string& name );

  /*!
   * \brief Return pointer to const view with given name or path.
   *
   * Thie method requires that all groups in the path exist if a path is given.
   *
   * If no such view exists, ATK_NULLPTR is returned.
   */
  const DataView * getView( const std::string& name ) const;

  /*!
   * \brief Return pointer to non-const view with given index.
   *
   * If no such view exists, ATK_NULLPTR is returned.
   */
  DataView * getView( IndexType idx )
  {
    SLIC_CHECK_MSG( hasView(idx),
                    "Group " << this->getName() 
                    << " has no view with index " << idx);

    return m_view_coll.getItem(idx);
  }

  /*!
   * \brief Return pointer to const view with given index.
   *
   * If no such view exists, ATK_NULLPTR is returned.
   */
  const DataView * getView( IndexType idx ) const
  {
    SLIC_CHECK_MSG( hasView(idx),
                    "Group " << this->getName() 
                    << " has no view with index " << idx);

    return m_view_coll.getItem(idx);
  }

  /*!
   * \brief Return first valid view index index in group object
   *        (i.e., smallest index over all views).
   *
   * sidre::InvalidIndex is returned if group has no views.
   */
  IndexType getFirstValidViewIndex() const
  {
    return m_view_coll.getFirstValidIndex();
  }

  /*!
   * \brief Return next valid view index in group object after given index 
   *        (i.e., smallest index over all view indices larger than given one).
   *
   * sidre::InvalidIndex is returned if there is no valid index greater
   * than given one.
   */
  IndexType getNextValidViewIndex(IndexType idx) const
  {
    return m_view_coll.getNextValidIndex(idx);
  }

//@}


//@{
//!  @name View creation methods.

  /*!
   * \brief Create and attach an undescribed (i.e., empty) DataView object
   *  with given name or path.
   *
   * IMPORTANT: To do anything useful with the view, it has to be described
   * and associated with data; for example, attach it to a data buffer and
   * apply a data description, describe it and allocate it, etc.
   *
   * If name is an empty string or group already has a view with given
   * name method does nothing.
   *
   * \return pointer to created DataView object or ATK_NULLPTR if new
   * view is not created.
   */
  DataView * createView( const std::string& name );

  /*!
   * \brief Create DataView object with given name and described by data type
   *        and number of elements, and attach new view to this group object.
   *
   * IMPORTANT: This method does not allocate data or associate the view
   * with data.
   *
   * If name is an empty string, or group already has a view with given
   * name, or given number of elements is < 0 method does nothing.
   *
   * \return pointer to created DataView object or ATK_NULLPTR if new
   * view is not created.
   */
  DataView * createView( const std::string& name,
                         TypeID type,
                         SidreLength num_elems );

  /*!
   * \brief createView(name, type, num_elems)->attachBuffer(buff)
   *
   * \return pointer to created DataView object or ATK_NULLPTR if new
   * view is not created.
   */
  DataView * createView( const std::string& name,
                         TypeID type,
                         SidreLength num_elems,
                         DataBuffer * buff )
  {
    DataView * view = createView(name, type, num_elems);
    if (view != ATK_NULLPTR)
    {
      view->attachBuffer(buff);
    }
    return view;
  }

  /*!
   * \brief createView(name, type, num_elems)->setExternalDataPtr(external_ptr)
   *
   * \return pointer to created DataView object or ATK_NULLPTR if new
   * view is not created.
   */
  DataView * createView( const std::string& name,
                         TypeID type,
                         SidreLength num_elems,
                         void * external_ptr )
  {
    DataView * view = createView(name, type, num_elems);
    if (view != ATK_NULLPTR)
    {
      view->setExternalDataPtr(external_ptr);
    }
    return view;
  }

  /*!
   * \brief Create DataView object with given name and described by data type,
   *        number of dimensions and number of elements per dimension,
   *        and attach new view to this group object.
   *
   * IMPORTANT: This method does not allocate data or associate the view
   * with data.
   *
   * If name is an empty string, or group already has a view with given
   * name, or given number of dimensions is < 0, or total number of elements
   * is < 0 method does nothing.
   *
   * \return pointer to created DataView object or ATK_NULLPTR if new
   * view is not created.
   */
  DataView * createView( const std::string& name,
                         TypeID type,
                         int ndims,
                         SidreLength * shape );

  /*!
   * \brief createView(name, type, ndims, shape)->attachBuffer(buff)
   *
   * \return pointer to created DataView object or ATK_NULLPTR if new
   * view is not created.
   */
  DataView * createView( const std::string& name,
                         TypeID type,
                         int ndims,
                         SidreLength * shape,
                         DataBuffer * buff )
  {
    DataView * view = createView(name, type, ndims, shape);
    if (view != ATK_NULLPTR)
    {
      view->attachBuffer(buff);
    }
    return view;
  }

  /*!
   * \brief createView(name, type, ndims, shape)->setExternalDataPtr(external_ptr)
   *
   * \return pointer to created DataView object or ATK_NULLPTR if new
   * view is not created.
   */
  DataView * createView( const std::string& name,
                         TypeID type,
                         int ndims,
                         SidreLength * shape,
                         void * external_ptr )
  {
    DataView * view = createView(name, type, ndims, shape);
    if (view != ATK_NULLPTR)
    {
      view->setExternalDataPtr(external_ptr);
    }
    return view;
  }

  /*!
   * \brief Create DataView object with given name and described by
   *        Conduit DataType, and attach new view to this group object.
   *
   * IMPORTANT: This method does not allocate data or associated the view
   * with data.
   *
   * If name is an empty string, or group already has a view with given
   * name, method does nothing.
   *
   * \return pointer to created DataView object or ATK_NULLPTR if new
   * view is not created.
   */
  DataView * createView( const std::string& name,
                         const DataType& dtype);

  /*!
   * \brief Create DataView object with given name, attach it to given buffer,
   *        and attach new view to this group object.
   *
   * This is equivalent to calling: createView(name)->attachBuffer(buff);
   *
   * IMPORTANT: The view cannot be used to access data in buffer until it
   * is described by calling a DataView::apply() method.
   *
   * If name is an empty string, or group already has a view with given
   * name, or given buffer pointer is null, method does nothing.
   *
   * \return pointer to created DataView object or ATK_NULLPTR if new
   * view is not created.
   */
  DataView * createView( const std::string& name,
                         DataBuffer * buff );

  /*!
   * \brief Create DataView object with given name to hold external data
   *        and attach new view to this group object.
   *
   * This is equivalent to calling:
   * createView(name)->setExternalDataPtr(external_ptr);
   *
   * IMPORTANT: Note that the view is "opaque" (it has no knowledge of
   * the type or structure of the data) until a DataView::apply() method
   * is called.
   *
   * If name is an empty string, or group already has a view with given
   * name, or given data pointer is null, method does nothing.
     //
     // RDH -- If a null data ptr is passed, should this be the same as creating
     //        an empty view?
     //
   *
   * \return pointer to created DataView object or ATK_NULLPTR if new
   * view is not created.
   */
  DataView * createView( const std::string& name,
                         void * external_ptr );

//@}


//@{
//!  @name Methods that create Views and allocate their data.

  /*!
   * \brief Create DataView object with given name, data type, and
   *        number of elements, allocate the data, and attach new
   *        view to this group object.
   *
   * This is equivalent to calling: createView(name)->allocate(type, num_elems);
   *
   * If name is an empty string, or group already has a view with given
   * name, or given number of elements is < 0 method does nothing.
   *
   * \return pointer to created DataView object or ATK_NULLPTR if new
   * view is not created.
   */
  DataView * createViewAndAllocate( const std::string& name,
                                    TypeID type,
                                    SidreLength num_elems );

  /*!
   * \brief Create DataView object with given name, data type,
   *        number of dimensions, and number of elements per dimension,
   *        allocate the data, and attach new view to this group object.
   *
   * This is equivalent to calling: createView(name)->allocate(type, ndims, shape);
   *
   * If name is an empty string, or group already has a view with given
   * name, or given number of elements is < 0 method does nothing.
   *
   * \return pointer to created DataView object or ATK_NULLPTR if new
   * view is not created.
   */
  DataView * createViewAndAllocate( const std::string& name,
                                    TypeID type,
                                    int ndims,
                                    SidreLength * shape );

  /*!
   * \brief Create DataView object with given name and Conduit DataType,
   *        allocate the data, and attach new view to this group object.
   *
   * This is equivalent to calling: createView(name)->allocate(dtype);
   *
   * If name is an empty string, or group already has a view with given
   * name, method does nothing.
   *
   * \return pointer to created DataView object or ATK_NULLPTR if new
   * view is not created.
   */
  DataView * createViewAndAllocate( const std::string& name,
                                    const DataType& dtype);

  /*!
   * \brief Create DataView object with given name for a scalar value.
   *
   * This is equivalent to calling: createView(name)->setScalar(value);
   *
   * If name is an empty string, or group already has a view with given
   * name, method does nothing.
   *
   * \return pointer to created DataView object or ATK_NULLPTR if new
   * view is not created.
   */
  template<typename ScalarType>
  DataView * createViewScalar( const std::string& name, ScalarType value)
  {
    DataView * view = createView(name);
    if (view != ATK_NULLPTR)
    {
      view->setScalar(value);
    }

    return view;
  }

  /*!
   * \brief Create DataView object with given name for a string value.
   *
   * This is equivalent to calling: createView(name)->setString(value);
   *
   * If name is an empty string, or group already has a view with given
   * name, method does nothing.
   *
   * \return pointer to created DataView object or ATK_NULLPTR if new
   * view is not created.
   */
  DataView * createViewString( const std::string& name,
                               const std::string& value);

//@}




//@{
//!  @name DataView destruction methods.

  /*!
   * \brief Destroy view in this DataGroup and leave its associated
   *        data intact.
   */
  void destroyView(DataView * view);

  /*!
   * \brief Destroy view in this DataGroup with given name and leave its
   *        associated data intact.
   */
  void destroyView(const std::string& name);

  /*!
   * \brief Destroy view in this DataGroup with given index and leave its
   *        associated data intact.
   */
  void destroyView(IndexType idx);

  /*!
   * \brief Destroy all views in this DataGroup and leave all associated
   *        data intact.
   */
  void destroyViews();

  /*!
   * \brief Destroy view in this DataGroup.  Destroy it's data also,
   *        if this is the only view referencing that data.
   *
   *        Data will not be destroyed as long as a view still exists that
   *        references it.
   */
  void destroyViewAndData( DataView * view );

  /*!
   * \brief Destroy view in this DataGroup with given name.  Destroy it's data
   *        also, if this is the only view referencing that data.
   *
   *        Data will not be destroyed as long as a view still exists that
   *        references it.
   */
  void destroyViewAndData(const std::string& name);

  /*!
   * \brief Destroy view in this DataGroup with given index.  Destroy it's data
   *        also, if this is the only view referencing that data.
   *
   *        Data will not be destroyed as long as a view still exists that
   *        references it.
   */
  void destroyViewAndData(IndexType idx);

  /*!
   * \brief Calls destroyViewAndData on all views in this group.
   */
  void destroyViewsAndData();

//@}


//@{
//!  @name View move and copy methods.

  /*!
   * \brief Remove given view object from its owning group and move it
   *        to this group.
   *
   * If given view pointer is null or group already has a view with
   * same name as given view, method is a no-op.
   *
   * \return pointer to given argument view object or ATK_NULLPTR if view
   * is not moved into this group.
   */
  DataView * moveView(DataView * view);

  /*!
   * \brief Create a copy of given view object and add it to this group.
   *
   * Note that view copying is a "shallow" copy; the data associated with
   * the view is not copied. The new view object is associated with
   * the same data as the original.
   *
   * If given group pointer is null or group already has a child group with
   * same name as given group, method is a no-op.
   *
   * \return pointer to given argument group object or ATK_NULLPTR if group
   * is not moved into this group.
   */
  DataView * copyView(DataView * view);

//@}


//@{
//!  @name Child Group query methods.

  /*!
   * \brief Return true if group has an immediate child group with given 
   * name; else false.
   */
  bool hasGroup( const std::string& name ) const
  {
    return m_group_coll.hasItem(name);
  }

  /*!
   * \brief Return true if group has an immediate child group
   *        with given index; else false.
   */
  bool hasGroup( IndexType idx ) const
  {
    return m_group_coll.hasItem(idx);
  }

  /*!
   * \brief Return the index of immediate child group with given name.
   *
   *        If no such child group exists, return sidre::InvalidIndex;
   */
  IndexType getGroupIndex(const std::string& name) const
  {
    SLIC_CHECK_MSG(hasGroup(name),
                   "Group " << this->getName() << 
                   " has no child group with name '" << name << "'");

    return m_group_coll.getItemIndex(name);
  }

  /*!
   * \brief Return the name of immediate child group with given index.
   *
   *        If no such child group exists, return sidre::InvalidName.
   */
  const std::string& getGroupName(IndexType idx) const
  {
    SLIC_CHECK_MSG(hasGroup(idx),
                   "Group " << this->getName() << 
                   " has no child group with index " << idx);

    return m_group_coll.getItemName(idx);
  }

//@}


//@{
//!  @name Group access and iteration methods.

  /*!
   * \brief Return pointer to non-const child group with given name or path.
   *
   * Thie method requires that all groups in the path exist if a path is given.
   *
   * If no such group exists, ATK_NULLPTR is returned.
   */
  DataGroup * getGroup( const std::string& name );

  /*!
   * \brief Return pointer to const child group with given name or path.
   *
   * Thie method requires that all groups in the path exist if a path is given.
   *
   * If no such group exists, ATK_NULLPTR is returned.
   */
  DataGroup const * getGroup( const std::string& name ) const;

  /*!
   * \brief Return pointer to non-const immediate child group with given index.
   * 
   * If no such group exists, ATK_NULLPTR is returned.
   */
  DataGroup * getGroup( IndexType idx )
  {
    SLIC_CHECK_MSG(hasGroup(idx),
                   "Group " << this->getName() << 
                   " has no child group with index " << idx);

    return m_group_coll.getItem(idx);
  }

  /*!
   * \brief Return pointer to const immediate child group with given index.
   * 
   * If no such group exists, ATK_NULLPTR is returned.
   */
  const DataGroup * getGroup( IndexType idx ) const
  {
    SLIC_CHECK_MSG(hasGroup(idx),
                   "Group " << this->getName() << 
                   " has no child group with index " << idx);

    return m_group_coll.getItem(idx);
  }

  /*!
   * \brief Return first valid child group index (i.e., smallest
   *        index over all child groups).
   *
   * sidre::InvalidIndex is returned if group has no child groups.
   */
  IndexType getFirstValidGroupIndex() const
  {
    return m_group_coll.getFirstValidIndex();
  }

  /*!
   * \brief Return next valid child group index after given index
   *        (i.e., smallest index over all child group indices larger
   *        than given one).
   *
   * sidre::InvalidIndex is returned if there is no valid index greater
   * than given one.
   */
  IndexType getNextValidGroupIndex(IndexType idx) const
  {
    return m_group_coll.getNextValidIndex(idx);
  }

//@}


//@{
//!  @name Child Group creation and destruction methods.

  /*!
   * \brief Create a child group within this group with given name or path.
   *
   * If name is an empty string or group already has a child group with
   * given name or path, method is a no-op.
   *
   * \return pointer to created DataGroup object or ATK_NULLPTR if new
   * group is not created.
   */
  DataGroup * createGroup( const std::string& name );

  /*!
   * \brief Destroy child group in this group with given name or path.
   * 
   * If no such group exists, method is a no-op.
   */
  void destroyGroup(const std::string& name);

  /*!
   * \brief Destroy child group within this group with given index.
   * 
   * If no such group exists, method is a no-op.
   */
  void destroyGroup(IndexType idx);

  /*!
   * \brief Destroy all child groups in this group.
   *
   * Note that this will recrusively destroy entire group sub-tree below 
   * this group.
   */
  void destroyGroups();

//@}


//@{
//!  @name Group move and copy methods.

  /*!
   * \brief Remove given group object from its parent group and make it 
   *        a child of this group.
   *
   * If given group pointer is null or group already has a child group with
   * same name as given group, method is a no-op.
   *
   * \return pointer to given argument group object or ATK_NULLPTR if group
   * is not moved into this group.
   */
  DataGroup * moveGroup(DataGroup * group);

  /*!
   * \brief Create a copy of group hierarchy rooted at given group and make it 
   *        a child of this group.
   *
   * Note that all views in the group hierarchy are copied as well.
   *
   * Note that group copying is a "shallow" copy; the data associated
   * with views in a group are not copied. In particular, the new group
   * hierachy and all its views is associated with the same data as the 
   * given group.
   *
   * If given group pointer is null or group already has a child group with
   * same name as given group, method is a no-op.
   *
   * \return pointer to given argument group object or ATK_NULLPTR if group
   * is not moved into this group.
   */
  DataGroup * copyGroup(DataGroup * group);

//@}


//@{
//!  @name Group print methods.

  /*!
   * \brief Print JSON description of data group to stdout.
   *
   * Note that this will recursively print entire group (sub) tree
   * starting at this group object.
   */
  void print() const;

  /*!
   * \brief Print JSON description of data group to an ostream.
   *
   * Note that this will recursively print entire group (sub) tree
   * starting at this group object.
   */
  void print(std::ostream& os) const;


  /*!
   * \brief Print given number of levels of group (sub) tree
   *        starting at this group object to an output stream.
   */
  void printTree( const int nlevels, std::ostream& os ) const;

//@}

  /*!
   * \brief Copy description of group hierarchy rooted at this group to 
   * given Conduit node.
   *
   * The description includes views of this group and all of its children
   * recursively.
   */
  void info(Node& n) const;

  /*!
   * \brief Return true if this group is equivalent to given group; else false. 
   *
   * Two groups are equivalent if they are the root groups of identical 
   * group hierarchy structures with the same names for all views and 
   * groups in the hierarchy, and the views are also equivalent.
   *
   * \sa DataView::isEquivalentTo
   */
  bool isEquivalentTo(const DataGroup * other) const;

private:

//@{
//!  @name Private DataGroup ctors and dtors
//!        (callable only by DataStore and DataGroup methods).

  /*!
   *  \brief Private ctor that creates a Group with given name
   *         in given parent group.
   */
  DataGroup(const std::string& name, DataGroup * parent);

  /*!
   *  \brief Private ctor that creates a Group with given name
   *         in the given DataStore root group.
   */
  DataGroup(const std::string& name, DataStore * datastore);

  //
  // Unimplemented copy ctors and copy-assignment operators.
  //
#ifdef USE_CXX11
  DataGroup( const DataGroup& source ) = delete;
  DataGroup( DataGroup&& source ) = delete;

  DataGroup& operator=( const DataGroup& rhs ) = delete;
  DataGroup& operator=( const DataGroup&& rhs ) = delete;
#else
  DataGroup( const DataGroup& source );
  DataGroup& operator=( const DataGroup& rhs );
#endif

  /*!
   * \brief Destructor destroys all views and child groups.
   */
  ~DataGroup();

//@}


//@{
//!  @name Private DataGroup view and buffer manipulation methods.

  /*!
   * \brief Private methods to attach/detach DataView object to DataGroup.
   */
  DataView * attachView(DataView * view);
  ///
  DataView * detachView(const DataView * view)
  {
    return detachView(view->getName());
  }
  //
  DataView * detachView(const std::string& name);
  ///
  DataView * detachView(IndexType idx);

//@}


//@{
//!  @name DataGroup (child) group manipulation methods.

  /*!
   * \brief Private methods to attach/detach DataGroup object to DataGroup.
   */
  DataGroup * attachGroup(DataGroup * group);
  ///
  DataGroup * detachGroup(const std::string& name);
  ///
  DataGroup * detachGroup(IndexType idx);

//@}

  /*!
   * \brief Private method to retrieve the next-to-last entry in a path.  This
   * entry is usually the one that needs to perform an action, such as creating or
   * retrieving a group or view.
   *
   * path - The path to traverse.  This parameter is modified during algorithm
   * execution.  Upon completion it will contain the last entry in the path.  This
   * is typically a name of a group or view that needs to be created or retrieved.
   *
   * create_on_demand - This controls whether any missing groups should be created
   * while traversing a path.
   */
  DataGroup * walkPath(std::string& path, bool create_on_demand );

//@{
//!  @name Private DataGroup methods for interacting with Conduit Nodes.

  /*!
   * \brief Private methods to copy DataGroup to Conduit Node.
   *
   * \param buffer_indices Used to track what buffers are referenced
   * by the views in this group and sub-groups.
   */
  void exportTo(conduit::Node& data_holder,
                std::set<IndexType>& buffer_indicies) const;

  /*!
   * \brief Private methods to copy DataGroup from Conduit Node.
   *
   * Map of buffer indices tracks old buffer ids in the file to the
   * new buffer ids in the datastore.  Buffer ids are not guaranteed
   * to remain the same when a tree is restored.
   *
   */
  void importFrom(conduit::Node& node,
                  const std::map<IndexType, IndexType>& buffer_id_map);

//@}


  /// Name of this DataGroup object.
  std::string m_name;

  /// Parent DataGroup of this DataGroup object.
  DataGroup * m_parent;

  /// This DataGroup object lives in the tree of this DataStore object.
  DataStore * m_datastore;

  /// Character used to denote a path string passed to get/create calls.
  static const char s_path_delimiter;

  ///
  /// Typedefs for view and shild group containers. They are here to
  /// avoid propagating specific type names in the DataGroup class
  /// implementation when we experiment with different containers.
  ///
  ///////////////////////////////////////////////////////////////////
  //
  // Associative container options
  //
  // To try a different container, set the "MapType" typedef to
  // what you want.  Note: only one typedef should be active!!!
  //
  // Current options are std::map and boost/std::unordered_map
  ///
  // typedef std::map<std::string, IndexType> MapType;
  ///
#if defined(USE_UNORDERED_MAP)
#if defined(USE_CXX11)
  typedef std::unordered_map<std::string, IndexType> MapType;
#else
  typedef boost::unordered_map<std::string, IndexType> MapType;
#endif
#else
#if defined(USE_DENSE_HASH_MAP)
  typedef google::dense_hash_map<std::string, IndexType> MapType;
#endif
#endif
  //
  ///////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////
  // Collection options (API between DataGroup and assoc. container)
  ///////////////////////////////////////////////////////////////////
  //
#if defined(USE_NEW_MAP_COLLECTION)
  ///////////////////////////////////////////////////////////////////
  // Improved implementation (index-item association constant as long
  // as item is in collection, but holes in index sequence)
  ///////////////////////////////////////////////////////////////////
  //
  typedef NewMapCollection<DataView, MapType> DataViewCollection;
  //
  typedef NewMapCollection<DataGroup, MapType> DataGroupCollection;
  ///////////////////////////////////////////////////////////////////
#else
  ///////////////////////////////////////////////////////////////////
  // Original implementation (no holes in index sequence)
  ///////////////////////////////////////////////////////////////////
  //
  typedef MapCollection<DataView, MapType> DataViewCollection;
  //
  typedef MapCollection<DataGroup, MapType> DataGroupCollection;
  //
  ///////////////////////////////////////////////////////////////////
  //
#endif

  /// Collection of DataViews
  DataViewCollection m_view_coll;

  /// Collection of child DataGroups
  DataGroupCollection m_group_coll;

};


} /* end namespace sidre */
} /* end namespace asctoolkit */

#endif /* DATAGROUP_HPP_ */

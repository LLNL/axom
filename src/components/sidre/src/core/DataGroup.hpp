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
//!  @name Methods to create a View that has no associated data.
//!
//! IMPORTANT: These methods do not allocate data or associate a view
//! with data. Thus, to do anything useful with a view created by one
//! of these methods, the view should be allocated, attached to a buffer
//! or attached to externally-owned data.
//! 
//! Each of these methods is a no-op if the given view name is an
//! empty string or the group already has a view with given name or path.
//!
//! Additional conditions under which a method can be a no-op are described
//! for each method.

  /*!
   * \brief Create an undescribed (i.e., empty) View object with given name 
   * or path in this group.
   *
   * \return pointer to new View object or ATK_NULLPTR if one is not created.
   */
  DataView * createView( const std::string& name );

  /*!
   * \brief Create View object with given name or path in this group that
   *  has a data description with data type and number of elements.
   *
   * If given data type is undefined, or given number of elements is < 0,
   * method is a no-op.
   *
   * \return pointer to new View object or ATK_NULLPTR if one is not created.
   */
  DataView * createView( const std::string& name,
                         TypeID type,
                         SidreLength num_elems );

  /*!
   * \brief Create View object with given name or path in this group that
   *  has a data description with data type and shape.
   *
   * If given data type is undefined, or given number of dimensions is < 0,
   * or given shape ptr is null, method is a no-op.
   *
   * \return pointer to new View object or ATK_NULLPTR if one is not created.
   */
  DataView * createView( const std::string& name,
                         TypeID type,
                         int ndims,
                         SidreLength * shape );

  /*!
   * \brief Create View object with given name or path in this group that
   *  is described by a Conduit DataType object.
   *
   * \return pointer to new View object or ATK_NULLPTR if one is not created.
   */
  DataView * createView( const std::string& name,
                         const DataType& dtype);

//@}


//@{
//!  @name Methods to create a View with a data buffer attached.
//!
//! IMPORTANT: The buffer passed to each of these methods may or may not 
//! be allocated. Thus, to do anything useful with a view created by one
//! of these methods, the buffer must be allocated and it must be compatible 
//! with the view data description.
//! 
//! Each of these methods is a no-op if the given view name is an
//! empty string or the group already has a view with given name or path.
//! 
//! Also, calling one of these methods with a null buffer pointer is 
//! similar to creating a view with no data association.
//!
//! Additional conditions under which a method can be a no-op are described
//! for each method.

  /*!
   * \brief Create an undescribed View object with given name or path in 
   * this group and attach given buffer to it.
   *
   * IMPORTANT: The view cannot be used to access data in buffer until it
   * is described by calling a DataView::apply() method.
   *
   * This method is equivalent to: 
   * group->createView(name)->attachBuffer(buff).
   *
   * \return pointer to new View object or ATK_NULLPTR if one is not created.
   *
   * \sa DataView::attachBuffer
   */
  DataView * createView( const std::string& name,
                         DataBuffer * buff );

  /*!
   * \brief Create View object with given name or path in this group that
   * has a data description with data type and number of elements and
   * attach given buffer to it.
   *
   * If given data type is undefined, or given number of elements is < 0,
   * method is a no-op.
   *
   * This method is equivalent to: 
   * group->createView(name, type, num_elems)->attachBuffer(buff), or
   * group->createView(name)->attachBuffer(buff)->apply(type, num_elems).
   *
   * \return pointer to new View object or ATK_NULLPTR if one is not created.
   *
   * \sa DataView::attachBuffer
   */
  DataView * createView( const std::string& name,
                         TypeID type,
                         SidreLength num_elems,
                         DataBuffer * buff );

  /*!
   * \brief Create View object with given name or path in this group that
   * has a data description with data type and shape and attach given 
   * buffer to it.
   *
   * If given data type is undefined, or given number of dimensions is < 0,
   * or given shape ptr is null, method is a no-op.
   *
   * This method is equivalent to: 
   * group->createView(name, type, ndims, shape)->attachBuffer(buff), or
   * group->createView(name)->attachBuffer(buff)->apply(type, ndims, shape).
   *
   * \return pointer to new View object or ATK_NULLPTR if one is not created.
   *
   * \sa DataView::attachBuffer
   */
  DataView * createView( const std::string& name,
                         TypeID type,
                         int ndims,
                         SidreLength * shape,
                         DataBuffer * buff );

  /*!
   * \brief Create View object with given name or path in this group that
   *  is described by a Conduit DataType object and attach given buffer to it.
   *
   * This method is equivalent to: 
   * group->createView(name, dtype)->attachBuffer(buff), or
   * group->createView(name)->attachBuffer(buff)->apply(dtype).
   *
   * \return pointer to new View object or ATK_NULLPTR if one is not created.
   *
   * \sa DataView::attachBuffer
   */
  DataView * createView( const std::string& name,
                         const DataType& dtype, 
                         DataBuffer * buff );

//@}


//@{
//!  @name Methods to create a View with externally-owned data attached.
//!
//! IMPORTANT: To do anything useful with a view created by one of these 
//! methods, the external data must be allocated and compatible with the
//! view description. 
//! 
//! Each of these methods is a no-op if the given view name is an
//! empty string or the group already has a view with given name or path.
//!
//! Additional conditions under which a method can be a no-op are described
//! for each method.

  /*!
   * \brief Create View object with given name with given name or path in
   * this group and attach external data ptr to it.
   *
   * IMPORTANT: Note that the view is "opaque" (it has no knowledge of
   * the type or structure of the data) until a DataView::apply() method
   * is called.
   *
   * This method is equivalent to: 
   * group->createView(name)->setExternalDataPtr(external_ptr).
   *
   * \return pointer to new View object or ATK_NULLPTR if one is not created.
   *
   * \sa DataView::setExternalDataPtr
   */
  DataView * createView( const std::string& name,
                         void * external_ptr );

  /*!
   * Create View object with given name or path in this group that
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
   * \return pointer to new View object or ATK_NULLPTR if one is not created.
   *
   * \sa DataView::setExternalDataPtr
   */
  DataView * createView( const std::string& name,
                         TypeID type,
                         SidreLength num_elems,
                         void * external_ptr );


  /*!
   * \brief Create View object with given name or path in this group that
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
   * \return pointer to new View object or ATK_NULLPTR if one is not created.
   *
   * \sa DataView::setExternalDataPtr
   */
  DataView * createView( const std::string& name,
                         TypeID type,
                         int ndims,
                         SidreLength * shape,
                         void * external_ptr );
  /*!
   * \brief Create View object with given name or path in this group that
   * is described by a Conduit DataType object and attach externally-owned 
   * data to it.
   *
   * This method is equivalent to:
   * group->createView(name, dtype)->setExternalDataPtr(external_ptr), or
   * group->createView(name)->setExternalDataPtr(external_ptr)->apply(dtype).
   *
   * \return pointer to new View object or ATK_NULLPTR if one is not created.
   *
   * \sa DataView::attachBuffer
   */
  DataView * createView( const std::string& name,
                         const DataType& dtype,
                         void * external_ptr );

//@}


//@{
//!  @name Methods to create a View and allocate data for it.
//! 
//! Each of these methods is a no-op if the given view name is an
//! empty string or the group already has a view with given name or path.
//!
//! Additional conditions under which a method can be a no-op are described
//! for each method.

  /*!
   * \brief Create View object with given name or path in this group that
   * has a data description with data type and number of elements and
   * allocate data for it.
   *
   * If given data type is undefined, or given number of elements is < 0,
   * method is a no-op.
   * 
   * This is equivalent to: createView(name)->allocate(type, num_elems), or
   * createView(name, type, num_elems)->allocate()
   *
   * \return pointer to new View object or ATK_NULLPTR if one is not created.
   *
   * \sa DataView::allocate
   */
  DataView * createViewAndAllocate( const std::string& name,
                                    TypeID type,
                                    SidreLength num_elems );

  /*!
   * \brief Create View object with given name or path in this group that
   * has a data description with data type and shape and allocate data for it.
   *
   * If given data type is undefined, or given number of dimensions is < 0,
   * or given shape ptr is null, method is a no-op.
   *
   * This method is equivalent to:
   * group->createView(name)->allocate(type, ndims, shape), or
   * createView(name, type, ndims, shape)->allocate().
   *
   * \return pointer to new View object or ATK_NULLPTR if one is not created.
   *
   * \sa DataView::allocate
   */
  DataView * createViewAndAllocate( const std::string& name,
                                    TypeID type,
                                    int ndims,
                                    SidreLength * shape );

  /*!
   * \brief Create View object with given name or path in this group that
   * is described by a Conduit DataType object and allocate data for it.
   *
   * This method is equivalent to:
   * group->createView(name)->allocate(dtype), or 
   * group->createView(name, dtype)->allocate().
   *
   * If given data type object is empty, data will not be allocated.
   *
   * \return pointer to new View object or ATK_NULLPTR if one is not created.
   *
   * \sa DataView::allocate
   */
  DataView * createViewAndAllocate( const std::string& name,
                                    const DataType& dtype);

  /*!
   * \brief Create View object with given name or path in this group 
   * set its data to given scalar value. 
   *
   * This is equivalent to: createView(name)->setScalar(value);
   *
   * If given data type object is empty, data will not be allocated.
   *
   * \return pointer to new View object or ATK_NULLPTR if one is not created.
   *
   * \sa DataView::setScalar
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
   * \brief Create View object with given name or path in this group
   * set its data to given string.
   *
   * This is equivalent to: createView(name)->setScalar(value);
   *
   * If given data type object is empty, data will not be allocated.
   *
   * \return pointer to new View object or ATK_NULLPTR if one is not created.
   *
   * \sa DataView::setString
   */
  DataView * createViewString( const std::string& name,
                               const std::string& value);

//@}


//@{
//!  @name View destruction methods.
//!
//! Each of these methods is a no-op if the specified view does not exist.

  /*!
   * \brief Destroy view with given name owned by this group, but leave
   * its data intect.
   */
  void destroyView(const std::string& name);

  /*!
   * \brief Destroy view with given index owned by this group, but leave
   * its data intect.
   */
  void destroyView(IndexType idx);

  /*!
   * \brief Destroy all views owned by this group, but leave all their
   *        data intact.
   */
  void destroyViews();

  /*!
   * \brief Destroy view with given name owned by this group and deallocate
   * its data if it's the only view associated with that data.
   */
  void destroyViewAndData(const std::string& name);

  /*!
   * \brief Destroy view with given index owned by this group and deallocate
   * its data if it's the only view associated with that data.
   */
  void destroyViewAndData(IndexType idx);

  /*!
   * \brief Destroy all views owned by this group and deallocate
   * data for each view when it's the only view associated with that data.
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

  /*!
   *  Unimplemented ctors and copy-assignment operators.
   */
#ifdef USE_CXX11
  DataGroup( const DataGroup& source ) = delete;
  DataGroup( DataGroup&& source ) = delete;

  DataGroup& operator=( const DataGroup& rhs ) = delete;
  DataGroup& operator=( const DataGroup&& rhs ) = delete;
#else
  DataGroup( const DataGroup& source );
  DataGroup& operator=( const DataGroup& rhs );
#endif

//@{
//!  @name Private Group ctors and dtors
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

  /*!
   * \brief Destructor destroys all views and child groups.
   */
  ~DataGroup();

//@}


//@{
//!  @name Private Group view manipulation methods.

  /*!
   * \brief Attach View object to this Group.
   */
  DataView * attachView(DataView * view);

  /*!
   * \brief Detach View object from this Group.
   */
  DataView * detachView(const DataView * view)
  {
    return detachView(view->getName());
  }

  /*!
   * \brief Detach View with given name from this Group.
   */
  DataView * detachView(const std::string& name);

  /*!
   * \brief Detach View with given index from this Group.
   */
  DataView * detachView(IndexType idx);

  /*!
   * \brief Destroy view and its data if its data is not shared with any
   * other view.
   *
   * Data will not be destroyed as long as a view still exists that
   * references it.
   *
   * IMPORTANT: this method assumes view is owned by this group.
   */
  void destroyViewAndData( DataView * view );

//@}


//@{
//!  @name Private (child) Group manipulation methods.

  /*!
   * \brief Attach Group to this Group as a child.
   */
  DataGroup * attachGroup(DataGroup * group);

  /*!
   * \brief Detaich child Group with given name from this Group.
   */
  DataGroup * detachGroup(const std::string& name);

  /*!
   * \brief Detaich child Group with given index from this Group.
   */
  DataGroup * detachGroup(IndexType idx);

//@}


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

  /*!
   * \brief Private method that returns the Group that is the next-to-last 
   * entry in a slash-deliminated ("/") path string.
   *
   * The string before the last "/" character, if there is one, is the 
   * next-to-last path entry. In this case, the return value is that group
   * in the path. 
   *
   * If there is no "/" in the given path, the entire string is considered 
   * the next-to-last path entry. In this case, the erturn value is this
   * group.
   *
   * The path argument is modified while walking the path. Its value when 
   * the method returns is the last entry in the path, either the string
   * following the last "/" in the input (if there is one) or the entire
   * input path string if it contains no "/".
   */
  DataGroup * walkPath(std::string& path, bool create_groups_in_path );


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

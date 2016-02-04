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
   * \brief Return const reference to name of DataGroup instance.
   */
  const std::string& getName() const
  {
    return m_name;
  }

  /*!
   * \brief Return pointer to non-const DataGroup parent of group.
   */
  DataGroup * getParent()
  {
    return m_parent;
  }

  /*!
   * \brief Return pointer to const DataGroup parent of group.
   */
  const DataGroup * getParent() const
  {
    return m_parent;
  }

  /*!
   * \brief Return pointer to non-const DataStore that owns group hierarchy
   *        to which DataGroup instance belongs.
   */
  DataStore * getDataStore()
  {
    return m_datastore;
  }

  /*!
   * \brief Return pointer to const DataStore that owns group hierarchy
   *        to which DataGroup instance belongs.
   */
  const DataStore * getDataStore() const
  {
    return m_datastore;
  }

  /*!
   * \brief Return number of DataViews contained in this DataGroup.
   */
  size_t getNumViews() const
  {
    return m_view_coll.getNumItems();
  }

  /*!
   * \brief Return number of (child) DataGroups contained in this DataGroup.
   */
  size_t getNumGroups() const
  {
    return m_group_coll.getNumItems();
  }

//@}


//@{
//!  @name DataView query methods.

  /*!
   * \brief Return true if DataGroup owns a DataView with given name;
   *        else false.
   */
  bool hasView( const std::string& name ) const
  {
    return m_view_coll.hasItem(name);
  }

  /*!
   * \brief Return true if DataGroup owns a DataView with given index;
   *        else false.
   */
  bool hasView( IndexType idx ) const
  {
    return m_view_coll.hasItem(idx);
  }

  /*!
   * \brief Return the index of DataView with given name.
   *
   *        If none, return sidre::InvalidIndex;
   */
  IndexType getViewIndex(const std::string& name) const
  {
    SLIC_CHECK_MSG( !name.empty() && hasView(name),
                    "Group " << this->getName() << " does not have view with name '" << name <<
                    "'");

    return m_view_coll.getItemIndex(name);
  }

  /*!
   * \brief Return the name of DataView with given index.
   *
   *        If none, return sidre::InvalidName.
   */
  const std::string& getViewName(IndexType idx) const
  {
    SLIC_CHECK_MSG( hasView(idx),
                    "Group " << this->getName() << " does not have view with index " <<
                    idx);

    return m_view_coll.getItemName(idx);
  }

//@}


//@{
//!  @name DataView access and iteration methods.
//!
//!  Some of these methods support a path syntax for input.  When a path is
//!  provided the method will retrieve the last item in the path.  All
//!  items in the path must exist.
//!
//!  Example:
//!  getView("foo/bar/baz")
//!  is equivalent to
//!  getGroup("foo")->getGroup("bar")->getView("baz")
//!


  /*!

   * \brief Return (non-const) pointer to DataView with given name or path.
   * This algorithm requires that all groups in the path are already created.
   * Example:
   * getView("foo/bar/baz")
   * is equivalent to
   * getGroup("foo")->getGroup("bar")->getView("baz")
   */
  DataView * getView( const std::string& name )
  {
    std::string path = name;
    DataGroup * group = walkPath( path, false );

    SLIC_CHECK_MSG( !path.empty() && group->hasView(
                      name),
                    "Group " << group->getName() << " can't retrieve view with name '" << path <<
                    "'");

    return group->m_view_coll.getItem(path);
  }

  /*!
   * \brief Return (const) pointer to DataView with given name.
   */
  // TODO - Add path support to const function version.
  const DataView * getView( const std::string& name ) const
  {
    SLIC_CHECK_MSG( !name.empty() && hasView(name),
                    "Group " << this->getName() << " can't retrieve view with name '" << name <<
                    "'");

    return m_view_coll.getItem(name);
  }

  /*!
   * \brief Return (non-const) pointer to DataView with given index.
   */
  DataView * getView( IndexType idx )
  {
    SLIC_CHECK_MSG( hasView(idx),
                    "Group " << this->getName() << " does not have view with index " <<
                    idx);

    return m_view_coll.getItem(idx);
  }

  /*!
   * \brief Return (const) pointer to DataView with given index.
   */
  // TODO - Add path support to const function version.
  const DataView * getView( IndexType idx ) const
  {
    SLIC_CHECK_MSG( hasView(idx),
                    "Group " << this->getName() << " does not have view with index " <<
                    idx);

    return m_view_coll.getItem(idx);
  }

  /*!
   * \brief Return first valid DataView index (i.e., smallest index
   *        over all DataViews).
   *
   * sidre::InvalidIndex is returned if group has no views.
   */
  IndexType getFirstValidViewIndex() const
  {
    return m_view_coll.getFirstValidIndex();
  }

  /*!
   * \brief Return next valid DataView index after given index (i.e.,
   *        smallest index over all view indices larger than given one).
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
  //!  @name Methods to create a DataView
  //!
  //!  These methods support a path syntax for input.  When a path is provided
  //!  the method will traverse down to the last group in the path before
  //!  executing the called operation.  Creation methods will automatically
  //!  create a chain of groups if they are not already present.
  //!
  //!  Example:
  //!  createView("foo/bar/baz")
  //!  is equivalent to
  //!  createGroup("foo")->createGroup("bar")->createView("baz")
  //!

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
   * \brief Create DataView object with given name and described by
   *        Conduit Schema, and attach new view to this group object.
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
                         const Schema& schema);

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
//!  @name DataView creation methods that also describe and allocate data.
//!
//!  These methods support a path syntax for input.  When a path is provided
//!  the method will traverse down to the last group in the path before
//!  executing the called operation.  Creation methods will automatically
//!  create a chain of groups if they are not already present.
//!
//!  Example:
//!  createView("foo/bar/baz")
//!  is equivalent to
//!  createGroup("foo")->createGroup("bar")->createView("baz")


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
   * \brief Create DataView object with given name and Conduit Schema,
   *        allocate the data, and attach new view to this group object.
   *
   * This is equivalent to calling: createView(name)->allocate(schema);
   *
   * If name is an empty string, or group already has a view with given
   * name, method does nothing.
   *
   * \return pointer to created DataView object or ATK_NULLPTR if new
   * view is not created.
   */
  DataView * createViewAndAllocate( const std::string& name,
                                    const Schema& schema);

//@}




//@{
//!  @name DataView destruction methods.

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
   * \brief Destroy view in this DataGroup with given name AND destroy
   *        its associated data, if it owns the data.
   */
  void destroyViewAndData(const std::string& name);

  /*!
   * \brief Destroy view in this DataGroup with given index AND destroy
   *        its associated data, if it owns the data.
   */
  void destroyViewAndData(IndexType idx);

  /*!
   * \brief Destroy all views in this DataGroup AND destroy their
   *        associated data.
   */
  void destroyViewsAndData();

//@}


//@{
//!  @name DataView move and copy methods.

  /*!
   * \brief Remove DataView object from its owning group and attach
   *        to this DataGroup object.
   *
   * If given view pointer is null or group already contains a view with
   * same name as given view, method does nothing.
   *
   * \return pointer to given DataView object or ATK_NULLPTR if new
   * view is not moved into this group.
   */
  DataView * moveView(DataView * view);

  /*!
   * \brief Create a copy of given DataView object and attach
   *        to this DataGroup object.
   *
   * Note that this is a "shallow" copy; the data associated with
   * the view is not copied. The new view object is associated with
   * the same data.
   *
   * If given view pointer is null or group already contains a view with
   * same name as given view, method does nothing.
   *
   * \return pointer to given DataView object or ATK_NULLPTR if new
   * view is not copied into this group.
   */
  DataView * copyView(DataView * view);

//@}


//@{
//!  @name (child) DataGroup query methods.

  bool hasGroup( const std::string& name ) const
  {
    return m_group_coll.hasItem(name);
  }

  /*!
   * \brief Return true if DataGroup has an (immediate) child DataGroup
   *        with given index; else false.
   */
  bool hasGroup( IndexType idx ) const
  {
    return m_group_coll.hasItem(idx);
  }

  /*!
   * \brief Return the index of child DataGroup with given name.
   *
   *        If none, return sidre::InvalidIndex;
   */
  IndexType getGroupIndex(const std::string& name) const
  {
    SLIC_CHECK_MSG( !name.empty() && hasGroup(name),
                    "Group " << this->getName() << " does not have group with name '" << name <<
                    "'");

    return m_group_coll.getItemIndex(name);
  }

  /*!
   * \brief Return the name of child DataGroup with given index.
   *
   *        If none, return sidre::InvalidName.
   */
  const std::string& getGroupName(IndexType idx) const
  {
    SLIC_CHECK_MSG( hasGroup(idx),
                    "Group " << this->getName() << " does not have group with index " <<
                    idx);

    return m_group_coll.getItemName(idx);
  }

//@}


//@{
//!  @name (child) DataGroup access and iteration methods.
//!
//!  Some of these methods support a path syntax for input.  When a path is
//!  provided the method will retrieve the last item in the path.  All
//!  items in the path must exist.
//!
//!  Example:
//!  getGroup("foo/bar/baz")
//!  is equivalent to
//!  getGroup("foo")->getGroup("bar")->getGroup("baz")
//!

  /*!
   * \brief Return (non-const) pointer to child DataGroup with given name.
   */
  DataGroup * getGroup( const std::string& name )
  {
    std::string path = name;
    DataGroup * group = walkPath( path, false );

    SLIC_CHECK_MSG( !path.empty() && group->hasGroup(
                      path),
                    "Group " << group->getName() << " does not have group with name '" << path <<
                    "'");

    return group->m_group_coll.getItem(path);
  }

  /*!
   * \brief Return (const) pointer to to child DataGroup with given name.
   */
  DataGroup const * getGroup( const std::string& name ) const
  {
    SLIC_CHECK_MSG( !name.empty() && hasGroup(name),
                    "Group " << getName() << " does not have group with name '" << name <<
                    "'");

    return m_group_coll.getItem(name);
  }

  /*!
   * \brief Return (non-const) pointer to child DataGroup with given index.
   */
  DataGroup * getGroup( IndexType idx )
  {
    SLIC_CHECK_MSG( hasGroup(idx),
                    "Group " << this->getName() << " does not have group with index " <<
                    idx);

    return m_group_coll.getItem(idx);
  }

  /*!
   * \brief Return (const) pointer to child DataGroup with given index.
   */
  const DataGroup * getGroup( IndexType idx ) const
  {
    SLIC_CHECK_MSG( hasGroup(idx),
                    "Group " << this->getName() << " does not have group with index " <<
                    idx);

    return m_group_coll.getItem(idx);
  }

  /*!
   * \brief Return first valid (child) DataGroup index (i.e., smallest
   *        index over all (child) DataGroups).
   *
   * sidre::InvalidIndex is returned if group has no child groups.
   */
  IndexType getFirstValidGroupIndex() const
  {
    return m_group_coll.getFirstValidIndex();
  }

  /*!
   * \brief Return next valid (child) DataGroup index after given index
   *        (i.e., smallest index over all chid group indices larger
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
//!  @name (child) DataGroup create and destroy methods.
//!  The creation methods support a path syntax for input.  When a path is
//!  provided the method will create all groups in the method (if not already
//!  present).
//!  Example:
//!  createGroup("foo/bar/baz")
//!  is equivalent to
//!  createGroup("foo")->createGroup("bar")->createGroup("baz")

  /*!
   * \brief Create a DataGroup object with given name and attach as a
   * child group of this group.
   *
   * If name is an empty string or group already has a child group with
   * given name, method does nothing.
   *
   * \return pointer to created DataGroup object or ATK_NULLPTR if new
   * group is not created.
   */
  DataGroup * createGroup( const std::string& name );

  /*!
   * \brief Destroy child group in this DataGroup with given name.
   */
  void destroyGroup(const std::string& name);

  /*!
   * \brief Destroy child group in this DataGroup with given index.
   */
  void destroyGroup(IndexType idx);

  /*!
   * \brief Destroy all DataGroups in this DataGroup.
   *
   * Note that this will recrusively destroy all child groups of the
   * child groups in this DataGroup.
   */
  void destroyGroups();

//@}


//@{
//!  @name (child) DataGroup move and copy methods.

  /*!
   * \brief Remove DataGroup object from its parent group and attach
   *        to this DataGroup object.
   *
   * If given group pointer is null or group already contains a group with
   * same name as given group, method does nothing.
   *
   * \return pointer to given DataGroup object or ATK_NULLPTR if group
   * is not moved into this group.
   */
  DataGroup * moveGroup(DataGroup * group);

  /*!
   * \brief Create a copy of given DataGroup object (including all of its
   *        DataViews and child DataGroups) and attach to this DataGroup
   *        object.
   *
   * Note that DataGroup copying is a "shallow" copy; the data associated
   * with views in a group are not copied. The new DataGroup is associated
   * with the same data as the given group.
   *
   * If given group pointer is null or group already contains a group with
   * same name as given group, method does nothing.
   *
   * \return pointer to new DataGroup object in this group or ATK_NULLPTR if
   * new group is not created.
   */
  DataGroup * copyGroup(DataGroup * group);

//@}


//@{
//!  @name DataGroup print methods.

  /*!
   * \brief Print JSON description of data group to stdout.
   *
   * Note that this will recursively print entire group (sub) tree
   * starting at this DataGroup object.
   */
  void print() const;

  /*!
   * \brief Print JSON description of data group to an ostream.
   *
   * Note that this will recursively print entire group (sub) tree
   * starting at this DataGroup object.
   */
  void print(std::ostream& os) const;


  /*!
   * \brief Print given number of levels of group (sub) tree
   *        starting at this DataGroup object to an output stream.
   */
  void printTree( const int nlevels, std::ostream& os ) const;

//@}


//@{
//!  @name DataGroup save and load methods (basic file I/O).

  /*!
   * \brief Save this DataGroup object (including data views and child
   *        groups) to a file set named "obase".
   *
   * \warning Currently, only valid protocol is "conduit".
   */
  void save(const std::string& obase,
            const std::string& protocol) const;

  /*!
   * \brief Load data group (including data views and child groups)
   *        from a file set named "obase" into this DataGroup object.
   *
   * \warning Currently, only valid protocol is "conduit".
   */
  void load(const std::string& obase,
            const std::string& protocol);

//@}

  /*!
   * \brief Copy data group description to given Conduit node.
   */
  void info(Node& n) const;

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
   * \brief Private methods to copy DataGroup to/from Conduit Node.
   */
  void copyToNode(Node& n) const;
  ///
  void copyFromNode(Node& n);

  /*!
   * \brief Private methods to copy DataGroup to Conduit Node.
   *
   * Vector of ids is used to maintain correct association of DataBuffers
   * to DataViews......???? punt!
   */
  void copyToNode(Node& n,
                  std::vector<IndexType>& buffer_ids) const;

  /*!
   * \brief Private methods to copy DataGroup from Conduit Node.
   *
   * Vector of ids is used to maintain correct association of DataBuffers
   * to DataViews......???? punt!
   */
  void copyFromNode(Node& n,
                    std::map<IndexType, IndexType>& id_map);

//@}


  /// Name of this DataGroup object.
  std::string m_name;

  /// Parent DataGroup of this DataGroup object.
  DataGroup * m_parent;

  /// This DataGroup object lives in the tree of this DataStore object.
  DataStore * m_datastore;

  /// Character used to denote a path string passed to get/create calls.
  static const char m_path_delimiter;

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

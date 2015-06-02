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

#if defined(USE_CXX11)
#include <unordered_map>
#else
#include "boost/unordered_map.hpp"
#endif

// Other CS Toolkit headers
#include "common/Utilities.hpp"
#include "conduit/conduit.hpp"

// SiDRe project headers
#include "Collections.hpp"
#include "DataView.hpp"
#include "SidreTypes.hpp"



namespace asctoolkit
{
namespace sidre
{

// using directives to make Conduit usage easier and less visible
using conduit::Node;
using conduit::Schema;
using conduit::DataType;

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
 *    - DataView objects are (only) created by DataGroup objects. The DataGroup
 *      that creates a DataView owns it.
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
//!  @name Basic Accessor methods

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
    DataGroup* getParent()
    {
       return m_parent;
    }

    /*!
     * \brief Return pointer to const DataGroup parent of group.
     */
    const DataGroup* getParent() const
    {
       return m_parent;
    }

    /*!
     * \brief Return pointer to non-const DataStore that owns group hierarchy
     *        to which DataGroup instance belongs.
     */
    DataStore* getDataStore()
    {
       return m_datastore;
    }

    /*!
     * \brief Return pointer to const DataStore that owns group hierarchy
     *        to which DataGroup instance belongs.
     */
    const DataStore* getDataStore() const
    {
       return m_datastore;
    }

//@}


//@{
//!  @name DataView accessor methods

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
     * \brief Return (non-const) pointer to DataView with given name.
     */
    DataView* getView( const std::string& name )
    {
        ATK_ASSERT_MSG( hasView(name), "no view found with name == " << name);
       
        return m_view_coll.getItem(name);     
    }

    /*!
     * \brief Return (const) pointer to DataView with given name.
     */
    const DataView* getView( const std::string& name ) const
    {
        ATK_ASSERT_MSG( hasView(name), "no view found with name == " << name);

        return m_view_coll.getItem(name);     
    }

    /*!
     * \brief Return (non-const) pointer to DataView with given index.
     */
    DataView* getView( IndexType idx )
    {
        ATK_ASSERT_MSG( hasView(idx), "no view found with idx == " << idx );

        return m_view_coll.getItem(idx);     
    }

    /*!
     * \brief Return (const) pointer to DataView with given index.
     */
    const DataView* getView( IndexType idx ) const
    {
        ATK_ASSERT_MSG( hasView(idx), "no view found with idx == " << idx );

        return m_view_coll.getItem(idx);     
    }

    /*!
     * \brief Return the index of DataView with given name.
     *
     *        If none, return sidre::InvalidIndex;
     */
    IndexType getViewIndex(const std::string& name) const
    {  
        ATK_ASSERT_MSG( hasView(name), "no view found with name == " << name);

        return m_view_coll.getItemIndex(name);
    }

    /*!
     * \brief Return the name of DataView with given index.
     *
     *        If none, return InvalidName.
     */
    const std::string& getViewName(IndexType idx) const
    {
        ATK_ASSERT_MSG( hasView(idx), "no view found with idx == " << idx );

        return m_view_coll.getItemName(idx);
    }

    /*!
     * \brief Return number of DataViews contained in this DataGroup.
     */
    size_t getNumViews() const
    {
        return m_view_coll.getNumItems();
    }

//@}


//@{
//!  @name DataView manipulation methods (create, destroy, copy, move, etc.)

    /*!
     * \brief Create a DataView object (and buffer) with given name and 
     *        attach to this DataGroup object.
     *
     * Note that created DataBuffer will be owned by associated DataStore.
     *
     * \return pointer to created DataView object.
     */ 
    DataView* createViewAndBuffer( const std::string& name );

    /*!
     * \brief Create a DataView object (and buffer) with given name,
     *        and use the Conduit DataType to allocate the 
     *        underlying buffer and initialize the view.
     *        Attach the new view to this DataGroup object.
     *
     * Note that created DataBuffer will be owned by associated DataStore.
     *
     * \return pointer to created DataView object.
     */ 
    DataView* createViewAndBuffer( const std::string& name, 
                                   const DataType& dtype);

    /*!
     * \brief Create a DataView object (and buffer) with given name,
     *        and use the Conduit Schema to allocate the 
     *        underlying buffer and initialize the view.
     *        Attach the new view to this DataGroup object.
     *
     * Note that created DataBuffer will be owned by associated DataStore.
     *
     * \return pointer to created DataView object.
     */ 
    DataView* createViewAndBuffer( const std::string& name, 
                                   const Schema& schema);

    /*!
     * \brief Create a DataView object (for view into given buffer) with 
     *        given name, and attach to this DataGroup object.
     *
     * \return pointer to created DataView object.
     */ 
    DataView* createView( const std::string& name, 
                          DataBuffer* buff );


    /*!
     * \brief Create a DataView object (for view into given buffer) with 
     *        given name and use the Conduit DataType to initialize the 
     *        view into this buffer. Attach the new view to this 
     *        DataGroup object.
     *
     * \return pointer to created DataView object.
     */ 
    DataView* createView( const std::string& name,
                          DataBuffer* buff, 
                          const DataType& dtype);


    /*!
     * \brief Create a DataView object (for view into given buffer) with 
     *        given name and use the Conduit Schema to initialize the 
     *        view into this buffer. Attach the new view to this 
     *        DataGroup object.
     *
     * \return pointer to created DataView object.
     */ 
    DataView* createView( const std::string& name,
                          DataBuffer* buff, 
                          const Schema& schema);

    /*!
     * \brief Create an opaque DataView with given name, holding data
     *        referenced with given pointer, and attach to this DataGroup 
     *        object.
     *
     * \return pointer to created DataView object.
     */ 
    DataView* createOpaqueView( const std::string& name, 
                                void* opaque_ptr);

    /*!
     * \brief Create an DataView into externally-owned data with given name, 
     *        using the given Conduit DataType to describe the data. Attach 
     *        the new view to this DataGroup object.
     *
     * \return pointer to created DataView object.
     */
    DataView* createExternalView( const std::string& name, 
                                  void* external_data,
                                  const DataType& dtype );
   
    /*!
     * \brief Create an DataView into externally-owned data with given name,
     *        using the given Conduit Schema to describe the data. Attach 
     *        the new view to this DataGroup object.
     *
     * \return pointer to created DataView object.
     */
    DataView* createExternalView( const std::string& name,
                                  void* external_data,
                                  const Schema& schema );
 
    /*!
     * \brief Destroy view in this DataGroup with given name and leave its
     *        associated DataBuffer intact.
     */
    void destroyView(const std::string& name);

    /*!
     * \brief Destroy view in this DataGroup with given index and leave its
     *        associated DataBuffer intact.
     */
    void destroyView(IndexType idx);

    /*!
     * \brief Destroy all views in this DataGroup and leave all associated
     *        DataBuffers intact.
     */
    void destroyViews();

    /*! 
     * \brief Destroy views in this DataGroup with given name AND destroy 
     *        its associated DataBuffer object.
     */
    void destroyViewAndBuffer(const std::string& name);

    /*! 
     * \brief Destroy view in this DataGroup with given index AND destroy 
     *        its associated DataBuffer object.
     */
    void destroyViewAndBuffer(IndexType idx);
 
    /*!
     * \brief Destroy all views in this DataGroup AND destroy their 
     *        associated DataBuffer objects.
     */
    void destroyViewsAndBuffers();

    /*!
     * \brief Remove DataView object from its owning group and attach
     *        to this DataGroup object.
     *
     * \return pointer to given DataView object.
     */ 
    DataView* moveView(DataView* view);

    /*!
     * \brief Create a copy of given DataView object and attach
     *        to this DataGroup object.
     *
     * Note that this is a "shallow" copy; the DataBuffer associated with 
     * the view is not copied. The new DataView is associated with the same
     * buffer object.
     *
     * \return pointer to created DataView object.
     */ 
    DataView* copyView(DataView* view);
    
//@}


//@{
//!  @name (child) DataGroup accessor methods

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
     * \brief Return (non-const) pointer to child DataGroup with given name.
     */
    DataGroup* getGroup( const std::string& name )
    {
        ATK_ASSERT_MSG( hasGroup(name), "no group found with name == " << name);

        return m_group_coll.getItem(name);
    }

    /*!
     * \brief Return (const) pointer to to child DataGroup with given name.
     */
    DataGroup const* getGroup( const std::string& name ) const
    {
        ATK_ASSERT_MSG( hasGroup(name), "no group found with name == " << name);

        return m_group_coll.getItem(name);
    }

    /*!
     * \brief Return (non-const) pointer to child DataGroup with given index.
     */
    DataGroup* getGroup( IndexType idx ) 
    {
        ATK_ASSERT_MSG( hasGroup(idx), "no group found with idx == " << idx );

        return m_group_coll.getItem(idx);
    }

    /*!
     * \brief Return (const) pointer to child DataGroup with given index.
     */
    const DataGroup* getGroup( IndexType idx ) const
    {
        ATK_ASSERT_MSG( hasGroup(idx), "no group found with idx == " << idx );

        return m_group_coll.getItem(idx);
    }

    /*!
     * \brief Return the index of child DataGroup with given name.
     *
     *        If none, return sidre::InvalidIndex;
     */
    IndexType getGroupIndex(const std::string& name) const
    {
        ATK_ASSERT_MSG( hasGroup(name), "no group found with name == " << name);

        return m_group_coll.getItemIndex(name);
    }

    /*!
     * \brief Return the name of child DataGroup with given index.
     *
     *        If none, return InvalidName.
     */
    const std::string& getGroupName(IndexType idx) const
    {
        ATK_ASSERT_MSG( hasGroup(idx), "no group found with idx == " << idx );

        return m_group_coll.getItemName(idx);
    }

    /*!
     * \brief Return number of child DataGroups contained in this DataGroup.
     */
    size_t getNumGroups() const
    {
        return m_group_coll.getNumItems();
    }

//@}


//@{
//!  @name DataGroup manipulation methods (create, destroy, copy, move, etc.)

    /*!
     * \brief Create a child DataGroup object with given name in this group.
     *
     * \return pointer to created DataGroup object.
     */
    DataGroup* createGroup( const std::string& name );
  
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

    /*!
     * \brief Remove DataGroup object from its parent group and attach
     *        to this DataGroup object.
     *
     * \return pointer to given DataGroup object.
     */
    DataGroup* moveGroup(DataGroup* group);

    /*!
     * \brief Create a copy of given DataGroup object (including all of its 
     *        DataViews and child DataGroups) and attach to this DataGroup 
     *        object.
     *
     * Note that DataView copying is a "shallow" copy; the DataBuffer 
     * associated with a view is not copied. The new DataView is associated 
     * with the same buffer object.
     *
     * \return pointer to created DataGroup object.
     */
    DataGroup* copyGroup(DataGroup* group);

//@}

 
    /*!
     * \brief Copy data group description to given Conduit node.
     */  
    void info(Node& n) const;

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
     *        starting at this DataGroup object to stdout.
     */
    void printTree( const int nlevels ) const;
 

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
 
private:

    /*!
     *  \brief Private ctor that creates a Group with given name
     *         in given parent group.
     */
    DataGroup(const std::string& name, DataGroup* parent);

    /*!
     *  \brief Private ctor that creates a Group with given name
     *         in the given DataStore root group.
     */
    DataGroup(const std::string& name, DataStore* datastore);

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
    
    /*!
     * \brief Private methods to attach/detach DataView object to DataGroup.
     */ 
    DataView* attachView(DataView* view);
    ///
    DataView* detachView(const std::string& name);
    ///
    DataView* detachView(IndexType idx);

    /*!
     * \brief Private methods to attach/detach DataGroup object to DataGroup.
     */ 
    DataGroup* attachGroup(DataGroup* group);
    ///
    DataGroup* detachGroup(const std::string& name);
    ///
    DataGroup* detachGroup(IndexType idx);

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


   
    /// Name of this DataGroup object.
    std::string  m_name;

    /// Parent DataGroup of this DataGroup object.
    DataGroup* m_parent;

    /// This DataGroup object lives in the tree of this DataStore object.
    DataStore* m_datastore;


    ///
    /// Typedefs for view and shild group containers. They are here to
    /// avoid propagating specific type names in the DataGroup class
    /// implementation when we experiment with different containers.
    ///
    ///////////////////////////////////////////////////////////////////
    // 
    // Map container options
    //
    // To try a different container, set the "MapType" typedef to 
    // what you want.  Note: only one typedef should be active!!!
    //
    // Current options are std::map and boost/std::unordered_map
    //
       typedef std::map<std::string, IndexType> MapType;
    //
#if defined(USE_CXX11)
    // typedef std::unordered_map<std::string, IndexType> MapType;
#else
    // typedef boost::unordered_map<std::string, IndexType> MapType;
#endif
    ///
    typedef MapCollection<DataView, MapType> DataViewCollection;
    ///
    typedef MapCollection<DataGroup, MapType> DataGroupCollection;
    // 
    ///////////////////////////////////////////////////////////////////

 
    /// Collection of DataViews
    DataViewCollection m_view_coll;

    /// Collection of child DataGroups
    DataGroupCollection m_group_coll;

};


} /* end namespace sidre */
} /* end namespace asctoolkit */

#endif /* DATAGROUP_HPP_ */

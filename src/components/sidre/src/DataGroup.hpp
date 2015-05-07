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

// Other CS Toolkit headers
#include "common/Types.hpp"
#include "common/Utilities.hpp"

// SiDRe project headers
#include "Collections.hpp"
#include "DataView.hpp"


// using directives to make Conduit usage easier and less visible
using conduit::index_t;



namespace asctoolkit
{
namespace sidre
{

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
    DataGroup const* getParent() const
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
    DataStore const* getDataStore() const
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
    bool hasView( common::IDType idx ) const
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
    DataView const* getView( const std::string& name ) const
    {
        ATK_ASSERT_MSG( hasView(name), "no view found with name == " << name);

        return m_view_coll.getItem(name);     
    }

    /*!
     * \brief Return (non-const) pointer to DataView with given index.
     */
    DataView* getView( common::IDType idx )
    {
        ATK_ASSERT_MSG( hasView(idx), "no view found with idx == " << idx );

        return m_view_coll.getItem(idx);     
    }

    /*!
     * \brief Return (const) pointer to DataView with given index.
     */
    DataView const* getView( common::IDType idx ) const
    {
        ATK_ASSERT_MSG( hasView(idx), "no view found with idx == " << idx );

        return m_view_coll.getItem(idx);     
    }

#if 0  // 
       // RDH -- These should probably be removed since they are error-prone.
       //      
       // I compiled them out for now since they are not used anywhere. 
       //
       // For example, using the map "at" method throws an exception if no 
       // item is found with a given key (here, a name). In that case, the 
       // returned index value is undefined. We could replace the "at" usage 
       // with "find" and then return an "undefined" index value that a caller 
       // could test for. However, common::IDType we are using is unsigned. So 
       // there is no clear choice for an "undefined" index.  
       //
       // Similarly, if there is no item with a given index, then what should 
       // the returned name be? An empty string? Again, a caller would have to 
       // check against this.  
       //
       // IMO, users should always know whether they should use a name or 
       // index for lookup. We have routines to check whether an item exists 
       // with a name or index.
       //
    /*!
     * \brief Return the index of DataView with given name.
     */
    common::IDType getViewIndex(const std::string &name) const
    {  
        ATK_ASSERT_MSG( hasView(name), "no view found with name == " << name);

        return m_viewsNameMap.at(name);
    }

    /*!
     * \brief Return the name of DataView with given index.
     */
    const std::string& getViewName(common::IDType idx) const
    {
        ATK_ASSERT_MSG( hasView(idx), "no view found with idx == " << idx );

        const DataView* view = getView(idx);
        return view->getName();
    }
#endif

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
     * \brief Create an opaque DataView with given name, holding data
     *        referenced with given pointer, and attach to this DataGroup 
     *        object.
     *
     * \return pointer to created DataView object.
     */ 
    DataView* createOpaqueView( const std::string& name, void* opaque_ptr);
    
    /*!
     * \brief Create a DataView object (for view into given buffer) with 
     *        given name, and attach to this DataGroup object.
     *
     * \return pointer to created DataView object.
     */ 
    DataView* createView( const std::string& name, DataBuffer* buff );

    /*!
     * \brief Destroy view in this DataGroup with given name and leave its
     *        associated DataBuffer intact.
     */
    void destroyView(const std::string &name);

    /*!
     * \brief Destroy view in this DataGroup with given index and leave its
     *        associated DataBuffer intact.
     */
    void destroyView(common::IDType idx);

    /*!
     * \brief Destroy all views in this DataGroup and leave all associated
     *        DataBuffers intact.
     */
    void destroyViews();

    /*! 
     * \brief Destroy views in this DataGroup with given name AND destroy 
     *        its associated DataBuffer object.
     */
    void destroyViewAndBuffer(const std::string &name);

    /*! 
     * \brief Destroy view in this DataGroup with given index AND destroy 
     *        its associated DataBuffer object.
     */
    void destroyViewAndBuffer(common::IDType idx);
 
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
    DataView* copyView(DataView *view);
    
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
    bool hasGroup( common::IDType idx ) const
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
    DataGroup* getGroup( common::IDType idx ) 
    {
        ATK_ASSERT_MSG( hasGroup(idx), "no group found with idx == " << idx );

        return m_group_coll.getItem(idx);
    }

    /*!
     * \brief Return (const) pointer to child DataGroup with given index.
     */
    DataGroup const* getGroup( common::IDType idx ) const
    {
        ATK_ASSERT_MSG( hasGroup(idx), "no group found with idx == " << idx );

        return m_group_coll.getItem(idx);
    }

#if 0  // 
       //  RDH -- These should probably be removed since they are error-prone.
       //
       // See comments above for similar view methods.
       //
    /*!
     * \brief Return the index of child DataGroup with given name.
     */
    common::IDType getGroupIndex(const std::string &name) const
    {
        ATK_ASSERT_MSG( hasGroup(name), "no group found with name == " << name);

        return m_groupsNameMap.at(name);
    }

    /*!
     * \brief Return the name of child DataGroup with given index.
     */
    const std::string& getGroupName(common::IDType idx) const
    {
        ATK_ASSERT_MSG( hasGroup(idx), "no group found with idx == " << idx );

        const DataGroup* group = getGroup(idx);
        return group->getName();
    }
#endif

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
    void destroyGroup(const std::string &name);

    /*!
     * \brief Destroy child group in this DataGroup with given index.
     */
    void destroyGroup(common::IDType idx); 

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
    void info(Node &n) const;

    /*!
     * \brief Print JSON description of data group to stdout.
     *
     * Note that this will recursively print entire group (sub) tree
     * starting at this DataGroup object.
     */
    void print() const;

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
    DataView* detachView(common::IDType idx);

    /*!
     * \brief Private methods to attach/detach DataGroup object to DataGroup.
     */ 
    DataGroup* attachGroup(DataGroup* group);
    ///
    DataGroup* detachGroup(const std::string& name);
    ///
    DataGroup* detachGroup(common::IDType idx);

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
                    std::vector<common::IDType>& buffer_ids) const;

    /*!
     * \brief Private methods to copy DataGroup from Conduit Node.
     *
     * Vector of ids is used to maintain correct association of DataBuffers
     * to DataViews......???? punt!
     */  
    void copyFromNode(Node& n,
                      std::map<common::IDType, common::IDType>& id_map);


   
    /// Name of this DataGroup object.
    std::string  m_name;

    /// Parent DataGroup of this DataGroup object.
    DataGroup* m_parent;

    /// This DataGroup object lives in the tree of this DataStore object.
    DataStore* m_datastore;


    ///
    /// Typedefs for view and shild group containers. They are here to
    /// avoid propagating specific type names in the group implementation.
    ///
    typedef MapCollection<DataView> DataViewCollection;
    ///
    typedef MapCollection<DataGroup> DataGroupCollection;
 
    /// Collection of DataViews
    DataViewCollection m_view_coll;

    /// Collection of child DataGroups
    DataGroupCollection m_group_coll;

};


} /* end namespace sidre */
} /* end namespace asctoolkit */

#endif /* DATAGROUP_HPP_ */

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

// SiDRe project headers
#include "DataView.hpp"
#include "Types.hpp"
#include "Utilities.hpp"


// using directives to make Conduit usage easier and less visible
using conduit::index_t;


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
    bool hasView( const std::string& name );

    /*!
     * \brief Return (non-const) pointer to DataView with given name.
     */
    DataView* getView( const std::string& name )
    {
        ATK_ASSERT_MSG( m_viewsNameMap.find(name) != m_viewsNameMap.end(), \
                        "no view found with name == " << name);
            
        const IDType idx = m_viewsNameMap.at(name);
        return m_views[idx];
    }

    /*!
     * \brief Return (const) pointer to DataView with given name.
     */
    DataView const* getView( const std::string& name ) const
    {
        ATK_ASSERT_MSG( m_viewsNameMap.find(name) != m_viewsNameMap.end(), \
                        "no view found with name == " << name);

        const IDType idx = m_viewsNameMap.at(name);
        return m_views[idx];
    }

    /*!
     * \brief Return (non-const) pointer to DataView with given index.
     */
    DataView* getView( const IDType idx )
    {
        ATK_ASSERT_MSG( idx >= 0 && idx < m_views.size(), \
                        "no view found with idx == ");
        // TODO: add "idx" to error message

        return m_views[idx];
    }

    /*!
     * \brief Return (const) pointer to DataView with given index.
     */
    DataView const* getView( const IDType idx ) const
    {
        ATK_ASSERT_MSG( idx >= 0 && idx < m_views.size(), \
                        "no view found with idx == ");

        return m_views[idx];
    }

    /*!
     * \brief Return the index of DataView with given name.
     */
    IDType getViewIndex(const std::string &name) const
    {  
        ATK_ASSERT_MSG( m_viewsNameMap.find(name) != m_viewsNameMap.end(), \
                        "no view found with name == " << name);

        return m_viewsNameMap.at(name);
    }

    /*!
     * \brief Return the name of DataView with given index.
     */
    const std::string& getViewName(IDType idx) const
    {
        const DataView* view = getView(idx);
        return view->getName();
    }
  
    /*!
     * \brief Return number of DataViews contained in this DataGroup.
     */
    size_t getNumberOfViews() const
    {
        return m_views.size();
    }

//@}


//@{
//!  @name DataView create, destroy, copy, move methods

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
    void destroyView(IDType idx);

    /*!
     * \brief Destroy all views in this DataGroup and leave all associated
     *        DataBuffers intact.
     */
    void destroyViews();

    /*! \brief Destroy views in this DataGroup with given name AND destroy 
     *         its associated DataBuffer object.
     */
    void destroyViewAndBuffer(const std::string &name);

    /*! \brief Destroy view in this DataGroup with given index AND destroy 
     *         its associated DataBuffer object.
     */
    void destroyViewAndBuffer(IDType idx);
 
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
     * \return pointer to given DataView object.
     */ 
    DataView* copyView(DataView *view);
    
//@}


//@{
//!  @name (child) DataGroup accessor methods

    /// -----  DataGroup Children ---- /// 
    bool hasGroup( const std::string& name );

    /*!
    * @param name Name of DataGroup to find.
    * \brief Return pointer to DataGroup.
    */
    DataGroup const * getGroup( const std::string& name ) const
    {
      const IDType idx = m_groupsNameMap.at(name);
      return m_groups[idx];
    }

    DataGroup * getGroup( const std::string& name )
    {
      const IDType idx = m_groupsNameMap.at(name);
      return m_groups[idx];
    }

    /*!
    * @param idx Index of DataGroup to find.
    * \brief Return pointer to DataGroup.
    */
    DataGroup const * getGroup(IDType idx) const
    {
     return m_groups[idx];
    }

    DataGroup * getGroup( IDType idx)
    {
     return m_groups[idx];
    }

    /*!
    * \brief Return the index of the DataGroup with the given name
    */
    IDType getGroupIndex(const std::string &name) const
    {
       return m_groupsNameMap.at(name);
    }

    /*!
    * \brief Return the name of the DataGroup at the given index
    */
    std::string GetGroupName(IDType idx) const
    {
      return m_views[idx]->getName();
    }

    /*!
    * \brief Return number of DataGroups contained in this DataGroup.
    */
    size_t getNumberOfGroups() const
    {
    return m_groups.size();
    }

//@}


//@{
//!  @name DataGoup create, destroy, copy, move methods

    /*!
    * @param name Name of DataGroup to create.
    * \brief Create a new DataGroup within this DataGroup.
    */
    DataGroup* createGroup( const std::string& name );
    
    // removes a group from another group into this group
    // returns `grp`
    DataGroup *moveGroup(DataGroup *grp);
    // creates a copy of the given group into this group
    // this will also copy all sub groups and views. 
    // Recall:copying the views does not imply copying the buffers.
    // returns the new group
    DataGroup *copyGroup(DataGroup *grp);


    void destroyGroup(const std::string &name);
    void destroyGroup(IDType idx);

    /*!
    * \brief Remove all DataViews from this DataGroup.
    */
    void destroyGroups();

//@}


    void info(Node &n) const;
    void print() const;

    void printTree( const int level ) const;
 
 
    /// ---------------------------------------------------------------
    ///  Save + Restore Prototypes (ATK-39)
    /// ---------------------------------------------------------------
    /// saves "this", associated views and buffers to a file set. 
    void save(const std::string &obase,
              const std::string &protocol) const;

    /// restores as "this"
    void load(const std::string &obase,
              const std::string &protocol);
 
private:

    /// these are private b/c we want folks to create groups
    /// from another group or a  datastore
    DataGroup(const std::string &name, DataGroup *parent);
    DataGroup(const std::string &name, DataStore *datastore);

    /*!
    * @param source
    * \brief default copy constructor
    */
    DataGroup( const DataGroup& source );

    /*!
    *
    * @param rhs the DataView to be copied
    * @return *this
    */
    DataGroup& operator=( const DataGroup& rhs );

#ifdef USE_CXX11
  /*!
    * @param source
    * \brief default move constructor
    */
    DataGroup( DataGroup&& source );

    /*!
    *
    * @param rhs the DataView to be moved into *this
    * @return *this
    */
    DataGroup& operator=( const DataGroup&& rhs );
#endif
    
    
    /*!
    * \brief destructor
    */
    ~DataGroup();
    
    
    /// Attach + Detach are private since they have scary 
    /// bookkeeping side affects.
    
    /// Our use cases should be supported by:
    ///  CreateView|Group()
    ///  MoveView|Group()
    ///  CopyView|Group()
    ///  DestroyView|Group()
 
    DataView *attachView(DataView *view);
    DataView *detachView(const std::string &name);
    DataView *detachView(IDType idx);


    DataGroup *attachGroup(DataGroup *grp);
    DataGroup *detachGroup(const std::string &name);
    DataGroup *detachGroup(IDType idx);

    ///
    /// there may be value to make these public
    ///
    
    void copyToNode(Node &n) const;
    void copyFromNode(Node &n);
    
    ///
    /// these should stay private
    ///
    void copyToNode(Node &n,
                    std::vector<IDType> &buffer_ids) const;

    /// we could use an unordered map to track the id mapping
    void copyFromNode(Node &n,
                      std::map<IDType,IDType> &id_map);

    
    std::string  m_name;
    DataGroup   *m_parent;
    DataStore   *m_datastore;

    std::vector<DataView*>       m_views;
    std::map<std::string,IDType> m_viewsNameMap;

    std::vector<DataGroup*>      m_groups;
    std::map<std::string,IDType> m_groupsNameMap;


};

} /* namespace sidre */
#endif /* DATAGROUP_HPP_ */

/*
 * DataGroup.hpp
 *
 *  Created on: Dec 2, 2014
 *      Author: settgast
 */

#ifndef DATAGROUP_HPP_
#define DATAGROUP_HPP_

#include "DataView.hpp"
#include <memory>
#include <map>
#include <vector>
#include "Types.hpp"

#include "conduit/conduit.h"

#include "Utilities.hpp"


using conduit::index_t;

namespace sidre
{

/**
 * \class DataGroup
 *
 * \brief Class to access collections of DataViews.
 *  The DataGroup will name each DataView as it is added.
 */
class DataGroup
{
public:
    friend class DataStore;

    /// -----  Basic Members  ---- /// 
    std::string getName() const
    {return m_name; }

    DataGroup  *getParent()
    {return m_parent;}

    DataGroup const *getParent() const
    {return m_parent;}

    DataStore *getDataStore()
    {return m_datastore;}

    DataStore const *getDataStore() const
    {return m_datastore;}


    /// -----  DataView Children ---- /// 
    bool hasView( const std::string& name );

    /*!
    * @param name Name for created DataView.
    * \brief Create a DataView and add to this DataGroup.
    */
    
    DataView *createViewAndBuffer( const std::string& name );
    DataView *createOpaqueView( const std::string& name, void *);
    
    DataView *createView( const std::string& name, DataBuffer *buff);


    // removes a view from another group into this group
    // returns `view`
    DataView *moveView(DataView *view);
    // creates a copy of the given view for this group
    // Recall:copying the view does not imply copying the buffer.
    // returns the new view
    DataView *copyView(DataView *view);
    
    
    void destroyViewAndBuffer(const std::string &name);
    void destroyViewAndBuffer(IDType idx);
    
    void destroyView(const std::string &name);
    void destroyView(IDType idx);
 
    /*!
    * @param name Name of DataView to find.
    * \brief Return pointer to DataView.
    */
    DataView *getView( const std::string& name )
    {
        ATK_ASSERT_MSG( m_viewsNameMap.find(name) != m_viewsNameMap.end(), "GetView() tried to fetch invalid view named ");
        // TODO: add "name" to error message, I had problems doing this with the macro
            
        const IDType idx = m_viewsNameMap.at(name);
        return m_views[idx];
    }

    DataView const * getView( const std::string& name ) const
    {
        ATK_ASSERT_MSG( m_viewsNameMap.find(name) != m_viewsNameMap.end(), "GetView() tried to fetch invalid view named ");
        // TODO: add "name" to error message, I had problems doing this with the macro
        const IDType idx = m_viewsNameMap.at(name);
        return m_views[idx];
    }

    /*!
    * @param idx Index of DataView within this DataGroup.
    * \brief Return pointer to DataView.
    */
    DataView *getView( const IDType idx )
    {
        ATK_ASSERT_MSG( idx >= 0 && idx < m_views.size(), "GetView() tried to fetch view at invalid index ");
        // TODO: add "idx" to error message, I had problems doing this with the macro
        return m_views[idx];
    }

    /*!
    * @param idx Index of DataView within this DataGroup.
    * \brief Return pointer to DataView.
    */
    DataView const *getView( const IDType idx ) const
    {
        ATK_ASSERT_MSG( idx >= 0 && idx < m_views.size(), "GetView() tried to fetch view at invalid index ");
        // TODO: add "idx" to error message, I had problems doing this with the macro
        return m_views[idx];
    }

    /*!
    * \brief Return the index of the DataView with the given name
    */
    IDType getViewIndex(const std::string &name) const
    {  
      return m_viewsNameMap.at(name);
    }

    /*!
    * \brief Return the name of the DataView at the given index
    */
    std::string getViewName(IDType idx) const
    {
       return m_views[idx]->getName();
    }
  
    /*!
    * \brief Return number of DataViews contained in this DataGroup.
    */
    size_t getNumberOfViews() const
    {
      return m_views.size();
    }

    /*!
    * \brief Remove all view from this group.
    */
    void destroyViews();

    /*!
    * \brief Remove all views from this group and destroy their buffers.
    */
    void destroyViewsAndBuffers();

    /// -----  DataGroup Children ---- /// 
    bool hasGroup( const std::string& name );

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

    /*!
    * \brief Remove all DataViews from this DataGroup.
    */
    void destroyGroups();

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

#ifdef USECXX11
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

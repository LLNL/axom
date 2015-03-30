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

using conduit::index_t;

namespace DataStoreNS
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
    /*!
    * \brief destructor
    */
    ~DataGroup();


    /// -----  Basic Members  ---- /// 
    std::string GetName() const
    {return m_name; }

    DataGroup  *GetParent()
    {return m_parent;}

    DataGroup const *GetParent() const
    {return m_parent;}

    DataStore *GetDataStore()
    {return m_datastore;}

    DataStore const *GetDataStore() const
    {return m_datastore;}


    /*!
    * @param name Name to check.
    * \brief Return true if the name exists in this DataGroup.
    */
    
    /// we can have child groups and views
    bool HasChild( const std::string& name );

    /// -----  DataView Children ---- /// 
    bool HasView( const std::string& name );

    /*!
    * @param name Name for created DataView.
    * \brief Create a DataView and add to this DataGroup.
    */
    DataView *CreateView( const std::string& name );
    DataView *CreateView( const std::string& name, DataBuffer *buff);


    // removes a view from another group into this group
    // returns `view`
    DataView *MoveView(DataView *view);
    // creates a copy of the given view for this group
    // Recall:copying the view does not imply copying the buffer.
    // returns the new view
    DataView *CopyView(DataView *view);
    
    
    void DestroyView(const std::string &name);
    void DestroyView(IDType idx);
    void DestroyView(DataView *view);
    
    

    /*!
    * @param name Name of DataView to find.
    * \brief Return pointer to DataView.
    */
    DataView *GetView( const std::string& name )
    {
        const IDType idx = m_viewsNameMap.at(name);
        return m_views[idx];
    }

    DataView const * GetView( const std::string& name ) const
    {
        const IDType idx = m_viewsNameMap.at(name);
        return m_views[idx];
    }

    /*!
    * @param idx Index of DataView within this DataGroup.
    * \brief Return pointer to DataView.
    */
    DataView *GetView( const IDType idx )
    {
          return m_views[idx];
    }

    /*!
     *
     * @return
     */
    std::map<std::string,IDType> const& GetViewsNameMap() const
    {
      return m_viewsNameMap;
    }

    /*!
    * @param idx Index of DataView within this DataGroup.
    * \brief Return pointer to DataView.
    */
    DataView const *GetView( const IDType idx ) const
    {
      return m_views[idx];
    }

    /*!
     *
     * @return
     */
    std::vector<DataView*>& GetViews()
    {
      return m_views;
    }

    /*!
     *
     * @return
     */
    std::vector<DataView*> const & GetViews() const
    {
      return m_views;
    }

    /*!
    * \brief Return the index of the DataView with the given name
    */
    IDType GetViewIndex(const std::string &name) const
    {  
      return m_viewsNameMap.at(name);
    }


    /*!
    * \brief Return the name of the DataView at the given index
    */
    std::string GetViewName(IDType idx) const
    {
       return m_views[idx]->GetName();
    }
  
    /*!
    * \brief Return number of DataViews contained in this DataGroup.
    */
    size_t CountViews() const
    {
      return m_views.size();
    }

    /*!
    * \brief Remove all DataViews from this DataGroup.
    */
    void DestroyViews();

    /// -----  DataGroup Children ---- /// 
    bool HasGroup( const std::string& name );
  
    /*!
    * @param name Name of DataGroup to create.
    * \brief Create a new DataGroup within this DataGroup.
    */
    DataGroup* CreateGroup( const std::string& name );
    
    // removes a group from another group into this group
    // returns `grp`
    DataGroup *MoveGroup(DataGroup *grp);
    // creates a copy of the given group into this group
    // this will also copy all sub groups and views. 
    // Recall:copying the views does not imply copying the buffers.
    // returns the new group
    DataGroup *CopyGroup(DataGroup *grp);


    void DestroyGroup(const std::string &name);
    void DestroyGroup(IDType idx);

    /*!
    * @param name Name of DataGroup to find.
    * \brief Return pointer to DataGroup.
    */
    DataGroup const * GetGroup( const std::string& name ) const
    {
      const IDType idx = m_groupsNameMap.at(name);
      return m_groups[idx];
    }

    DataGroup * GetGroup( const std::string& name )
    {
      const IDType idx = m_groupsNameMap.at(name);
      return m_groups[idx];
    }

    /*!
    * @param idx Index of DataGroup to find.
    * \brief Return pointer to DataGroup.
    */
    DataGroup const * GetGroup(IDType idx) const
    {
     return m_groups[idx];
    }

    DataGroup * GetGroup( IDType idx)
    {
     return m_groups[idx];
    }

    /*!
    * \brief Return the index of the DataGroup with the given name
    */
    IDType GetGroupIndex(const std::string &name) const
    {
       return m_groupsNameMap.at(name);
    }

    /*!
    * \brief Return the name of the DataGroup at the given index
    */
    std::string GetGroupName(IDType idx) const
    {
      return m_views[idx]->GetName();
    }


    /*!
    * \brief Return number of DataGroups contained in this DataGroup.
    */
    size_t CountGroups() const
    {
    return m_groups.size();
    }

    /*!
    * \brief Remove all DataViews from this DataGroup.
    */
    void DestroyGroups();

    void Info(Node &n) const;
    void Print() const;

    void PrintTree( const int level ) const;
 
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
    
    
    /// Attach + Detach are private since they have scary 
    /// bookkeeping side affects.
    
    /// Our use cases should be supported by:
    ///  CreateView|Group()
    ///  MoveView|Group()
    ///  CopyView|Group()
    ///  DestroyView|Group()
 
    DataView *AttachView(DataView *view);
    DataView *DetachView(const std::string &name);
    DataView *DetachView(IDType idx);


    DataGroup *AttachGroup(DataGroup *grp);
    DataGroup *DetachGroup(const std::string &name);
    DataGroup *DetachGroup(IDType idx);

    
    std::string  m_name;
    DataGroup   *m_parent;
    DataStore   *m_datastore;

    std::vector<DataView*>       m_views;
    std::map<std::string,IDType> m_viewsNameMap;

    std::vector<DataGroup*>      m_groups;
    std::map<std::string,IDType> m_groupsNameMap;


};

} /* namespace DataStore */
#endif /* DATAGROUP_HPP_ */

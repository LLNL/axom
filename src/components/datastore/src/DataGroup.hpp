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

namespace DataStoreNS
{

/**
 * \class DataGroup
 *
 * \brief Class to access collections of DataView.
 *  The DataGroup will name each DataView as it is added.
 */
class DataGroup
{
public:
  /*!
   * \brief vector of DataView pointers.
   */
  typedef std::vector<DataView*> dataArrayType;

  /*!
   * \brief map of name to index of DataView within this DataGroup.
   */
  typedef std::map<std::string, IDType> lookupType;

  /*!
   * \brief map of name to DataGroup pointer.
   */
  typedef std::map<std::string, DataGroup*> lookupGroup;

private:
  DataGroup *m_parent;
  DataStore *m_datastore;
  dataArrayType m_DataViews;  // DataViews by index
  lookupType m_DataViewLookup;      // DataViews name to View pointer
  lookupGroup m_childGroups;  // child Groups: name->Group pointer
  std::string m_name;

#if 0
  DataShape m_dataShape;
#endif

public:

  /*!
   * @param parent name Pointer to DataGroup which contains this Group.
   * @param datastore Pointer to DataStore container.
   * \brief Constructor.
   */
  DataGroup( DataGroup *parent, DataStore *datastore ) :
      m_parent(parent), m_datastore(datastore)
  {
  }



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








  /*!
   * @param name Name to check.
   * \brief Return true if the name exists in this DataGroup.
   */
  bool HasName( const std::string& name );

  /*!
   * @param name Name for created DataView.
   * \brief Create a DataView and add to this DataGroup.
   */
  DataView *CreateDataView( const std::string& name );

  /*!
   * @param name Name of DataView to add.
   * @param obj  Pointer to an existing DataView.
   * \brief Add existing DataView to this DataGroup.
   */
  DataView *AttachDataView( const std::string& name, DataView *obj );

  DataView* DetatchDataView( const std::string& name );


  /*!
   * @param name Name of DataView to find.
   * \brief Return pointer to DataView.
   */
  DataView *GetDataView( const std::string& name )
  {
    const IDType indx = m_DataViewLookup.at(name);
    return m_DataViews[indx];
  }
  DataView const * GetDataView( const std::string& name ) const
  {
    const IDType indx = m_DataViewLookup.at(name);
    return m_DataViews[indx];
  }

  /*!
   * @param indx Index of DataView within this DataGroup.
   * \brief Return pointer to DataView.
   */
  DataView *GetDataView( const IDType indx )
  {
    DataView *obj = m_DataViews[indx];
    if( obj == NULL )
    {
      // View has been deleted and index is a hole in the table.
      throw std::exception();
    }
    return obj;
  }

  /*!
   * @param name Name of DataView to find.
   * \brief Return index of DataView in this DataGroup.
   */
  IDType IndexDataView( const std::string& name )
  {
    return m_DataViewLookup.at(name);
  }

  /*!
   * @param obj Name of DataView to find.
   * \brief Return name of DataView in this DataGroup.
   */
  std::string const & NameDataView( DataView *obj );

  /*!
   * @param name Name of DataView to remove.
   * \brief Remove named DataView from the index.
   *   The DataView still exists in the DataStore.
   */
  void RemoveDataView( const std::string& name );

  /*!
   * @param name Name of DataGroup to create.
   * \brief Create a new DataGroup within this DataGroup.
   */
  DataGroup* CreateDataGroup( const std::string& name );

  /*!
   * @param name Name of DataGroup to find.
   * \brief Return pointer to DataGroup.
   */
  DataGroup const * GetDataGroup( const std::string& name ) const
  {
    return m_childGroups.at(name);
  }

  DataGroup * GetDataGroup( const std::string& name )
  {
    return m_childGroups.at(name);
  }

  /*!
   * \brief Return number of DataViews contained in this DataGroup.
   */
  size_t CountViews()
  {
    return m_DataViews.size();
  }

  /*!
   * \brief Return number of DataGroups contained in this DataGroup.
   */
  size_t CountGroups()
  {
    return m_childGroups.size();
  }

  /*!
   * \brief Return DataViews contained in this DataGroup.
   */
  lookupType const & GetDataViewLookup() const
  {
    return m_DataViewLookup;
  }

  dataArrayType const & GetDataViews() const
  {
    return m_DataViews;
  }

  /*!
   * \brief Return DataGroups contained in this DataGroup.
   */
  lookupGroup& GetDataGroups()
  {
    return m_childGroups;
  }




  /**
   * @name members that will be deprecated by convenience layer
   */
  ///@{

  DataShape m_dataShape;
  /*!
   *
   * @param dataDescriptor
   * @return
   */
  DataGroup* SetDataShape( const DataShape& dataShape )
  {
    m_dataShape = dataShape;
    return this;
  }


  const DataShape& GetDataShape() const
  {
    return m_dataShape;
  }



  ///@}


};

} /* namespace DataStore */
#endif /* DATAGROUP_HPP_ */

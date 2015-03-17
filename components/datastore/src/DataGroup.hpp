/*
 * DataGroup.hpp
 *
 *  Created on: Dec 2, 2014
 *      Author: settgast
 */

#ifndef DATAGROUP_HPP_
#define DATAGROUP_HPP_

#include "DataObject.hpp"
#include <memory>
#include <map>
#include <vector>
#include "Types.hpp"

namespace DataStoreNS
{

/**
 * \class DataGroup
 *
 * \brief Class to access collections of DataObject.
 *  The DataGroup will name each DataObject as it is added.
 */
class DataGroup
{
public:
  /*!
   * \brief vector of DataObject pointers.
   */
  typedef std::vector<DataObject*> dataArrayType;

  /*!
   * \brief map of name to index of DataObject within this DataGroup.
   */
  typedef std::map<std::string, IDType> lookupType;

  /*!
   * \brief map of name to DataGroup pointer.
   */
  typedef std::map<std::string, DataGroup*> lookupGroup;

private:
  DataGroup *m_parent;
  DataStore *m_datastore;
  dataArrayType m_DataObjects;  // DataObjects by index
  lookupType m_DataObjectLookup;      // DataObjects name to Object pointer
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
   * @param rhs the DataObject to be copied
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
   * @param rhs the DataObject to be moved into *this
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
   * @param name Name for created DataObject.
   * \brief Create a DataObject and add to this DataGroup.
   */
  DataObject *CreateDataObject( const std::string& name )
  {
    return AddDataObject(name, NULL);
  }

  /*!
   * @param name Name of DataObject to add.
   * @param obj  Pointer to an existing DataObject.
   * \brief Add existing DataObject to this DataGroup.
   */
  DataObject *AddDataObject( const std::string& name, DataObject *obj );

  /*!
   * @param name Name of DataObject to find.
   * \brief Return pointer to DataObject.
   */
  DataObject *GetDataObject( const std::string& name )
  {
    const IDType indx = m_DataObjectLookup.at(name);
    return m_DataObjects[indx];
  }
  DataObject const * GetDataObject( const std::string& name ) const
  {
    const IDType indx = m_DataObjectLookup.at(name);
    return m_DataObjects[indx];
  }

  /*!
   * @param indx Index of DataObject within this DataGroup.
   * \brief Return pointer to DataObject.
   */
  DataObject *GetDataObject( const IDType indx )
  {
    DataObject *obj = m_DataObjects[indx];
    if( obj == NULL )
    {
      // Object has been deleted and index is a hole in the table.
      throw std::exception();
    }
    return obj;
  }

  /*!
   * @param name Name of DataObject to find.
   * \brief Return index of DataObject in this DataGroup.
   */
  IDType IndexDataObject( const std::string& name )
  {
    return m_DataObjectLookup.at(name);
  }

  /*!
   * @param obj Name of DataObject to find.
   * \brief Return name of DataObject in this DataGroup.
   */
  std::string const & NameDataObject( DataObject *obj );

  /*!
   * @param name Name of DataObject to remove.
   * \brief Remove named DataObject from the index.
   *   The DataObject still exists in the DataStore.
   */
  DataObject *RemoveDataObject( const std::string& name );

  /*!
   * @param obj Pointer to DataObject to remove.
   * \brief Remove obj from the index.
   *   The DataObject still exists in the DataStore.
   */
  DataObject *RemoveDataObject( DataObject *obj );

  /*!
   * \brief Remove all DataObjects from this DataGroup.
   */
  void ClearDataObjects();

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
   * \brief Return number of DataObjects contained in this DataGroup.
   */
  size_t CountObjects()
  {
    return m_DataObjects.size();
  }

  /*!
   * \brief Return number of DataGroups contained in this DataGroup.
   */
  size_t CountGroups()
  {
    return m_childGroups.size();
  }

  /*!
   * \brief Return DataObjects contained in this DataGroup.
   */
  lookupType const & GetDataObjectLookup() const
  {
    return m_DataObjectLookup;
  }

  dataArrayType const & GetDataObjects() const
  {
    return m_DataObjects;
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

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
#include "Utilities.hpp"
#include <vector>
#include "Types.hpp"

namespace DataStoreNS
{



/**
 * \class DataGroup
 *
 * \brief class to access (and potentially own) collections of DataObjects
 */
class DataGroup
{
public:
  // change to regular pointer for now.
  typedef std::vector< DataObject* > dataArrayType;
  typedef std::unordered_map<std::string, sizet> lookupType;


private:
  DataGroup* m_parent;

  DataStore* m_dataStore;

  /// array of dataObject*
  dataArrayType m_DataObjects;

  /// hash table for index lookup by name
  lookupType m_DataObjectLookup;

  std::unordered_map<std::string, DataGroup* > m_childGroups;

  std::string m_name;

  DataShape m_dataShape;


public:

  /**
   * @name Constructor, Destructor, Assignment, Move, Copy
   */
  ///@{

  /// non-callable default constructor
  DataGroup() = delete;

  /*!
   *
   * @param name name of the group
   * @param path path
   * \brief constructor
   */
  DataGroup( const std::string& name,
             const std::string& path,
             DataGroup* const parent,
             DataStore* const dataStore );

  /*!
   * @param source
   * \brief default copy constructor
   */
  DataGroup( const DataGroup& source );

  /*!
   * @param source
   * \brief default move constructor
   */
  DataGroup( DataGroup&& source );

  /*!
   * \brief default destructor
   */
  virtual ~DataGroup();


  /*!
   *
   * @param rhs the DataObject to be copied
   * @return *this
   */
  DataGroup& operator=( const DataGroup& rhs );

  /*!
   *
   * @param rhs the DataObject to be moved into *this
   * @return *this
   */
  DataGroup& operator=( const DataGroup&& rhs );

  ///@}




  /**
   *
   * @name creation, allocation, insertion
   *
   */
  ///@{

  DataObject* CreateDataObject( const std::string& name );
  DataObject* DetachDataObject( const std::string& name );

  void RemoveDataObject( const std::string& name, const bool removeFromDataStore );

  void ReleaseDataObject( const std::string& name );

  DataGroup* CreateDataGroup( const std::string& name );

//  virtual DataGroup* Allocate() override

  ///@}






  /**
   * @name access functions
   */
  ///@{

  /*
  template< typename T >
  DataObject* SetParameter( const std::string& name , const T& value );

  template< typename T >
  const T& GetParameter( const std::string& name ) const
  {
    return (*(GetData<T>( name )));
  }

  template< typename T >
  T& GetParameter( const std::string& name )
  {
    return (*(GetData<T*>( name )));
  }
*/
  std::string& Name()  {return m_name;}
  const std::string& Name() const {return m_name;}


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


  DataGroup* GetDataGroup( const std::string& name )
  {
    return m_childGroups.at(name);
  }

  const DataGroup* GetDataGroup( const std::string& name ) const
  {
    return m_childGroups.at(name);
  }


  DataGroup* GetDataGroup( const Attribute& Attribute );

  DataObject* GetDataObject(const std::string& name )
  {
    const std::size_t lookup = m_DataObjectLookup.at(name);
    return m_DataObjects[ lookup ];
  }

  DataObject const * GetDataObject(const std::string& name ) const
  {
    const std::size_t lookup = m_DataObjectLookup.at(name);
    return m_DataObjects[ lookup ];
  }

  const dataArrayType GetDataObjects() const
  {
    return m_DataObjects;
  }

  template< typename DATATYPE >
  typename std::enable_if<std::is_pointer<DATATYPE>::value,DATATYPE>::type GetData( const std::string& name )
  {
    const std::size_t lookup = m_DataObjectLookup.at(name);
	  return m_DataObjects[ lookup ]->GetData<DATATYPE>();
  }

  template< typename DATATYPE >
  typename std::enable_if<std::is_pointer<DATATYPE>::value,const DATATYPE>::type GetData( const std::string& name ) const
  {
    const std::size_t lookup = m_DataObjectLookup.at(name);
    return m_DataObjects[ lookup ]->GetData<DATATYPE>();
  }

  const DataShape& GetDataShape() const
  {
    return m_dataShape;
  }

  const DataShape& GetDataShape( const std::string& name ) const
  {
    return GetDataObject(name)->GetDataShape();
  }

  ///@}

  std::size_t Lookup( const std::string& name )
  {
    return m_DataObjectLookup.at(name);
  }

};

/*
template< typename T >
DataObject* DataGroup::SetParameter( const std::string& name , const T& value )
{
  DataObject* const newParam = this->CreateDataObject( name );
  const std::size_t one = 1;
  DataShape desc(1,&one );
  newParam->SetType<T>();
  newParam->SetDataShape(desc)->Allocate();
  *(newParam->GetData<T*>()) = value;
  return newParam;
}
*/

} /* namespace DataStore */
#endif /* DATAGROUP_HPP_ */

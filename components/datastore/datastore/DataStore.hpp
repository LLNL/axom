/**
 * @name DataStore.hpp
 *
 *  @date Dec 2, 2014
 *  @author settgast
 */

#ifndef DATASTORE_HPP_
#define DATASTPRE_HPP_

#include "DataObject.hpp"
#include <memory>
#include "Utilities.hpp"
#include <vector>
#include "Types.hpp"

namespace DataStoreNS
{


/**
 * \class DataStore
 *
 * \brief class to manage collections of DataObjects, and attribute DataGroups
 */
class DataStore
{
public:
  // change to regular pointer for now.
  typedef std::vector< DataObject* > dataObjectContainerType;
  typedef std::unordered_map<std::string, std::size_t> lookupType;

  typedef std::vector< DataGroup* > dataGroupContainerType;


private:

  /// container of dataObject*
  dataObjectContainerType m_DataObjects;

  /// hash table for index lookup of m_DataObjects by name
//  lookupType m_DataObjectLookup;

  std::unordered_map<std::string, DataGroup* > m_childGroups;
  lookupType m_DataGroupLookup;


public:
  /**
   * @name Constructor, Destructor, Assignment, Move, Copy
   */
  ///@{

  /*!
   *
   * @param name name of the group
   * @param path path
   * \brief constructor
   */
  DataStore();


  /*!
   * @param source
   * \brief default copy constructor
   */
  DataStore( const DataStore& source );

  /*!
   * @param source
   * \brief default move constructor
   */
  DataStore( DataStore&& source );

  /*!
   * \brief default destructor
   */
  virtual ~DataStore();


  /*!
   *
   * @param rhs the DataObject to be copied
   * @return *this
   */
  DataStore& operator=( const DataStore& rhs );

  /*!
   *
   * @param rhs the DataObject to be moved into *this
   * @return *this
   */
  DataStore& operator=( const DataStore&& rhs );

  ///@}




  /**
   *
   * @name creation, allocation, insertion
   *
   */
  ///@{

  DataObject* CreateDataObject( const std::string& name,
                                const DataGroup* const parent  );

  DataGroup* CreateDataGroup( const std::string& name );

//  virtual DataGroup* Allocate() override


  void RemoveDataObject( const std::size_t& index );
  void RemoveDataObject( DataObject*& dobj );




  ///@}






  /**
   * @name access functions
   */
  ///@{


  /*!
   *
   * @param dataDescriptor
   * @return
   */


  DataGroup* GetDataGroup( const std::string& name )
  {
    return m_childGroups.at(name);
  }

  const DataGroup* GetDataGroup( const std::string& name ) const
  {
    return m_childGroups.at(name);
  }

/*
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
  */

  const dataObjectContainerType& GetDataObjects() const
  {
    return m_DataObjects;
  }

  /*
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



  ///@}

  std::size_t Lookup( const std::string& name )
  {
    return m_DataObjectLookup.at(name);
  }
*/
};


} /* namespace DataStoreNS */
#endif /* DATASTORE_HPP_ */

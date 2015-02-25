/*
 * Dataset.hpp
 *
 *  Created on: Dec 1, 2014
 *      Author: settgast
 */

#ifndef DATASET_HPP_
#define DATASET_HPP_

#include<map>
#include<unordered_map>
#include "Attributes.hpp"
#include<vector>
#include<typeindex>

namespace DataStoreNS
{

class DataGroup;
class DataStore;

/**
 * \class DataObject
 *
 * \brief A class to manages interface to an actual piece of data, whether it be owned and allocated by the DataObject
 * or passed in from the outside.
 *
 * Requirements for this class are:
 *    - one-to-one relation with data. This means that a DataObject is the only DataObject that refers to a specific
 *      data address.
 *    - contain pointer to the owning entity of this DataObject.
 *    - contain collection of attributes
 *    - description of shape of data (if applicable)
 *    - description of type of data (if applicable)
 *
 */
class DataObject
{
public:

  /**
   * @name Constructor, Destructor, Assignment, Move, Copy
   */
  ///@{

  /// non-callable default constructor
  DataObject() = delete;

  /**
   *
   * @param name name of the object
   * @param path path
   * \brief constructor
   */
  DataObject( const std::string& name,
              const DataGroup* const parent );

  /**
   * @param source
   * \brief default copy constructor
   */
  DataObject( const DataObject& source );


  /**
   * @param source
   * \brief default move constructor
   */
  DataObject( DataObject&& source );

  /**
   * \brief default destructor
   */
  virtual ~DataObject();

  /**
   *
   * @param rhs the DataObject to be copied
   * @return *this
   * \brief copy assignment
   */
  DataObject& operator=( const DataObject& rhs );

  /**
   *
   * @param rhs the DataObject to be moved into *this
   * @return *this
   * \brief move assignment
   */
  DataObject& operator=( const DataObject&& rhs );

  ///@}





  /**
   *
   * @name creation, allocation, insertion
   *
   */
  ///@{

  /**
   *
   * @return
   */
  virtual DataObject* Allocate();

  /**
   *
   * @param pathandName
   * @param data
   */
  template< typename TYPE >
  void insert( const std::string& pathandName, TYPE& data )
  {}


  ///@}





  /**
   * @name access functions
   */
  ///@{

  /**
   * @param name the attribute name
   * @return the attribute value
   */
  /*
  const Attribute& GetAttribute( const std::string& name ) const
  {
    return m_attributes.at(name);
  }
  */

  /**
   * @param newAttribute the attribute that is to be added
   * @return *this
   */
//  DataObject* SetAttribute( const Attribute& newAttribute );

  std::map< std::string, DataGroup* >& GetGroups()
  {
    return m_groups;
  }




  /**
   *
   * @return m_dataShape
   */
  const DataShape& GetDataShape() const
  { return m_dataShape; }

  /**
   *
   * @param dataShape
   * @return
   */
  virtual DataObject* SetDataShape( const DataShape& dataShape )
  {
    // check to see what conditions m_dataDescriptor can be set.
    m_dataShape = dataShape;
    m_data = dataShape.m_dataPtr;
    return this;
  }

  /**
   * @param ptr
   * @return
   */
  template< typename T >
  typename std::enable_if<std::is_pointer<T>::value,DataObject*>::type SetDataPointer( T const ptr )
  {
    m_data = static_cast<void*>(ptr);
    return this;
  }


  rtTypes::TypeID GetType() const
  {
    return m_dataType;
  }

  DataObject* SetType( const rtTypes::TypeID type )
  {
    m_dataType = type;
    return this;
  }

  template< typename T >
  DataObject* SetType()
  {
    m_dataType = rtTypes::GetTypeID<T>();
    return this;
  }




  virtual DataObject* SetLength( const std::size_t newsize );

  std::size_t length()
  { return m_dataShape.m_dimensions[0]; }

  /**
   *
   * @return m_name
   */
  std::string& Name()  {return m_name;}
  const std::string& Name() const {return m_name;}


  /**
   *
   * @param attributeKey
   * @return
   */
  bool HasAttribute( const std::string& attributeKey ) const {}

  /**
   *
   * @param attributeKey
   * @return
   */
  bool DeleteAttribute( const std::string& attributeKey );


  /**
   *
   * @return casted pointer to m_data
   * \brief there is a lot more that has to happen to ensure that the cast is legal
   */
  template< typename TYPE >
  typename std::enable_if<std::is_pointer<TYPE>::value,TYPE>::type GetData()
  { return static_cast<TYPE>(m_data); }


  const DataGroup* GetParent() const
  { return m_parent; }

  ///@}














private:
  /// pointer to DataGroup that "owns" this DataObject
  const DataGroup* m_parent;

  /// unique name of this DataObject
  std::string m_name;

  /// hash table of attributes keyed by the attribute name.
//  std::unordered_map<std::string,Attribute> m_attributes;

  /// pointer to the data. This is intended to be a one-to-one relationship (i.e. One and only one DataObjects m_data are equivalent to this->m_data.
  void* m_data;

  /// a complete description of the data, assuming one exists. Mainly for simple data arrays, and potentially
  /// pod structure.
  DataShape m_dataShape;

  ///
  //  std::type_index m_dataType;
  rtTypes::TypeID m_dataType;

  /// list of groups/attributes that hold pointers to "this"
  std::map< std::string, DataGroup* > m_groups;


  /// use a vector to allocate data until we implement an appropriate allocator interface.
  std::vector<char> m_memblob;

};




} /* namespace DataStore */
#endif /* DATASET_HPP_ */

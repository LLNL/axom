/*
 * Dataset.hpp
 *
 *  Created on: Dec 1, 2014
 *      Author: settgast
 */

#ifndef DATASET_HPP_
#define DATASET_HPP_

#include <map>
#include <set>
#include <vector>
#include "Types.hpp"

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
  /// brief container of DataGroup pointers.
  typedef std::set< DataGroup* > GroupContainerType;

private:


  /// universal identification - unique within a DataStore
  IDType m_uid;

  /// a string that describes what is in the DataObject
  std::string m_stringDescriptor;

  /// container of groups that contain this DataObject
  GroupContainerType m_GroupSet;

  /// pointer to the data. This is intended to be a one-to-one relationship (i.e. One and only one DataObjects m_data are equivalent to this->m_data.
  void* m_data;

  /// a complete description of the data, assuming one exists. Mainly for simple data arrays, and potentially
  /// pod structure.
  DataShape m_dataShape;

  ///
  rtTypes::TypeID m_dataType;

  /// use a vector to allocate data until we implement an appropriate allocator interface.
  std::vector<char> m_memblob;

public:
  /// sample attribute, would need a more generic system
  bool m_dump;

  /*!
   *
   * @param uid
   * @param m_stringDescriptor
   */
  DataObject( const IDType uid,
              const std::string& m_stringDescriptor );

  /*!
   *
   * @param uid
   */
  DataObject( const IDType uid );

  /*!
   *
   * @param source
   */
  DataObject(const DataObject& source );


  /*!
   * destructor
   */
  ~DataObject();

  /*!
   * \brief Return the univeral id for this DataObject.
   */
  IDType GetUID() { return m_uid; }

  /*!
   * @param grp Pointer to DataGroup to associate with this DataObject.
   * \brief Associate a DataGroup with this DataObject.
   */
  void AttachGroup( DataGroup *grp ) { m_GroupSet.insert(grp); }

  /*!
   * @param grp Pointer to DataGroup to disassociate with this DataObject.
   * \brief Disassociate a DataGroup with this DataObject.
   */
  void DetachGroup( DataGroup *grp ) { m_GroupSet.erase(grp); }

  /*!
   * @param grp Pointer to DataGroup to test for membership.
   * \brief return true is grp is attached.
   */
  bool IsAttachedGroup( DataGroup *grp )
  {
    GroupContainerType::iterator it = m_GroupSet.find(grp);
    if (it == m_GroupSet.end())
      return false;
    else
      return true;
  }

  /*!
   * \brief Return DataGroups attached to this DataObject.
   */
  GroupContainerType *GetDataGroups() { return &m_GroupSet; }

  /**
   *
   * @return casted pointer to m_data
   * \brief there is a lot more that has to happen to ensure that the cast is legal
   */
  template< typename TYPE >
#ifdef USECXX11
  typename std::enable_if<std::is_pointer<TYPE>::value,TYPE>::type
#else
  TYPE
#endif
  GetData()
  { return static_cast<TYPE>(m_data); }

  template< typename TYPE >
#ifdef USECXX11
  typename std::enable_if<std::is_pointer<TYPE>::value,TYPE>::type
#else
  const TYPE
#endif
  GetData() const
  { return static_cast<const TYPE>(m_data); }




  /**
   * @name members that will be deprecated by conduit
   */
  ///@{

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
   *
   * @return
   */
  DataObject* Allocate();

  DataObject* SetLength(const std::size_t newsize);


  ///@}


};




} /* namespace DataStore */
#endif /* DATASET_HPP_ */

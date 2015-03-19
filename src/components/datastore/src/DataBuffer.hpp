/**
 * DataBuffer.hpp
 *
 */

#ifndef DATABUFFER_HPP_
#define DATABUFFER_HPP_

#include <map>
#include <set>
#include <vector>
#include "Types.hpp"

namespace DataStoreNS
{

class DataGroup;
class DataStore;
class DataView;

/**
 * \class DataBuffer
 *
 * \brief A class to manages interface to an actual piece of data, whether it be owned and allocated by the DataBuffer
 * or passed in from the outside.
 *
 * Requirements for this class are:
 *    - one-to-one relation with data. This means that a DataBuffer is the only DataBuffer that refers to a specific
 *      data address.
 *    - contain pointer to the owning entity of this DataBuffer.
 *    - contain collection of attributes
 *    - description of shape of data (if applicable)
 *    - description of type of data (if applicable)
 *
 */
class DataBuffer
{
public:
  /// brief container of DataGroup pointers.
  typedef std::set< DataView* > ViewContainerType;

private:


  /// universal identification - unique within a DataStore
  IDType m_uid;

  /// a string that describes what is in the DataBuffer
  std::string m_stringDescriptor;

  /// container of groups that contain this DataBuffer
  ViewContainerType m_ViewSet;

  /// pointer to the data. This is intended to be a one-to-one relationship (i.e. One and only one DataBuffers m_data are equivalent to this->m_data.
  void* m_data;

  DataShape m_dataShape;

  ///
  rtTypes::TypeID m_dataType;

  /// use a vector to allocate data until we implement an appropriate allocator interface.
  std::vector<char> m_memblob;

public:

  /*!
   *
   * @param uid
   * @param m_stringDescriptor
   */
  DataBuffer( const IDType uid,
              const std::string& m_stringDescriptor );

  /*!
   *
   * @param uid
   */
  DataBuffer( const IDType uid );

  /*!
   *
   * @param source
   */
  DataBuffer(const DataBuffer& source );


  /*!
   * destructor
   */
  ~DataBuffer();

  /*!
   * \brief Return the univeral id for this DataBuffer.
   */
  IDType GetUID() { return m_uid; }

  /*!
   * @param grp Pointer to DataGroup to associate with this DataBuffer.
   * \brief Associate a DataGroup with this DataBuffer.
   */
  void AttachGroup( DataView *grp ) { m_ViewSet.insert(grp); }

  /*!
   * @param grp Pointer to DataGroup to disassociate with this DataBuffer.
   * \brief Disassociate a DataGroup with this DataBuffer.
   */
  void DetachGroup( DataView *grp ) { m_ViewSet.erase(grp); }

  /*!
   * @param grp Pointer to DataGroup to test for membership.
   * \brief return true is grp is attached.
   */
  bool IsAttachedGroup( DataView *grp )
  {
    ViewContainerType::iterator it = m_ViewSet.find(grp);
    if (it == m_ViewSet.end())
      return false;
    else
      return true;
  }

  /*!
   * \brief Return DataGroups attached to this DataBuffer.
   */
  ViewContainerType *GetDataViews() { return &m_ViewSet; }

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

  DataBuffer* SetType( const rtTypes::TypeID type )
  {
    m_dataType = type;
    return this;
  }

  template< typename T >
  DataBuffer* SetType()
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
  virtual DataBuffer* SetDataShape( const DataShape& dataShape )
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
  DataBuffer* Allocate();

  ///@}


};




} /* namespace DataStore */
#endif /* DATABUFFER_HPP_ */

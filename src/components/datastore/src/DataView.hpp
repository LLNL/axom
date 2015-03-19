/**
 * DataView.hpp
 *
 */

#ifndef DATAVIEW_HPP_
#define DATAVIEW_HPP_

#include <map>
#include <set>
#include <vector>
#include "Types.hpp"

namespace DataStoreNS
{

class DataGroup;
class DataStore;
class DataBuffer;
/**
 * \class
 *
 * \brief A class to manages interface to an actual piece of data, whether it be owned and allocated by the
 * or passed in from the outside.
 *
 * Requirements for this class are:
 *    - one-to-one relation with data. This means that a  is the only  that refers to a specific
 *      data address.
 *    - contain pointer to the owning entity of this .
 *    - contain collection of attributes
 *    - description of shape of data (if applicable)
 *    - description of type of data (if applicable)
 *
 */
class DataView
{
public:
  /// brief container of DataGroup pointers.
  typedef std::set< DataGroup* > GroupContainerType;

private:

  /// Index of the dataBuffer in the DataStore
  IDType m_uid;

  /// pointer to the BEGINNING of the DataBuffer??
  DataBuffer* m_dataBuffer;

  /// Set of groups that contain this view??? Not necessary anymore, as the group/view relationship may be 1-to-1
  GroupContainerType m_GroupSet;

  /// pointer to the location in the actual data buffer that the view has access to
  void* m_viewStart;

  /// a complete description of the data, assuming one exists. Mainly for simple data arrays, and potentially
  /// pod structure.
  DataShape m_dataShape;

  ///
  rtTypes::TypeID m_dataType;


public:

  /*!
   *
   * @param uid
   * @param m_stringDescriptor
   */
  DataView( const IDType uid,
              const std::string& m_stringDescriptor );

  /*!
   *
   * @param uid
   */
  DataView( const IDType uid );

  /*!
   *
   * @param source
   */
  DataView(const DataView& source );


  /*!
   * destructor
   */
  ~DataView();

  /*!
   * \brief Return the univeral id for this DataView.
   */
  IDType GetUID() { return m_uid; }

#if 0
  /*!
   * @param grp Pointer to DataGroup to associate with this DataView.
   * \brief Associate a DataGroup with this DataView.
   */
  void AttachGroup( DataGroup *grp ) { m_GroupSet.insert(grp); }

  /*!
   * @param grp Pointer to DataGroup to disassociate with this DataView.
   * \brief Disassociate a DataGroup with this DataView.
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
   * \brief Return DataGroups attached to this DataView.
   */
  GroupContainerType *GetDataGroups() { return &m_GroupSet; }
#endif
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
  { return static_cast<TYPE>(m_viewStart); }

  template< typename TYPE >
#ifdef USECXX11
  typename std::enable_if<std::is_pointer<TYPE>::value,TYPE>::type
#else
  const TYPE
#endif
  GetData() const
  { return static_cast<const TYPE>(m_viewStart); }




  /**
   * @name members that will be deprecated by conduit
   */
  ///@{

  rtTypes::TypeID GetType() const
  {
    return m_dataType;
  }

  DataView* SetType( const rtTypes::TypeID type )
  {
    m_dataType = type;
    return this;
  }

  template< typename T >
  DataView* SetType()
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
  virtual DataView* SetDataShape( const DataShape& dataShape )
  {
    // check to see what conditions m_dataDescriptor can be set.
    m_dataShape = dataShape;
    m_viewStart = dataShape.m_dataPtr;
    return this;
  }



  /**
   *
   * @return
   */
  DataView* Allocate();

  DataView* SetLength(const std::size_t newsize);


  ///@}


};




} /* namespace DataStore */
#endif /* DATAVIEW_HPP_ */

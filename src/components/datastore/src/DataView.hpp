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

  std::string m_name;

  DataGroup* m_parentGroup;

  /// pointer to the DataBuffer??
  DataBuffer* m_dataBuffer;

  /// pointer to the location in the actual data buffer that the view has access to
  void* m_viewStart;

  /// a complete description of the data, assuming one exists. Mainly for simple data arrays, and potentially
  /// pod structure.
  DataShape m_dataShape;

  ///
  rtTypes::TypeID m_dataType;


public:


  DataView( const std::string& name,
            DataGroup* const parentGroup,
            DataBuffer* const dataBuffer );

  DataView( const std::string& name,
            DataGroup* const parentGroup,
            DataStore* const dataStore );


  /// copy constructor
  DataView(const DataView& source );


  /// destructor
  ~DataView();


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

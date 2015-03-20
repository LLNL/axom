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

#include "conduit/conduit.h"

using conduit::Node;
using conduit::Schema;
using conduit::DataType;

namespace DataStoreNS
{

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
  /// container of DataView pointers.
  typedef std::set< DataView* > ViewContainerType;

  /// default constructor
  DataBuffer();

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



  void AddDataView( DataView* dataView );
  void RemoveDataView( DataView* dataView );

  void ReconcileDataViews();

  /*!
   * \brief Return the universal id for this DataBuffer.
   */
  IDType GetUID() { return m_uid; }

  /*!
   * \brief Return DataGroups attached to this DataBuffer.
   */
  ViewContainerType *GetDataViews() { return &m_ViewContainer; }


  void *GetData()
  { return m_data;}


  DataBuffer* SetDescriptor(const Schema &schema)
  {
    m_schema.set(schema);
    return this;
  }
  
  
  DataBuffer* SetDescriptor(const DataType &dtype)
  {
    m_schema.set(dtype);
    return this;
  }
  
  
  /**
   *
   * @return m_schema
   */
  const Schema &GetDescriptor() const
  { return m_schema; }

  /// TODO: dangerous const issue needs to be resolved!
  Node &GetNode()
  { return m_node; }  


  DataBuffer *ApplyDescriptor();


  /**
   *
   * @return
   */
  DataBuffer* Allocate();

  ///@}

private:


  /// universal identification - unique within a DataStore
  IDType m_uid;

  /// a string that describes what is in the DataBuffer
  std::string m_stringDescriptor;

  /// container of groups that contain this DataBuffer
  ViewContainerType m_ViewContainer;

  /// pointer to the data. This is intended to be a one-to-one relationship (i.e. One and only one DataBuffers m_data are equivalent to this->m_data.
  void* m_data;
  
  /// use a vector to allocate data until we implement an appropriate allocator interface.
  std::vector<char> m_memblob;

  Node   m_node;
  Schema m_schema;

};




} /* namespace DataStore */
#endif /* DATABUFFER_HPP_ */

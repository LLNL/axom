/**
 * @name DataStore.hpp
 *
 *  @date Dec 2, 2014
 *  @author settgast
 */

#ifndef DATASTORE_HPP_
#define DATASTORE_HPP_

#include <vector>
#include <stack>
#include "DataGroup.hpp"
#include "Types.hpp"

namespace DataStoreNS
{

class DataBuffer;

/**
 * \class DataStore
 *
 * \brief Class to own DataObject and a root DataGroup.
 */
class DataStore
{
public:

  /// typedef for container for DataBuffers
  typedef std::vector< DataBuffer* > dataBufferContainerType;

  /*!
   * \brief Constructor.
   */
  DataStore();

  /*!
   * \brief Destructor.
   */
  ~DataStore();

  /*!
   * \brief Create a DataBuffer.
   *    It is assigned a universal id and owned by the DataStore
   */
  DataBuffer* CreateBuffer();


  /*!
   * @param id  Universal id of the DataObject.
   * \brief Remove a DataObject from the DataStore.
   *   It is disassociated with all groups and returned to the free pool.
   */
  void DestroyBuffer( const IDType id );

  /*!
   *
   * @param id the index
   * @return the DataBuffer that was at m_DataBuffer[id]
   * \brief Remove a DataBuffer from container, and return
   */
  DataBuffer* DetatchBuffer( const IDType id );


  /*!
   * \brief Return pointer to the root DataGroup.
   */
  DataGroup* GetRoot() 
      { return &m_RootGroup; };

  /*!
   *
   * @param id
   * @return
   */
  DataBuffer* GetBuffer( const IDType id ) 
      { return m_DataBuffers[id]; }

  void DestroyBuffers();

  void Print() const;
  void Print(Node &) const;

private:


  /// Root data group, created automatically with datastore.
  DataGroup m_RootGroup;

  /// container of DataBuffers
  dataBufferContainerType m_DataBuffers;

  /// stack of unique ids that can be recycled
  std::stack< IDType > m_AvailableDataBuffers;

#ifdef USECXX11
  DataStore( const DataStore& ) = delete;
  DataStore( DataStore&& ) = delete;

  DataStore& operator=( const DataStore& ) = delete;
  DataStore& operator=( DataStore&& ) = delete;
#else
  DataStore( const DataStore& );
  DataStore& operator=( const DataStore& );
#endif
};



} /* namespace DataStoreNS */
#endif /* DATASTORE_HPP_ */

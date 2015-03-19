/**
 * @name DataStore.hpp
 *
 *  @date Dec 2, 2014
 *  @author settgast
 */

#ifndef DATASTORE_HPP_
#define DATASTORE_HPP_

#include "DataGroup.hpp"
#include "DataBuffer.hpp"
#include <vector>
#include <stack>
#include "Types.hpp"

namespace DataStoreNS
{


/**
 * \class DataStore
 *
 * \brief Class to own DataObject and a root DataGroup.
 */
class DataStore
{
public:
  /*!
   * \brief vector of DataObject pointers.
   */
  typedef std::vector< DataBuffer* > dataObjectContainerType;

  // constructor
  // creates empty root data group and names it "/"
  // calls reserve on data object vector to minimize later resizes
private:


  // Root data group, created automatically with datastore.
  DataGroup m_RootGroup;

  // container of data object pointers
  // as long as we recycle ids, this should not have many vacancies with NULL pointers
  // if it's an issue, change this to a std::map (or boost/std::unordered_map)
  dataObjectContainerType m_DataBuffers;

  // stack of unique ids that can be recycled
  std::stack< IDType > m_AvailableDataBuffers;

public:
  /*!
   * \brief Constructor.
   */
  DataStore() :
    m_RootGroup(nullptr, this) {};

  /*!
   * \brief Destructor.
   */
  ~DataStore();

  // copy constructor
  //DataStore( const DataStore* store );

  /*!
   * \brief Create an empty DataObject.
   *    It is assigned a universal id and owned by the DataStore
   */
  DataBuffer *CreateDataBuffer();

  /*!
   * @param obj DataObject to delete.
   * \brief Remove a DataObject from the DataStore.
   *   It is disassociated with all groups and returned to the free pool.
   */
  void DeleteDataBuffer( DataBuffer*& obj);

  /*!
   * @param id  Universal id of the DataObject.
   * \brief Remove a DataObject from the DataStore.
   *   It is disassociated with all groups and returned to the free pool.
   */
  void DeleteDataBuffer(IDType id);

  /*!
   * \brief Return pointer to the root DataGroup.
   */
  DataGroup* GetRootDataGroup() { return &m_RootGroup; };

};



} /* namespace DataStoreNS */
#endif /* DATASTORE_HPP_ */

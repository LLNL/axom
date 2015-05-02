/*!
 ******************************************************************************
 *
 * \file
 *
 * \brief   Header file containing definition of DataStore class.
 *
 ******************************************************************************
 */

#ifndef DATASTORE_HPP_
#define DATASTORE_HPP_

// Standard C++ headers
#include <vector>
#include <stack>

// SiDRe project headers
#include "Types.hpp"


// using directives to make Conduit usage easier and less visible
using conduit::Node;


namespace sidre
{

class DataBuffer;
class DataGroup;

/*!
 * \class DataStore
 *
 * \brief DataStore is the main interface for creating and accessing
 *        DataBuffer objects.
 *
 * It maintains a collection of DataBuffer objects and owns the "root" 
 * DataGroup, called "/". A DataGroup hierachy (a tree) is created by 
 * creating (sub) DataGroups withing the root group and (sub) DataGroups
 * within child DataGroups.
 *
 */
class DataStore
{
public:

  /*!
   * \brief Default ctor initializes datastore object and creates root group.
   */
  DataStore();

  /*!
   * \brief Dtor destroys all contents of the datastore, including data held
   *        in buffers if owned by the buffers.
   */
  ~DataStore();


//@{
//!  @name DataBuffer methods

  /*!
   * \brief Return (non-const) pointer to data buffer object with given id.
   */
  DataBuffer* getBuffer( const IDType id ) 
  { 
     return m_DataBuffers[id]; 
  }

  /*!
   * \brief Create a data buffer object and return a pointer to it.
   *
   *    The buffer object is assigned a unique id when created and the 
   *    buffer object is owned by the data store.
   */
  DataBuffer* createBuffer();

  /*!
   * \brief Remove data buffer with given id from the datastore and 
   *        destroy it (including its data if data buffer owns it).
   *
   *   Note that buffer destruction detaches it from all groups and views 
   *   it was associated with.
   */
  void destroyBuffer( const IDType id );

  /*!
   * \brief Remove all data buffers from the datastore and destroy them 
   *        (including data they own).
   *
   *   Note that buffer destruction detaches it from all groups and views
   *   it was associated with.
   */
  void destroyBuffers();

  /*!
   * \brief Remove data buffer with given id from the datastore, but leave
   *        it intact.
   *
   * \return pointer to DataBuffer object that was datached.
   */
  DataBuffer* detachBuffer( const IDType id );

  /*!
   * \brief Return number of buffers in the datastore.
   */
  IDType getNumberOfBuffers() const
  {
    return m_DataBuffers.size() - m_AvailableDataBuffers.size();
  }

//@}

  /*!
   * \brief Return pointer to the root DataGroup.
   */
  DataGroup* getRoot() 
  { 
      return m_RootGroup; 
  };


  /*!
   * \brief Copy buffer descriptions and group tree, starting at root,
   *        to given Conduit node. 
   */ 
  void info(Node& n) const;

  /*!
   * \brief Print JSON description of data buffers and group tree,
   *        starting at root, to stdout.
   */
  void print() const;


private:
  //
  // Unimplemented ctors and copy-assignment operators.
  //
#ifdef USECXX11
  DataStore( const DataStore& ) = delete;
  DataStore( DataStore&& ) = delete;

  DataStore& operator=( const DataStore& ) = delete;
  DataStore& operator=( DataStore&& ) = delete;
#else
  DataStore( const DataStore& );
  DataStore& operator=( const DataStore& );
#endif

  /// Root data group, created when datastore is created.
  DataGroup *m_RootGroup;

  /// container of DataBuffers
  std::vector<DataBuffer*> m_DataBuffers;

  /// stack of unique buffer ids that can be recycled
  std::stack< IDType > m_AvailableDataBuffers;

};



} /* namespace sidre */
#endif /* DATASTORE_HPP_ */

/*!
 ******************************************************************************
 *
 * \file
 *
 * \brief   Implementation file for DataStore class.
 *
 ******************************************************************************
 */ 

// Associated header file
#include "DataStore.hpp"

// SiDRe project headers
#include "DataBuffer.hpp"
#include "DataGroup.hpp"


namespace asctoolkit
{
namespace sidre
{


/*
*************************************************************************
*
* Datastore ctor creates root group.
*
*************************************************************************
*/
DataStore::DataStore()
{
    m_RootGroup = new DataGroup("/", this);
};


/*
*************************************************************************
*
* Datastore dtor destroys all contents.
*
*************************************************************************
*/
DataStore::~DataStore()
{
      // clean up groups and views before we destroy buffers
      delete m_RootGroup;
      destroyBuffers();
}
 

/*
*************************************************************************
*
* Create new data buffer and assign unique id.
*
*************************************************************************
*/
DataBuffer* DataStore::createBuffer()
{
  // TODO: implement pool, look for free nodes.  Allocate in blocks.
  common::IDType newIndex = m_DataBuffers.size();
  m_DataBuffers.push_back( nullptr );
  if( !m_AvailableDataBuffers.empty() )
  {
    newIndex = m_AvailableDataBuffers.top();
    m_AvailableDataBuffers.pop();
  }
  DataBuffer* const obj = new DataBuffer( newIndex );

  m_DataBuffers[newIndex] = obj ;

  return obj;
}
 

/*
*************************************************************************
*
* Remove data buffer with given id from the datastore and destroy it,
* recover its id for reuse.
*
*************************************************************************
*/
void DataStore::destroyBuffer( common::IDType id )
{
  delete m_DataBuffers[id];
  m_DataBuffers[id] = nullptr;
  m_AvailableDataBuffers.push(id);
}


/*
*************************************************************************
*
* Destroy all buffers in datastore.
*
* Should we recover the ids for potential reuse???
*
*************************************************************************
*/
void DataStore::destroyBuffers()
{
    for( std::vector<DataBuffer*>::iterator iter=m_DataBuffers.begin() ;
         iter!=m_DataBuffers.end() ; ++iter )
    {
      delete *iter;
    }
}


/*
*************************************************************************
*
* Remove data buffer with given id from the datastore leaving it intact,
* and return a pointer to it. Its id is recovered for reuse.
*
*************************************************************************
*/
DataBuffer* DataStore::detachBuffer( common::IDType id )
{
  DataBuffer* const rval = m_DataBuffers[id];
  m_DataBuffers[id] = nullptr;
  m_AvailableDataBuffers.push(id);

  return rval;
}

  
/*
*************************************************************************
*
* Copy buffer descriptions and group tree, starting at root, to given 
* Conduit node.
*
*************************************************************************
*/
void DataStore::info(Node &n) const
{
    m_RootGroup->info(n["DataStore/root"]);
    for( std::vector<DataBuffer*>::const_iterator iter=m_DataBuffers.begin() ;
         iter!=m_DataBuffers.end() ;
         ++iter )
    {
        Node &b = n["DataStore/buffers"].append();
        if(*iter != nullptr)
        {
            (*iter)->info(b);
        }
    }
}


/*
*************************************************************************
*
* Print JSON description of data buffers and group tree, starting at root, 
* to stdout.
*
*************************************************************************
*/
void DataStore::print() const
{
    Node n;
    info(n);
    n.print();
}


} /* end namespace sidre */
} /* end namespace asctoolkit */

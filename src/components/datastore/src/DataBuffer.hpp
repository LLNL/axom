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

namespace sidre
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
    friend class DataStore;
    friend class DataView;
    friend class DataGroup;

    /*!
     * \brief Return the universal id for this buffer.
     */
    IDType GetUID() const
    {
        return m_uid;
    }

    /*!
     * \brief Return a view attached to this buffer.
     */
    DataView *GetView(IDType idx)
    {
        return m_views[idx];
    }


    /*!
    * \brief Return number of views attached to this buffer.
    */
    size_t CountViews() const
    {
      return m_views.size();
    }
    

    void *GetData()
    { return m_data;}


    DataBuffer* Declare(const Schema &schema)
    {
        m_schema.set(schema);
        return this;
    }
    
    DataBuffer* Declare(const DataType &dtype)
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

    /// note: in most cases, we want to use the const version of the node
    Node &GetNode()
    {return m_node; }  

    const Node &GetNode() const
    { return m_node; }  

    DataBuffer* Allocate();
  
    DataBuffer* Allocate(const Schema &schema);
    DataBuffer* Allocate(const DataType &dtype);

    void Info(Node &n) const;
    void Print() const;
  ///@}

private:
    /// default constructor
    DataBuffer();

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

     /// TODO: the data store should be the only entity that can create a buffer

    void AttachView( DataView* dataView );
    void DetachView( DataView* dataView );



    /// universal identification - unique within a DataStore
    IDType m_uid;

    /// container of views that contain this DataBuffer
    std::vector<DataView *> m_views;

    /// pointer to the data. This is intended to be a one-to-one relationship (i.e. One and only one DataBuffers m_data are equivalent to this->m_data.
    void* m_data;
  
    /// use a vector to allocate data until we implement an appropriate allocator interface.
    std::vector<char> m_memblob;

    Node   m_node;
    Schema m_schema;

};




} /* namespace DataStore */
#endif /* DATABUFFER_HPP_ */

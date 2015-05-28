/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

/*!
 ******************************************************************************
 *
 * \file
 *
 * \brief   Header file containing definition of DataBuffer class.
 *
 ******************************************************************************
 */

#ifndef DATABUFFER_HPP_
#define DATABUFFER_HPP_

// Standard C++ headers
#include <map>
#include <set>
#include <vector>

// Other toolkit component headers
#include "conduit/conduit.h"

// SiDRe project headers
#include "SidreTypes.hpp"

// using directives to make Conduit usage easier and less visible
using conduit::DataType;
using conduit::Node;
using conduit::Schema;


namespace asctoolkit
{
namespace sidre
{

class DataStore;
class DataView;

/*!
 * \class DataBuffer
 *
 * \brief DataBuffer holds a data object, which it owns (and allocates!)
 *
 * The DataBuffer class has the following properties:
 *
 *    - DataBuffer objects can only be created via the DataStore interface,
 *      not directly. 
 *    - A DataBuffer object has a unique identifier within a DataStore,
 *      which is assigned by the DataStore when the buffer is created.
 *    - The data object owned by a DataBuffer is unique to that DataDuffer
 *      object; i.e.,  DataBuffers do not share data.
 *    - A DataBuffer object maintains a collection of DataViews that
 *      refer to its data.
 *
 */
class DataBuffer
{
public:

    //
    // Friend declarations to constrain usage via controlled access to
    // private members.
    //
    friend class DataStore;
    friend class DataGroup;
    friend class DataView;


    /*!
     * \brief Return the unique id of this buffer object.
     */
    IDType getUID() const
    {
        return m_uid;
    }

    /*!
     * \brief Return number of views attached to this buffer.
     */
    size_t getNumViews() const
    {
      return m_views.size();
    }
   
 
//@{
//!  @name Data declaration and allocation methods
 
    /*!
     * \brief Declare a data object as a Conduit schema.
     *
     * \return pointer to this DataBuffer object.
     */
    DataBuffer* declare(const Schema& schema)
    {
        m_schema.set(schema);
        return this;
    }
    
    /*!
     * \brief Declare a data object as a pre-defined Conduit data type.
     *
     * \return pointer to this DataBuffer object.
     */
    DataBuffer* declare(const DataType& dtype)
    {
        m_schema.set(dtype);
        return this;
    }

    /*!
     * \brief Allocate data previously declared using a Declare() method.
     *
     * \return pointer to this DataBuffer object.
     */
    DataBuffer* allocate();
  
    /*!
     * \brief Declare and allocate data described as a Conduit schema.
     *
     *        Equivalent to calling Declare(schema), then Allocate().
     *
     * \return pointer to this DataBuffer object.
     */
    DataBuffer* allocate(const Schema &schema);

    /*!
     * \brief Declare and allocate data described as a pre-defined 
     *        Conduit data type.
     *
     *        Equivalent to calling Declare(dtype), then Allocate().
     *
     * \return pointer to this DataBuffer object.
     */
    DataBuffer* allocate(const DataType& dtype);
  
    /*!
     * \brief Reallocate data described as a Conduit schema.
     *
     *        Leverages Conduit Node::update() semantics. 
     *
     * \return pointer to this DataBuffer object.
     */
    DataBuffer* reallocate(const Schema& schema);

    /*!
     * \brief Reallocate data described as a pre-defined 
     *        Conduit data type.
     *
     *        Leverages Conduit Node::update() semantics. 
     *
     * \return pointer to this DataBuffer object.
     */
    DataBuffer* reallocate(const DataType& dtype);


//@}


//@{
//!  @name Accessor methods
 
    /*!
     * \brief Return void-pointer to data held by DataBuffer.
     */
    void* getData()
    { 
       return m_data;
    }

    /*!
     * \brief Return non-const reference to Conduit node holding data.
     */
    Node& getNode()
    {
       return m_node; 
    }

    /*!
     * \brief Return const reference to Conduit node holding data.
     */
    const Node& getNode() const
    { 
       return m_node; 
    }

    /*!
     * \brief Return const reference to Conduit schema describing data.
     */
    const Schema& getSchema() const
    { 
       return m_schema; 
    }

    /*!
     * \brief Return pointer to view attached to this buffer identified
     *        by the given index.
     */
    DataView* getView(IDType idx)
    { 
       return m_views[idx]; 
    }

//@}


    /*!
     * \brief Copy data buffer description to given Conduit node.
     */
    void info(Node& n) const;

    /*!
     * \brief Print JSON description of data buffer to stdout.
     */
    void print() const;


private:

    /*!
     *  \brief Private ctor that assigns unique id.
     */
    DataBuffer( IDType uid );

    /*!
     * \brief Private copy ctor.
     */
    DataBuffer(const DataBuffer& source );

    /*!
     * \brief Private dtor.
     */
    ~DataBuffer();

    /*!
     * \brief Private methods to attach/detach data view to buffer.
     */
    void attachView( DataView* dataView );
    ///
    void detachView( DataView* dataView );

    
    

    /*!
     * \brief Private methods allocate and deallocate bytes data view to buffer.
     */
    void cleanup();
    ///
    void *allocateBytes(std::size_t numBytes);
    ///
    void  releaseBytes(void *);


    /// Identifier - unique within a dataStore.
    IDType m_uid;

    /// Container of DataViews attached to this buffer.
    std::vector<DataView *> m_views;

    /// Pointer to the data owned by DataBuffer.
    void* m_data;
  
    /// Conduit Node that holds buffer data.
    Node   m_node;

    /// Conduit Schema that describes buffer data.
    Schema m_schema;

    /*!
     *  Unimplemented ctors and copy-assignment operators.
     */
  #ifdef USE_CXX11
    DataBuffer() = delete;
    DataBuffer( DataBuffer&& ) = delete;

    DataBuffer& operator=( const DataBuffer& ) = delete;
    DataBuffer& operator=( DataBuffer&& ) = delete;
  #else
    DataBuffer();
    DataBuffer& operator=( const DataBuffer& );
  #endif

};


} /* end namespace sidre */
} /* end namespace asctoolkit */

#endif /* DATABUFFER_HPP_ */

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
 * \brief   Header file containing definition of DataView class.
 *
 ******************************************************************************
 */

#ifndef DATAVIEW_HPP_
#define DATAVIEW_HPP_

// Standard C++ headers
#include <map>
#include <set>
#include <string>
#include <vector>

// Other CS Toolkit headers
#include "common/CommonTypes.hpp"
#include "conduit/conduit.h"


// using directives to make Conduit usage easier and less visible
using conduit::Node;
using conduit::Schema;
using conduit::DataType;



namespace asctoolkit
{
namespace sidre
{

class DataBuffer;
class DataGroup;
class DataStore;

/*!
 * \class DataView
 *
 * \brief DataView provides a "view" into data, which may be owned by a
 *        DataBuffer object or owned externally.
 *
 * The DataView class has the following properties:
 *
 *    - DataView objects can only be created via the DataGroup interface,
 *      not directly. The view object is owned by the DataGroup object that
 *      creates it.
 *    - A DataView object has a unique name (string) within the DataGroup
 *      that owns it.
 *    - A DataView holds a pointer to the DataGroup that created it and which
 *      owns it.
 *    - A DataView object can refer to data in one of three ways:
 *        # It can decribe a view into data owned by a pre-existing DataBuffer.
 *          In this case, there can be multiple views into a single buffer.
 *        # It can be used to declare and allocate data, using semantics 
 *          similar to DataBuffer data declaration and allocation. When
 *          data is declared via a DataView object, there is a one-to-one
 *          relationship between the view and the buffer and the view owns
 *          the buffer. DataViews cannot share data in this case.  
 *        # It can hold a pointer to an "opaque" data object. In this case,
 *          there is no associated DataBuffer.
 *    - When a Dataview object is associated with a DataBuffer, the DataView
 *      object holds a pointer to the DataBuffer object. The data description 
 *      of the view is represented as a Conduit schema or a Conduit 
 *      pre-defined data type. 
 *
 */
class DataView
{
public:

    //
    // Friend declarations to constrain usage via controlled access to
    // private members.
    //
    friend class DataGroup;

//@{
//!  @name Data declaration and allocation methods

    ////////////////////////////////////////////////////////////////
    ///
    /// IMPORTANT: We need to clearly destiguish which of these 
    ///            methods apply to the case of a one-to-one 
    ///            buffer-view relationship and which apply to the
    ///            case of defining a view into an existing buffer 
    ///            object (possibly one of many). 
    ///
    ////////////////////////////////////////////////////////////////

    /*!
     * \brief Declare a data view as a Conduit schema.
     *
     * \return pointer to this DataView object.
     */
    DataView* declare(const Schema& schema);  

    /*!
     * \brief Declare a data view as a pre-defined Conduit data type.
     *
     * \return pointer to this DataView object.
     */
    DataView* declare(const DataType& dtype);
 
    /*!
     * \brief Allocate data for a DataView, previously declared.
     *
     * This is equivalent to calling the allocate() method on the
     * associated data buffer and then calling apply() on this DataView
     * object.
     *
     * \return pointer to this DataView object.
     */
    DataView* allocate();

    /*!
     * \brief Declare a data view as a Conduit schema and allocate the data.
     *
     * This is equivalent to calling declare(Schema), then allocate(),
     * and then calling apply() on this DataView object.  
     *
     * \return pointer to this DataView object.
     */
    DataView* allocate(const Schema& schema);

    /*!
     * \brief Declare a data view as a Conduit pre-defined data type
     *        and allocate the data.
     *
     * This is equivalent to calling declare(DataType), then allocate(),
     * and then calling apply() on this DataView object.
     *
     * \return pointer to this DataView object.
     */
    DataView* allocate(const DataType& dtype);
 
    /*!
     * \brief  Reallocate the views' underlying buffer using a Conduit 
     *         schema.
     *
     *  Uses conduit's update semantics();
     *
     * \return pointer to this DataView object.
     */
     DataView* reallocate(const Schema& schema);

    /*!
     * \brief  Reallocate the views' underlying buffer using a Conduit 
     *         data type.
     *
     *  Uses conduit's update semantics();
     *
     * \return pointer to this DataView object.
     */
     DataView* reallocate(const DataType& dtype);
 
   
    /*!
     * \brief Apply a previously declared data view to data held in
     *        the DataBuffer associated with this DataView object.
     *
     * \return pointer to this DataView object.
     */ 
    DataView* apply();

    /*!
     * \brief Declare a data view as a Conduit schema and apply the
     *        schema to the data in DataBuffer associated with this
     *        DataView object.
     *
     * This is equivalent to calling declare(Schema), then apply().
     *
     * \return pointer to this DataView object.
     */
    DataView* apply(const Schema& schema);

    /*!
     * \brief Declare a data object as a Conduit pre-defined data type
     *        and apply schema to the data in DataBuffer associated with 
     *        this DataView object
     *
     * This is equivalent to calling declare(DataType), then apply().
     *
     * \return pointer to this DataView object.
     */
    DataView* apply(const DataType& dtype);

//@}
  
  
//@{
//!  @name DataView query methods

    /*!
     * \brief Return true if DataView is associated with a DataBuffer
     *        (view is not opaque); false otherwise.
     */
    bool  hasBuffer() const
    { 
       return m_data_buffer != ATK_NULLPTR;
    }

    /*!
     * \brief Return true if DataView is opaque (view is not associated
     *        with a DataBuffer and view has no knowledge of data type, 
     *        structure, etc.); false otherwise.
     */
    bool  isOpaque() const
    {
       return m_is_opaque;
    }

    /*!
     * \brief Return true if data declaration has been applied to data 
     *        in DataBuffer associated with DataView; false otherwise.
     */
    bool isApplied() const
    {
       return m_is_applied;
    }

    /*!
     * \brief Return true if DataView describes externally-owned data (i.e.,
     *        not in an associated buffer); false otherwise.
     */
    bool isExternal() const
    {
       return m_is_external;
    }

//@}


//@{
//!  @name Accessor methods

    /*!
     * \brief Return const reference to name of DataView.
     */
    const std::string& getName() const
    {
       return m_name;
    }
    
    /*!
     * \brief Return void* pointer to data in view.
     */
    void* getOpaque() const;
 
 
    /*!
     * \brief Return pointer to non-const DataBuffer associated with DataView.
     */
    DataBuffer* getBuffer()
    { 
       return m_data_buffer; 
    }
 
    /*!
     * \brief Return pointer to const DataBuffer associated with DataView.
     */
    DataBuffer const* getBuffer() const
    { 
       return m_data_buffer; 
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
     * \brief Return pointer to non-const DataGroup that owns DataView object.
     */
    DataGroup* getOwningGroup()
    {
       return m_owning_group;
    }

    /*!
     * \brief Return pointer to const DataGroup that owns DataView object.
     */
    DataGroup const* getOwningGroup() const
    {
       return m_owning_group;
    }

//@}


    /*!
     * \brief Copy data view description to given Conduit node.
     */
    void info(Node& n) const;

    /*!
     * \brief Print JSON description of data view to stdout.
     */
    void print() const;


private:

    /*!
     *  \brief Private ctor that creates a DataView with given name
     *         in given parent group and which is associated with given
     *         DataBuffer object.
     */
    DataView( const std::string& name,
              DataGroup* const owning_group,
              DataBuffer* const data_buffer );

    /*!
     *  \brief Private ctor that creates an opaque DataView with given name
     *         in given parent group.
     */
    DataView( const std::string& name,
              DataGroup* const owning_group,
              void* opaque_ptr);

    /*!
     *  \brief Private ctor that creates an external DataView with given 
     *         name in given parent group and which is described by given
     *         Conduit data type.
     */
    DataView( const std::string& name,
              DataGroup* const owning_group,
              void* external_data,
              const DataType& dtype );

    /*!
     *  \brief Private ctor that creates an external DataView with given      
     *         name in given parent group and which is described by given
     *         Conduit schema.
     */
    DataView( const std::string& name,
              DataGroup* const owning_group,
              void* external_data,
              const Schema& schema );

    /*!
     * \brief Private copy ctor.
     */
    DataView(const DataView& source);

    /*!
     * \brief Private dtor.
     */ 
    ~DataView();


    /// Name of this DataView object.
    std::string m_name;

    /// DataGroup object that owns this DataView object.
    DataGroup*  m_owning_group;

    /// DataBuffer associated with this DataView object.
    DataBuffer* m_data_buffer;

    /// Conduit Schema that describes this DataView.
    Schema      m_schema;
    
    /// Conduit node used to access the data in this DataView.
    Node        m_node;
    
    /// Is this DataView opaque?
    bool        m_is_opaque;

    /// Is this a DataView into external data?
    bool        m_is_external;

    /// Has Schema been applied to buffer data?
    bool        m_is_applied;

    /*!
     *  Unimplemented ctors and copy-assignment operators.
     */
  #ifdef USE_CXX11
    DataView() = delete;
    DataView( DataView&& ) = delete;

    DataView& operator=( const DataView& ) = delete;
    DataView& operator=( DataView&& ) = delete;
  #else
    DataView();
    DataView& operator=( const DataView& );
  #endif

};


} /* end namespace sidre */
} /* end namespace asctoolkit */

#endif /* DATAVIEW_HPP_ */

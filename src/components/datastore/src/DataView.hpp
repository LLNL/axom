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

// SiDRe project headers
#include "Types.hpp"


// using directives to make Conduit usage easier and less visible
using conduit::Node;
using conduit::Schema;
using conduit::DataType;


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
    DataView* allocate(const DataType &dtype);
 
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
//!  @name Methods to query whether data is opaque or describes a buffer

    /*!
     * \brief Return true if DataView is associated with a DataBuffer
     *        (view is not opaque); false otherwise.
     */
    bool  hasBuffer() const
    { 
       return m_data_buffer != nullptr;
    }

    /*!
     * \brief Return true if DataView is opaque (view is not associated
     *        with a DataBuffer); false otherwise.
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
    const Schema& getDescriptor() const
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
     *  \brief Private ctor that creates a DataView with given name
     *         in given parent group and which is opaque.
     */
    DataView( const std::string& name,
              DataGroup* const owning_group,
              void* opaque_ptr);

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
    
    /// Has Schema been applied to buffer data?
    bool        m_is_applied;
    
    /// Is this DataView opaque?
    bool        m_is_opaque;
};


} /* namespace sidre */
#endif /* DATAVIEW_HPP_ */

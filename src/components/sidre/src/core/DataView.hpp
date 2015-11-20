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
#include <string>

// Other CS Toolkit headers
#include "common/CommonTypes.hpp"

// SiDRe project headers
#include "sidre/SidreTypes.hpp"


namespace asctoolkit
{
namespace sidre
{

// using directives to make Conduit usage easier and less visible
using conduit::Node;
using conduit::Schema;

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
 *      of the view may described by calling one of the apply() methods.  
 *
 */
class DataView
{
public:

  //
  // Friend declaration to constrain usage via controlled access to
  // private members.
  //
  friend class DataGroup;

//@{
//!  @name DataView declaration methods

// RDH TODO -- Remove declaration methods from the public interface, since 
//             we prefer to go through the datagroup interface?  
//
// RDH TODO -- Add offset, stride args to first method with defaults.
//

  /*!
   * \brief Declare a data view with given type and number of elements, and
   *        
   *
   * IMPORTANT: If view has been previously declared, this operation will
   *            re-declare the view. To have the new declaration take effect,
   *            the apply() method must be called.
   *
   * If given number of elements < 0, or view is opaque, method does nothing.
   *
   * \return pointer to this DataView object.
   */
  DataView * declare( TypeID type, SidreLength numelems);

  /*!
   * \brief Declare a data view with a Conduit data type object.
   *
   * IMPORTANT: If view has been previously declared, this operation will
   *            re-declare the view. To have the new declaration take effect,
   *            the apply() method must be called.
   *
   * If view is opaque, the method does nothing.
   *
   * \return pointer to this DataView object.
   */
  DataView * declare(const DataType& dtype);

  /*!
   * \brief Declare a data view with a Conduit schema object.
   *
   * IMPORTANT: If view has been previously declared, this operation will
   *            re-declare the view. To have the new declaration take effect,
   *            the apply() method must be called.
   *
   * If view is opaque, the method does nothing.
   *
   * \return pointer to this DataView object.
   */
  DataView * declare(const Schema& schema);

//@}


//@{
//!  @name DataView allocation methods

  /*!
   * \brief Allocate data for a view, previously declared.
   *
   * NOTE: Allocation from a view is only allowed if it is the only
   *       view associated with its buffer, or if it is not opaque.
   *       If either condition is true, this method does nothing.
   *
   * \return pointer to this DataView object.
   */
  DataView * allocate();

  /*!
   * \brief Allocate data for view with given type and number of elements.
   *
   * This is equivalent to calling: declare(type, numelems)->allocate();
   *
   * NOTE: Allocation from a view is only allowed if it is the only
   *       view associated with its buffer, or if it is not opaque.
   *       If either condition is true, or given number of elements is < 0,
   *       this method does nothing.
   *
   * \return pointer to this DataView object.
   */
  DataView * allocate( TypeID type, SidreLength numelems);

  /*!
   * \brief Allocate data for view described by a Conduit data type object.
   *
   * This is equivalent to calling: declare(dtype)->allocate().
   *
   * NOTE: Allocation from a view is only allowed if it is the only
   *       view associated with its buffer, or if it is not opaque.
   *       If either condition is true, this method does nothing.
   *
   * \return pointer to this DataView object.
   */
  DataView * allocate(const DataType& dtype);

  /*!
   * \brief Allocate data for view described by a Conduit schema object.
   *
   * This is equivalent to calling: declare(schema)->allocate().
   *
   * NOTE: Allocation from a view is only allowed if it is the only
   *       view associated with its buffer, or if it is not opaque.
   *       If either condition is true, this method does nothing.
   *
   * \return pointer to this DataView object.
   */
  DataView * allocate(const Schema& schema);

  /*!
   * \brief  Reallocate data for view to given number of elements (type
   *         stays the same).
   *
   * NOTE: Reallocation from a view is only allowed if it is the only
   *       view associated with its buffer, or if it is not opaque.
   *       If either condition is true, or given number of elements is < 0,
   *       this method does nothing.
   *
   * \return pointer to this DataView object.
   */
  DataView * reallocate(SidreLength numelems);

  /*!
   * \brief  Reallocate data for view as specified by Conduit data type object.
   *
   * NOTE: Reallocation from a view is only allowed if it is the only
   *       view associated with its buffer, or if it is not opaque.
   *       If either condition is true, this method does nothing.
   *
   * NOTE: The current type of the view data must match that of the given 
   *       data type object. If not, the method does nothing.
   *
   * \return pointer to this DataView object.
   */
  DataView * reallocate(const DataType& dtype);

  /*!
   * \brief  Reallocate data for view as specified by Conduit schema object.
   *
   * NOTE: Reallocation from a view is only allowed if it is the only
   *       view associated with its buffer, or if it is not opaque.
   *       If either condition is true, this method does nothing.
   *
   * NOTE: The current type of the view data must match that of the given 
   *       schema object. If not, the method does nothing.
   *
   * \return pointer to this DataView object.
   */
  DataView * reallocate(const Schema& schema);

//@}


//@{
//!  @name DataView (data description) apply methods

  /*!
   * \brief Apply data description of a previously declared to data view.
   *
   * If view is opaque, the method does nothing.
   *
   * \return pointer to this DataView object.
   */
  DataView * apply();

  /*!
   * \brief Apply data description defined by number of elements, and
   *        optionally offset and stride to data view (type remains the same).
   *
   * NOTE: The units for offset and stride are in number of elements, which
   *       is different than the DataType usage below where offset and stride
   *       are in number of bytes.
   *
   * IMPORTANT: If view has been previously declared (or applied), this 
   *            operation will apply the new data description to the view.
   *
   * IMPORTANT: If view has no data buffer object attached, this method 
   *            does nothing because it doesn't know type information.
   *
   * If given number of elements < 0, offset < 0, or view is opaque, the
   * method also does nothing.
   *
   * \return pointer to this DataView object.
   */
  DataView * apply( SidreLength numelems,
                    SidreLength offset = 0,
                    SidreLength stride = 1);

  /*!
   * \brief Apply data description defined by type and number of elements, and 
   *        optionally offset and stride to data view.
   *
   * NOTE: The units for offset and stride are in number of elements, which
   *       is different than the DataType usage below where offset and stride
   *       are in number of bytes.
   *
   * IMPORTANT: If view has been previously declared (or applied), this 
   *            operation will apply the new data description to the view.
   *
   * If given number of elements < 0, offset < 0, or view is opaque, the
   * method does nothing.
   *
   * \return pointer to this DataView object.
   */
  DataView * apply( TypeID type, SidreLength numelems,
                                 SidreLength offset = 0,
                                 SidreLength stride = 1);

  /*!
   * \brief Apply data description of given Conduit data type to data view.
   *
   * This is equivalent to calling: declare(dtype)->apply().
   *
   * If view is opaque, the method does nothing.
   *
   * \return pointer to this DataView object.
   */
  DataView * apply(const DataType& dtype);

  /*!
   * \brief Apply data description of given Conduit schema to data view.
   *        this DataView object
   *
   * This is equivalent to calling: declare(dtype)->apply().
   *
   * If view is opaque, the method does nothing.
   *
   * \return pointer to this DataView object.
   */
  DataView * apply(const Schema& schema);

//@}


//@{
//!  @name DataView query methods

  /*!
   * \brief Return true if view has a an associated DataBuffer object
   *        (e.g., view is not opaque, it has been created with a buffer,
   *         it has been allocated, etc.); false otherwise.
   */
  bool  hasBuffer() const
  {
    return m_data_buffer != ATK_NULLPTR;
  }

  /*!
   * \brief Return true if view is opaque (i.e., view has no knowledge 
   *        of the data it holds and can only return a void* pointer to it);
   *        false otherwise.
   */
  bool  isOpaque() const
  {
    return m_is_opaque;
  }

  /*!
   * \brief Return true if data declaration has been applied to data
   *        in buffer associated with view; false otherwise.
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
   * \brief Return void* pointer to opaque data in view, else ATK_NULLPTR
   *        if view has not been declared opaque.
   */
  void * getOpaque() const;

  /*!
   * \brief Return pointer to non-const DataBuffer associated with DataView.
   */
  DataBuffer * getBuffer()
  {
    return m_data_buffer;
  }

  /*!
   * \brief Return pointer to const DataBuffer associated with DataView.
   */
  const DataBuffer * getBuffer() const
  {
    return m_data_buffer;
  }

  /*!
   * \brief Return void-pointer to data associated with DataView.
   *
   * This will return the data pointer for all DataViews, including opaque
   * and external.
   */
  void * getDataPointer() const;

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
   * \brief Returns Value class instance that supports casting to the appropriate data return type.  This function
   * version does require enough type information for the compiler to know what to cast the Value class to.
   * Example:
   * int* myptr = getValue();
   * int myint = getValue();
   */
  Node::Value getValue()
  {
    return m_node.value();
  }

  /*!
   * \brief Lightweight templated wrapper around getValue that returns a Value class.  This function can be used in cases
   * were not enough information is provided to the compiler to cast the Value class based on the caller code line.  The
   * function template type must be explicitly provided on call.
   *
   * Example:
   * // will not work, compiler does not know what type to cast to for above getValue function.
   * assert( getValue() == 10 );
   * // use the templated version instead
   * assert (getValue<int>() == 10);
   */
  template<typename ValueType>
  ValueType getValue()
  {
    ValueType valueptr = m_node.value();
    return valueptr;
  }

  /*!
   * \brief Set value in conduit node.
   */
  template<typename ValueType>
  void setValue(ValueType value)
  {
    m_node.set(value);
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
  DataGroup * getOwningGroup()
  {
    return m_owning_group;
  }

  /*!
   * \brief Return pointer to const DataGroup that owns DataView object.
   */
  const DataGroup * getOwningGroup() const
  {
    return m_owning_group;
  }

  /*!
   * \brief Return type of data for this DataView object.
   */
  TypeID getTypeID() const
  {
    return static_cast<TypeID>(m_node.dtype().id());
  }

  /*!
   * \brief Return total number of bytes allocated by this DataView object.
   */
  size_t getTotalBytes() const
  {
    return m_schema.total_bytes();
  }

  /*!
   * \brief Return total number of elements allocated by this DataView object.
   */
  size_t getNumElements() const
  {
    return m_node.dtype().number_of_elements();
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

  /*!
   * \brief Print JSON description of data view to an ostream.
   */
  void print(std::ostream& os) const;

private:

  /*!
   *  \brief Private ctor that creates a DataView with given name
   *         in given parent group and which has no data associated with it.
   */
  DataView( const std::string& name,
            DataGroup * const owning_group );

  /*!
   *  \brief Private ctor that creates a DataView with given name
   *         in given parent group and which is associated with given
   *         DataBuffer object.
   */
  DataView( const std::string& name,
            DataGroup * const owning_group,
            DataBuffer * const data_buffer );

  /*!
   *  \brief Private ctor that creates an opaque DataView with given name
   *         in given parent group.
   */
  DataView( const std::string& name,
            DataGroup * const owning_group,
            void * opaque_ptr);

  /*!
   * \brief Private copy ctor.
   */
  DataView(const DataView& source);

  /*!
   * \brief Private dtor.
   */
  ~DataView();

  /*!
   *  \brief Private method returns true if data allocation on view is a
   *         valid operation; else false
   */
  bool allocationIsValid() const; 


  /// Name of this DataView object.
  std::string m_name;

  /// DataGroup object that owns this DataView object.
  DataGroup * m_owning_group;

  /// DataBuffer associated with this DataView object.
  DataBuffer * m_data_buffer;

  /// Conduit Schema that describes this DataView.
  Schema m_schema;

  /// Conduit node used to access the data in this DataView.
  Node m_node;

  /// Is this DataView opaque?
  bool m_is_opaque;

  /// Has Schema been applied data?
  bool m_is_applied;

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

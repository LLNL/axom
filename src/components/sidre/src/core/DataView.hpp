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
   * \brief Declare a data view from sidre type and length.
   *
   * If given length < 0 or view has previously been declared opaque,
   * the method does nothing.
   *
   * \return pointer to this DataView object.
   */
  DataView * declare( TypeID type, SidreLength len);

  /*!
   * \brief Declare a data view as a Conduit schema.
   *
   * If view has previously been declared opaque, the method does nothing.
   *
   * \return pointer to this DataView object.
   */
  DataView * declare(const Schema& schema);

  /*!
   * \brief Declare a data view as a pre-defined Conduit data type.
   *
   * If view has previously been declared opaque, the method does nothing.
   *
   * \return pointer to this DataView object.
   */
  DataView * declare(const DataType& dtype);

  /*!
   * \brief Allocate data for a DataView, previously declared.
   *
   * This is equivalent to calling the allocate() method on the
   * associated data buffer and then calling apply() on this DataView
   * object.
   *
   * NOTE: Allocation from a view only makes sense if this is the only
   *       view associated with its buffer. If this is not the case, the
   *       method does nothing.
   *
   * \return pointer to this DataView object.
   */
  DataView * allocate();

  /*!
   * \brief Declare a data view from a sidre type and length then allocate
   * the data.
   *
   * This is equivalent to calling declare(Schema), then allocate(),
   * and then calling apply() on this DataView object.
   *
   * NOTE: Allocation from a view only makes sense if this is the only
   *       view associated with its buffer. If this is not the case, the
   *       method does nothing.  Also, if the given length is < 0 or has been
   *       previously declared opaque, this method does nothing.
   *
   * \return pointer to this DataView object.
   */
  DataView * allocate( TypeID type, SidreLength len);

  /*!
   * \brief Declare a data view as a Conduit schema then allocate the data.
   *
   * This is equivalent to calling declare(Schema), then allocate(),
   * and then calling apply() on this DataView object.
   *
   * NOTE: Allocation from a view only makes sense if this is the only
   *       view associated with its buffer. If this is not the case, the
   *       method does nothing.
   *       If view has previously been declared opaque, the method does nothing.
   *
   * \return pointer to this DataView object.
   */
  DataView * allocate(const Schema& schema);

  /*!
   * \brief Declare a data view as a Conduit pre-defined data type
   *        then allocate the data.
   *
   * This is equivalent to calling declare(DataType), then allocate(),
   * and then calling apply() on this DataView object.
   *
   * NOTE: Allocation from a view only makes sense if this is the only
   *       view associated with its buffer. If this is not the case, the
   *       method does nothing.
   *       If view has previously been declared opaque, the method does nothing.
   *
   * \return pointer to this DataView object.
   */
  DataView * allocate(const DataType& dtype);

  /*!
   * \brief  Reallocate the view's underlying buffer using a Sidre type
   *         and length.
   *
   * This is equivalent to calling declare(TypeID), then allocate(),
   * and then calling apply() on this DataView object.
   *
   * NOTE: Re-allocation from a view only makes sense if this is the only
   *       view associated with its buffer. If this is not the case, the
   *       method does nothing.  Also, if the given length is < 0 or if
   *       this is an opaque view, this method does nothing.
   *
   * \return pointer to this DataView object.
   */
  DataView * reallocate(TypeID type, SidreLength len);

  /*!
   * \brief  Reallocate the view's underlying buffer using a Conduit
   *         schema.
   *
   * NOTE: Re-allocation from a view only makes sense if this is the only
   *       view associated with its buffer. If this is not the case, the
   *       method does nothing.  Also, if this is an opaque view, this
   *       method does nothing.
   *
   * \return pointer to this DataView object.
   */
  DataView * reallocate(const Schema& schema);

  /*!
   * \brief  Reallocate the view's underlying buffer using a Conduit
   *         data type.
   *
   * NOTE: Re-allocation from a view only makes sense if this is the only
   *       view associated with its buffer. If this is not the case, the
   *       method does nothing.  Also, if this is an opaque view, this
   *       method does nothing.
   *
   * \return pointer to this DataView object.
   */
  DataView * reallocate(const DataType& dtype);


  /*!
   * \brief Apply a previously declared data view to data held in
   *        the DataBuffer associated with this DataView object.
   *
   * If view has previously been declared opaque, the method does nothing.
   *
   * \return pointer to this DataView object.
   */
  DataView * apply();

  /*!
   * \brief Declare a data view as a Conduit schema and apply the
   *        schema to the data in DataBuffer associated with this
   *        DataView object.
   *
   * This is equivalent to calling declare(Schema), then apply().
   *
   * If view has previously been declared opaque, the method does nothing.
   *
   * \return pointer to this DataView object.
   */
  DataView * apply(const Schema& schema);

  /*!
   * \brief Declare a data object as a Conduit pre-defined data type
   *        and apply schema to the data in DataBuffer associated with
   *        this DataView object
   *
   * This is equivalent to calling declare(DataType), then apply().
   *
   * If view has previously been declared opaque, the method does nothing.
   *
   * \return pointer to this DataView object.
   */
  DataView * apply(const DataType& dtype);

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
   * \brief Return total number of bytes allocated by this DataView object.
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
  size_t getNumberOfElements() const
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
   *  \brief Private ctor that creates associates meta-buffer with given name
   *         in given parent group.
   */
  DataView( const std::string& name,
	    DataGroup * const owning_group,
	    void * buffer_context, void * metabuffer);

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

  /// Meta-buffer support
  void *m_buffer_context;

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

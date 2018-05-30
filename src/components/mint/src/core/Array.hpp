/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#ifndef MINT_ARRAY_HPP_
#define MINT_ARRAY_HPP_

#include "axom/Macros.hpp"           // for disable copy/assignment macro
#include "axom_utils/Utilities.hpp"  // for memory allocation functions
#include "mint/config.hpp"           // for mint::IndexType definition

#include "slic/slic.hpp"            // for slic logging macros

#ifdef MINT_USE_SIDRE
#include "sidre/sidre.hpp"          // for sidre::View definition
#endif

// C/C++ includes
#include <cstring>                  // for std::memcpy

namespace axom
{
namespace mint
{

/* Provided so that 0 doesn't convert to nullptr and lead to ambiguous
 * constructor calls. */
namespace internal
{
constexpr IndexType ZERO = 0;
}

/*!
 * \class Array
 *
 * \brief Provides a generic multi-component array container.
 *
 *  The Array class provides a generic multi-component array container with
 *  dynamic re-allocation and insertion. Each element in the array is a tuple
 *  consisting of 1 or more components, which are stored contiguously.
 *
 *  Depending on which constructor is used, the Array object can have three
 *  different underlying storage types:
 *
 *  * <b> Native Storage </b> <br />
 *
 *     When using native storage, the Array object manages all memory.
 *     Typically, the Array object will allocate extra space to facilitate
 *     the insertion of new elements and minimize the number of reallocations.
 *     At any given instance, the actual capacity of the array (i.e., total
 *     number of tuples that the Array can hold) can be queried by calling the
 *     capacity() function. When all extra memory is exhausted, inserting a
 *     new element triggers a re-allocation. Each time a re-allocation occurs
 *     extra space is allocated according to the <em> resize_ratio </em>
 *     parameter, which is set to 2.0 by default. To return all extra memory,
 *     an application can call shrink().
 *
 *     \note Once the Array goes out-of-scope, all memory associated with it
 *      is de-allocated and returned to the system.
 *
 *  * <b> External Storage </b> <br />
 *
 *    An Array object may be constructed from an external, user-supplied buffer
 *    consisting of the given number of tuples and specified number of
 *    components per tuple. In this case, the Array object does not own the
 *    memory. Instead, the Array object makes a shallow copy of the pointer.
 *
 *    \warning An Array object that points to an external buffer has a fixed
 *     size and cannot be dynamically resized.
 *
 *    \note When the Array object is deleted, it does not de-allocate the
 *     user-supplied buffer, since it does not manage any memory.
 *
 *  * <b> Sidre </b> <br />
 *
 *    An Array object may also be constructed from a sidre::View. From the
 *    application's perspective, other than the fact that a different
 *    constructor is used, all array operations can be performed transparently.
 *    The size of the Array can grow as needed, but instead all memory
 *    management is delegated to Sidre.
 *
 *    \note When the Array object is deleted, it does not delete the associated
 *     data in Sidre, since, Sidre owns the data.
 *
 * \warning Reallocations tend to be costly operations in terms of performance.
 *  Use `reserve()` when the number of nodes is known a priori, or opt to
 *  use a constructor that takes an actual size and capacity when possible.
 *
 * \tparam T the type of the values to hold.
 *
 * \see sidre::Group
 * \see sidre::View
 */
template< typename T >
class Array
{
public:
  static constexpr double DEFAULT_RESIZE_RATIO = 2.0;
  static constexpr IndexType MIN_DEFAULT_CAPACITY = 32;

public:

/// \name Native Storage Array Constructors
/// @{

  /*!
   * \brief Default constructor. Disabled.
   */
  Array() = delete;

  /*!
   * \brief Constructs an Array instance with the given number of tuples.
   *
   * \param [in] num_tuples the number of tuples the Array holds.
   * \param [in] num_components the number of components per tuple.
   * \param [in] capacity the number of tuples to allocate space for.
   *
   * \note The last argument is optional. If not specified, the
   *  capacity of the array will be initialized to
   *  num_tuples * DEFAULT_RESIZE_RATIO.
   *
   * \pre num_tuples >= 0
   * \pre num_components >= 1
   *
   * \post capacity() >= size()
   * \post size() == num_tuples
   * \post numComponents() == num_components
   * \post getResizeRatio() == DEFAULT_RESIZE_RATIO
   */
  Array( IndexType num_tuples, IndexType num_components,
         IndexType capacity=USE_DEFAULT );

/// @}

/// \name External Storage Array Constructors
/// @{

  /*!
   * \brief Constructs an Array instance with the given number of tuples from
   *  an external data buffer.
   *
   * \param [in] data the external data this Array will wrap.
   * \param [in] num_tuples the number of tuples in the Array.
   * \param [in] num_components the number of values per tuple.
   * \param [in] capacity the capacity of the external buffer.
   *
   * \pre data != AXOM_NULLPTR
   * \pre num_tuples > 0
   * \pre num_components >= 1
   *
   * \post numComponents() == num_components
   * \post getResizeRatio == 0.0
   *
   * \note If no capacity is specified then it will default to the number of
   *  tuples.
   *
   * \note This constructor wraps the supplied buffer and does not own the data.
   *  Consequently, the Array instance cannot be reallocated.
   */
  Array( T* data, IndexType num_tuples, IndexType num_components,
         IndexType capacity=USE_DEFAULT );

/// @}

/// \name Sidre Array constructors
/// @{

#ifdef MINT_USE_SIDRE

  /*!
   * \brief Creates an Array instance from a sidre::View that already has data.
   *
   * \param [in] view the sidre::View that holds this Array's data.
   *
   * \note The Sidre view shape has two dimensions. The first dimension
   *  corresponds to the max capacity of the array and the second corresponds to
   *  the number of components per tuple.
   *
   * \pre view != AXOM_NULLPTR
   * \pre view->isEmpty() == false.
   * \pre view->getNumDimensions() == 2
   *
   * \post capacity() == view->getDimension(0)
   * \post numComponents() == view->getDimension(1)
   * \post getResizeRatio() == DEFAULT_RESIZE_RATIO
   */
  Array( sidre::View* view );

  /*!
   * \brief Creates an Array instance of `num_tuples` size, where each
   *  tuple consists of `num_components` values and populates the associated
   *  sidre::View.
   *
   * \param [in] view the sidre::View that will hold this Array's data.
   * \param [in] num_tuples the number of tuples accounted for in the Array.
   * \param [in] num_components the number of values per tuple.
   * \param [in] capacity the number of tuples to allocate space for.
   *
   * \note The last argument is optional. If not specified, the
   *  capacity of the array will be initialized to
   *  num_tuples * DEFAULT_RESIZE_RATIO.
   *
   * \note The view is expected to be empty and will be populated to hold this
   *  Array's data.
   *
   * \note The Sidre view shape has two dimensions. The first dimension
   *  corresponds to the number of tuples and the second corresponds to
   *  the number of components per tuple.
   *
   * \pre view != AXOM_NULLPTR
   * \pre view->isEmpty() == true
   * \pre num_tuples >= 1
   * \pre num_components >= 1
   *
   * \post view->getNumDimensions() == 2
   * \post view->isEmpty() == false
   * \post size() == num_tuples.
   * \post numComponents() == num_components
   * \post getResizeRatio() == DEFAULT_RESIZE_RATIO
   */
  Array( sidre::View* view, IndexType num_tuples, IndexType num_components,
         IndexType capacity=USE_DEFAULT );

#endif

/// @}

  /*!
   * Destructor. Free's the associated buffer unless the memory is external
   * or owned by Sidre.
   */
  ~Array();

/// \name Array tuple access operators
/// @{

  /*!
   * \brief Accessor, returns a reference to the given component of the
   *  specified tuple.
   * \param [in] pos the tuple to querry.
   * \param [in] component the component to return.
   * \return a reference to the given component of the specified tuple.
   *
   * \pre 0 <= pos < size()
   * \pre 0 <= component < numComponents()
   */
  inline T & operator()( IndexType pos, IndexType component=0 )
  {
    SLIC_ASSERT( pos >= 0 );
    SLIC_ASSERT( pos < m_num_tuples );
    SLIC_ASSERT( component >= 0 );
    SLIC_ASSERT( component < m_num_components );
    return m_data[ pos * m_num_components + component ];
  }

  /*!
   * \brief Accessor, returns a reference to the given value.
   * \param [in] idx the position of the value to return.
   * \return a reference to the given value.
   * \note equivalent to *(array.getData() + idx).
   *
   * \pre 0 <= idx < m_num_tuples * m_num_components
   */
  inline T & operator[]( IndexType idx )
  {
    SLIC_ASSERT( idx >= 0 );
    SLIC_ASSERT( idx < m_num_tuples * m_num_components );
    return m_data[ idx ];
  }

  /*!
   * \brief Constant accessor, returns a reference to the given component of the
   *  specified tuple.
   * \param [in] pos the tuple to querry.
   * \param [in] component the component to return.
   * \return a reference to the given component of the specified tuple.
   */
  inline const T & operator()( IndexType pos, IndexType component=0 ) const
  { return m_data[ pos * m_num_components + component ]; }

  /*!
   * \brief Return a pointer to the array of data.
   * \return a pointer to the array of data.
   */
  inline T* getData() { return m_data; }

  /*!
   * \brief Return a constant pointer to the array of data.
   * \return a constant pointer to the array of data.
   */
  inline const T* getData() const { return m_data; }

/// @}

/// \name Array methods to modify the data.
/// @{

  /*!
   * \brief Set all the values of the array.
   * \param [in] value the value to set to.
   */
  inline void fill( const T& value )
  { std::fill_n( m_data, m_num_tuples * m_num_components, value ); }

  /*!
   * \brief Append a value to the end of the array.
   * \param [in] value the value to append.
   * \pre m_num_components == 1.
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  inline void append( const T& value );

  /*!
   * \brief Append tuples to the end of the array.
   * \param [in] tuples the tuples to append.
   * \param [in] n the number of tuples to append.
   * \note It's assumed that tuples is of length n * m_num_components.
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  inline void append( const T* tuples, IndexType n );

  /*!
   * \brief Write tuples into existing spots in the array.
   * \param [in] tuples the new tuples to be written.
   * \param [in] n the number of tuples to write.
   * \param [in] pos the position at which to begin writing.
   * \pre pos + n <= m_num_tuples.
   */
  inline void set( const T* tuples, IndexType n, IndexType pos );

  /*!
   * \brief Insert a value into the array at the given position.
   * \param [in] value the value to insert.
   * \param [in] pos the position at which to insert.
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  inline void insert( const T& value, IndexType pos );

  /*!
   * \brief Insert tuples into the array at the given position.
   * \param [in] tuples the tuples to insert.
   * \param [in] n the number of tuples to insert.
   * \param [in] pos the position at which to begin the insertion.
   * \pre pos < m_num_tuples.
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  inline void insert( const T* tuples, IndexType n, IndexType pos );

  /*!
   * \brief Insert multiple values at the given position.
   * \param [in] n the number of values to insert.
   * \param [in] pos the position to insert at.
   * \param [in] value the value to insert.
   */
  inline void emplace( IndexType n, IndexType pos, const T& value=T() );

/// @}

/// \name Array methods to query and set attributes
/// @{

  /*!
   * \brief Return the number of tuples allocated for the data array.
   * \return the amount of space allocated for the data array.
   */
  inline IndexType capacity() const { return m_capacity; }

  /*!
   * \brief Increase the capacity. Does nothing if the new capacity is less
   *  than the current capacity.
   * \param [in] capacity the new number of tuples to allocate.
   */
  inline void reserve( IndexType capacity )
  { if ( capacity > m_capacity ) setCapacity( capacity ); }

  /*!
   * \brief Shrink the capacity to be equal to the size.
   */
  inline void shrink() { setCapacity( m_num_tuples ); }

  /*!
   * \brief Returns true iff the Array stores no elements.
   * \return true iff the Array stores no elements.
   * \note If the Array is empty the capacity can still be greater than zero.
   */
  inline bool empty() const { return m_num_tuples == 0; }

  /*!
   * \brief Return the number of tuples stored in the data array.
   * \return the number of tuples stored in the data array.
   */
  inline IndexType size() const { return m_num_tuples; }

  /*!
   * \brief Update the number of tuples stored in the data array.
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  inline void resize( IndexType new_num_tuples );

  /*!
   * \brief Get the ratio by which the capacity increases upon dynamic resize.
   * \return The ration by which the capacity increases upon dynamic resize.
   */
  inline double getResizeRatio() const { return m_resize_ratio; }

  /*!
   * \brief Set the ratio by which the capacity increases upon dynamic resize.
   * \param [in] ratio the new resize ratio.
   */
  inline void setResizeRatio( double ratio ) { m_resize_ratio = ratio; }

  /*!
   * \brief Get the chunk size of all allocations.
   * \return the chunk size of all allocations.
   */
  inline IndexType numComponents() const { return m_num_components; }

  /*!
   * \brief Checks if this array instance points to an external buffer
   * \return status true iff the external buffer constructor was called.
   */
  inline bool isExternal() const { return m_is_external; }

  /*!
   * \brief Checks if this array instance is in sidre.
   * \return status true iff a sidre constructor was called.
   */
  inline bool isInSidre() const
  {
    #ifdef MINT_USE_SIDRE
    return m_view != AXOM_NULLPTR;
    #else
    return false;
    #endif
  }

#ifdef MINT_USE_SIDRE
  /*!
   * \brief Return a pointer to the sidre::View that this Array wraps.
   * \return a const pointer to a sidre::View.
   */
  inline const sidre::View* getView() const { return m_view; }
#endif

/// @}

private:

  /*!
   * \brief Make space for a subsequent insertion into the array.
   * \param [in] n the number of tuples to insert.
   * \param [in] pos the position at which to begin the insertion.
   * \return a pointer to the beginning of the insertion space.
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  inline T* reserveForInsert( IndexType n, IndexType pos );

  /*!
   * \brief Update the number of tuples.
   * \param [in] new_num_tuples the new number of tuples.
   */
  inline void updateNumTuples( IndexType new_num_tuples );

  /*!
   * \brief Set the number of tuples allocated for the data array.
   * \param [in] capacity the new number of tuples to allocate.
   */
  inline void setCapacity( IndexType new_capacity );

  /*!
   * \brief Reallocates the data array when the size exceeds the capacity.
   * \param [in] new_num_tuples the number of tuples which exceeds the current
   *  capacity.
   */
  inline void dynamicRealloc( IndexType new_num_tuples );



#ifdef MINT_USE_SIDRE

  /*!
   * \brief Return the sidre::TypeID corresponding to T. This function
   *  handles when T is not an enum.
   */
  template < typename U = T >
  static constexpr
  typename std::enable_if< !std::is_enum< U >::value, sidre::TypeID >::type
  sidreTypeId()
  { return sidre::detail::SidreTT< U >::id; }

  /*!
   * \brief Return the sidre::TypeID corresponding to T. This function handles
   *  when T is an enum.
   */
  template < typename U = T >
  static constexpr
  typename std::enable_if< std::is_enum< U >::value, sidre::TypeID >::type
  sidreTypeId()
  {
    return sidre::detail::SidreTT<
      typename std::underlying_type< U >::type >::id;
  }

  /*!
   * \brief Describes m_view as having dimensions
   * (m_num_tuples, m_num_components).
   */
  inline void describeView();

  /*!
   * \brief Given a non-empty sidre::View of dimension 2, returns the length
   *  of the given dimension.
   * \param [in] view the sidre::View to examine.
   * \param [in] dim the dimension (0 or 1) to return the length of.
   * \return The length of dimension dim of view.
   *
   * \pre 0 <= dim <= 1
   */
  inline IndexType getViewShape( int dim ) const;

  /*!
   * \brief Allocates space within the Array's sidre::View.
   */
  inline void reallocViewData();

  sidre::View* m_view;
#endif

  T* m_data;
  IndexType m_num_tuples;
  IndexType m_capacity;
  IndexType m_num_components;
  double m_resize_ratio;
  bool m_is_external;

  DISABLE_COPY_AND_ASSIGNMENT( Array );
  DISABLE_MOVE_AND_ASSIGNMENT( Array );
};


//------------------------------------------------------------------------------
//                            Array IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template< typename T >
Array< T >::Array( IndexType num_tuples, IndexType num_components,
                   IndexType capacity ) :
#ifdef MINT_USE_SIDRE
  m_view( AXOM_NULLPTR ),
#endif
  m_data( AXOM_NULLPTR ),
  m_num_tuples( num_tuples ),
  m_capacity( 0 ),
  m_num_components( num_components ),
  m_resize_ratio( DEFAULT_RESIZE_RATIO ),
  m_is_external( false )
{
  SLIC_ERROR_IF( m_num_tuples < 0,
                 "Number of tuples (" << m_num_tuples << ") " <<
                 "cannot be negative." );
  SLIC_ERROR_IF( m_num_components <= 0,
                 "Components per tuple (" << m_num_components << ") " <<
                 "must be greater than zero." );
  SLIC_ERROR_IF( capacity >= 0 && m_num_tuples > capacity,
                 "Number of tuples (" << m_num_tuples << ") " <<
                 "Number of components (" << m_num_components << ") " <<
                 "cannot be greater than the specified capacity " <<
                 "(" << capacity << ")." );

  if ( capacity == USE_DEFAULT )
  {
    capacity = ( m_num_tuples > MIN_DEFAULT_CAPACITY ) ?
               m_num_tuples : MIN_DEFAULT_CAPACITY;
  }
  setCapacity( capacity );

  // sanity checks
  SLIC_ASSERT( capacity >= 0 );
  if ( capacity > 0 )
  {
    SLIC_ASSERT( m_data != AXOM_NULLPTR );
  }
  SLIC_ASSERT( m_num_tuples >= 0 );
  SLIC_ASSERT( m_num_components >= 1 );
}

//------------------------------------------------------------------------------
template< typename T >
Array< T >::Array( T* data, IndexType num_tuples, IndexType num_components,
                   IndexType capacity ) :
#ifdef MINT_USE_SIDRE
  m_view( AXOM_NULLPTR ),
#endif
  m_data( data ),
  m_num_tuples( num_tuples ),
  m_capacity( 0 ),
  m_num_components( num_components ),
  m_resize_ratio( 0.0 ),
  m_is_external( true )
{
  if ( capacity == USE_DEFAULT )
  {
    m_capacity = num_tuples;
  }
  else
  {
    m_capacity = capacity;
  }

  SLIC_ERROR_IF( m_num_tuples < 0,
                 "Number of tuples (" << m_num_tuples << ") " <<
                 "cannot be negative." );
  SLIC_ERROR_IF( m_num_components < 1,
                 "Components per tuple (" << m_num_components << ") " <<
                 "must be greater than zero." );
  SLIC_ERROR_IF( m_num_tuples > m_capacity,
                 "Number of tuples (" << m_num_tuples << ") " <<
                 "cannot be greater than the tuple capacity " <<
                 "(" << m_capacity << ")." );
  SLIC_ERROR_IF( m_data == AXOM_NULLPTR, "specified array must not be NULL." );
}

#ifdef MINT_USE_SIDRE

//------------------------------------------------------------------------------
template< typename T >
Array< T >::Array( sidre::View* view ) :
  m_view( view ),
  m_data( AXOM_NULLPTR ),
  m_num_tuples(0),
  m_capacity(0),
  m_num_components(0),
  m_resize_ratio( DEFAULT_RESIZE_RATIO ),
  m_is_external( false )
{
  SLIC_ERROR_IF( m_view == AXOM_NULLPTR, "Provided View cannot be null." );
  SLIC_ERROR_IF( m_view->isEmpty(), "Provided View cannot be empty." );

  m_num_tuples       = getViewShape(0);
  m_num_components = getViewShape(1);

  IndexType buffer_size = m_view->getBuffer()->getNumElements();
  SLIC_ERROR_IF( buffer_size % m_num_components != 0,
                 "The buffer size (" << buffer_size << ") " <<
                 "is not a multiple of the number of componets " <<
                 "(" << m_num_components << ")." );
  m_capacity = buffer_size / m_num_components;

  SLIC_ERROR_IF( m_num_tuples < 0,
                 "Number of tuples (" << m_num_tuples << ") " <<
                 "cannot be negative." );

  SLIC_ERROR_IF( m_num_components <= 0,
                 "Number of components (" << m_num_components << ") " <<
                 "must be greater than 0." );

  SLIC_ERROR_IF( m_num_tuples > m_capacity,
                 "Number of tuples (" << m_num_tuples << ") " <<
                 "cannot be greater than the tuple capacity " <<
                 "(" << m_capacity << ")." );

  sidre::TypeID view_type = m_view->getTypeID();
  sidre::TypeID T_type = sidreTypeId();
  SLIC_ERROR_IF( view_type != T_type,
                 "View data type (" << view_type << ")" <<
                 "differs from this Array type (" << T_type << ")." );

  m_data = static_cast< T* >( m_view->getVoidPtr() );
  SLIC_ERROR_IF( m_data == AXOM_NULLPTR && m_capacity > 0,
                 "View returned a null pointer when the capacity " <<
                 "is greater than zero." );
}

//------------------------------------------------------------------------------
template< typename T >
Array< T >::Array( sidre::View* view, IndexType num_tuples,
                   IndexType num_components, IndexType capacity ) :
  m_view( view ),
  m_data( AXOM_NULLPTR ),
  m_num_tuples( num_tuples ),
  m_capacity( 0 ),
  m_num_components( num_components ),
  m_resize_ratio( DEFAULT_RESIZE_RATIO ),
  m_is_external( false )
{
  SLIC_ERROR_IF( m_view == AXOM_NULLPTR, "Provided View cannot be null." );
  SLIC_ERROR_IF( !m_view->isEmpty(), "View must be empty." );
  SLIC_ERROR_IF( m_num_tuples < 0,
                 "Number of tuples (" << m_num_tuples << ") " <<
                 "cannot be negative." );
  SLIC_ERROR_IF( m_num_components <= 0,
                 "Components per tuple (" << m_num_components << ") " <<
                 "must be greater than 0." );

  if ( capacity == USE_DEFAULT )
  {
    capacity = ( m_num_tuples > MIN_DEFAULT_CAPACITY ) ?
               m_num_tuples : MIN_DEFAULT_CAPACITY;
  }
  SLIC_ERROR_IF( m_num_tuples > capacity,
                 "Number of tuples (" << m_num_tuples << ") " <<
                 "cannot be greater than the tuple capacity " <<
                 "(" << capacity << ")." );

  setCapacity( capacity );

  // sanity checks
  SLIC_ASSERT( capacity >= 0 );
  if ( capacity > 0 )
  {
    SLIC_ASSERT( m_data != AXOM_NULLPTR );
  }
  SLIC_ASSERT( m_num_tuples >= 0 );
  SLIC_ASSERT( m_num_components >= 1 );
}
#endif

//------------------------------------------------------------------------------
template< typename T >
Array< T >::~Array()
{
#ifdef MINT_USE_SIDRE
  if ( m_view == AXOM_NULLPTR && m_data != AXOM_NULLPTR && !m_is_external )
  {
    utilities::free( m_data );
  }
  m_view = AXOM_NULLPTR;
#else
  if ( m_data != AXOM_NULLPTR && !m_is_external )
  {
    utilities::free( m_data );
  }
#endif
  m_data = AXOM_NULLPTR;
}


//------------------------------------------------------------------------------
template< typename T >
inline void Array< T >::append( const T& value )
{
  SLIC_ASSERT_MSG( m_num_components == 1, "Number of components must be 1." );

  IndexType new_size = m_num_tuples + 1;
  if ( new_size > m_capacity )
  {
    dynamicRealloc( new_size );
  }

  m_data[ m_num_tuples ] = value;
  updateNumTuples( new_size );
}

//------------------------------------------------------------------------------
template< typename T >
inline void Array< T >::append( const T* tuples, IndexType n )
{
  IndexType new_size = m_num_tuples + n;
  if ( new_size > m_capacity )
  {
    dynamicRealloc( new_size );
  }

  T* cur_end = m_data + m_num_tuples * m_num_components;
  std::memcpy( cur_end, tuples, n * m_num_components * sizeof(T) );
  updateNumTuples( new_size );
}

//------------------------------------------------------------------------------
template< typename T >
inline void Array< T >::set( const T* tuples, IndexType n, IndexType pos )
{
  SLIC_ASSERT( tuples != AXOM_NULLPTR );
  SLIC_ASSERT( pos >= 0 );
  SLIC_ASSERT( pos + n <= m_num_tuples );

  T* set_position = &m_data[ pos * m_num_components ];
  IndexType byte_size = n * m_num_components * sizeof(T);
  std::memcpy( set_position, tuples, byte_size );
}

//------------------------------------------------------------------------------
template< typename T >
inline void Array< T >::insert( const T& value, IndexType pos )
{
  SLIC_ASSERT_MSG( m_num_components == 1, "Number of components must be 1." );
  reserveForInsert( 1, pos );
  m_data[ pos ] = value;
}

//------------------------------------------------------------------------------
template< typename T >
inline void Array< T >::insert( const T* tuples, IndexType n, IndexType pos )
{
  SLIC_ASSERT( tuples != AXOM_NULLPTR );
  T* insert_pos = reserveForInsert( n, pos );
  std::memcpy( insert_pos, tuples, n * m_num_components * sizeof(T) );
}

//------------------------------------------------------------------------------
template< typename T >
inline void Array< T >::emplace( IndexType n, IndexType pos, const T& value )
{
  T* insert_pos = reserveForInsert( n, pos );
  std::fill_n( insert_pos, n * numComponents(), value );
}


//------------------------------------------------------------------------------
template< typename T >
inline void Array< T >::resize( IndexType new_num_tuples )
{
  SLIC_ASSERT( new_num_tuples >= 0 );

  if ( new_num_tuples > m_capacity )
  {
    dynamicRealloc( new_num_tuples );
  }

  updateNumTuples( new_num_tuples );
}

//------------------------------------------------------------------------------
template< typename T >
inline T* Array< T >::reserveForInsert( IndexType n, IndexType pos )
{
  SLIC_ASSERT( n >= 0 );
  SLIC_ASSERT( pos >= 0 );
  SLIC_ASSERT( pos <= m_num_tuples );

  if ( n == 0 )
  {
    return m_data + pos * m_num_components;
  }

  IndexType new_size = m_num_tuples + n;
  if ( new_size > m_capacity )
  {
    dynamicRealloc( new_size );
  }

  T* const insert_pos = m_data + pos * m_num_components;
  T* cur_pos = m_data + (m_num_tuples * m_num_components) - 1;
  for ( ; cur_pos >= insert_pos ; --cur_pos )
  {
    *(cur_pos + n * m_num_components) = *cur_pos;
  }

  updateNumTuples( new_size );
  return insert_pos;
}

//------------------------------------------------------------------------------
template< typename T >
inline void Array< T >::updateNumTuples( IndexType new_num_tuples )
{
  SLIC_ASSERT( new_num_tuples >= 0 );
  SLIC_ASSERT( new_num_tuples <= m_capacity );
  m_num_tuples = new_num_tuples;

#ifdef MINT_USE_SIDRE
  if ( m_view != AXOM_NULLPTR )
  {
    describeView();
  }
#endif
}

//------------------------------------------------------------------------------
template< typename T >
inline void Array< T >::setCapacity( IndexType new_capacity )
{
  SLIC_ASSERT( new_capacity >= 0 );

  if ( m_is_external && new_capacity <= m_capacity )
  {
    return;
  }

  SLIC_ERROR_IF( m_is_external, "Cannot change the capacity of external data.");

  m_capacity = new_capacity;
  if ( m_capacity < m_num_tuples )
  {
    updateNumTuples( m_capacity );
  }

#ifdef MINT_USE_SIDRE
  if ( m_view != AXOM_NULLPTR )
  {
    return reallocViewData();
  }
#endif

  m_data = utilities::realloc( m_data, m_capacity * m_num_components );
  SLIC_ERROR_IF( m_data == AXOM_NULLPTR && m_capacity > 0,
                 "Array reallocation failed." );
}

//------------------------------------------------------------------------------
template< typename T >
inline void Array< T >::dynamicRealloc( IndexType new_num_tuples )
{
  SLIC_ERROR_IF( m_is_external, "Cannot change the capacity of external data.");
  SLIC_ERROR_IF( m_resize_ratio < 1.0, "Resize ratio of " << m_resize_ratio <<
                 " doesn't support dynamic resizing");
  m_capacity = new_num_tuples * m_resize_ratio + 0.5;

#ifdef MINT_USE_SIDRE
  if ( m_view != AXOM_NULLPTR )
  {
    return reallocViewData();
  }
#endif

  m_data = utilities::realloc( m_data, m_capacity * m_num_components );
  SLIC_ERROR_IF( m_data == AXOM_NULLPTR && m_capacity > 0,
                 "Array reallocation failed." );
}

#ifdef MINT_USE_SIDRE

//------------------------------------------------------------------------------
template< typename T >
inline void Array< T >::describeView()
{
  SLIC_ASSERT( m_view != AXOM_NULLPTR );

  static constexpr sidre::TypeID T_type = sidreTypeId();
  sidre::SidreLength dims[2];
  dims[0] = m_num_tuples;
  dims[1] = m_num_components;

  m_view->apply( T_type, 2, dims );
}

//------------------------------------------------------------------------------
template< typename T >
inline IndexType Array< T >::getViewShape( int dim ) const
{
  SLIC_ERROR_IF( dim > 1, "Only two dimensional views supported." );
  SLIC_ERROR_IF( m_view->isEmpty(), "view cannot be empty." );
  SLIC_ERROR_IF( m_view->getNumDimensions() != 2,
                 "view must have dimension 2.");

  sidre::SidreLength dims[ 2 ];
  m_view->getShape( 2, dims );
  return dims[ dim ];
}

//------------------------------------------------------------------------------
template< typename T >
inline void Array< T >::reallocViewData()
{
  if ( m_view->isEmpty() )
  {
    constexpr sidre::TypeID T_type = sidreTypeId();
    m_view->allocate( T_type, m_capacity * m_num_components );
  }
  else
  {
    m_view->reallocate( m_capacity * m_num_components );
  }

  describeView();
  m_data = static_cast< T* >( m_view->getVoidPtr() );

  SLIC_ERROR_IF( m_data == AXOM_NULLPTR && m_capacity > 0,
                 "Array reallocation failed." );
}

#endif  /* MINT_USE_SIDRE */

} /* namespace mint */
} /* namespace axom */

#endif /* MINT_ARRAY_HPP_ */

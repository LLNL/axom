/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
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

#ifndef MINT_UTILS_ARRAY_HXX_
#define MINT_UTILS_ARRAY_HXX_

#ifndef MINT_USE_SIDRE
#define MINT_USE_SIDRE
#endif


#include "slic/slic.hpp"                // for slic macros
#include "mint/DataTypes.hpp"           // for localIndex
#include "axom_utils/Utilities.hpp"     // for allocation

#ifdef MINT_USE_SIDRE
#include "sidre/sidre.hpp"              // for View
#endif

#include <cstring>                      // for std::memcpy

namespace axom
{
namespace mint
{

/*!
 * \class Array
 * \brief Provides an array with dynamic reallocation and insertion.
 * \tparam T the type of the values to hold.
 */
template< typename T >
class Array
{
private:
  static constexpr double DEFAULT_RESIZE_RATIO = 2.0;

public:

//@{
//!  @name Array constructors for use without sidre.

  /*!
   * \brief Default constructor.
   */
  Array() = delete;

  /*!
   * \brief Constructs a Array instance with the give number of tuples.
   * \param [in] num_tuples the number of tuples the Array holds.
   * \param [in] num_components the number of components per tuple.
   */
  Array( localIndex num_tuples, localIndex num_components );

  /*!
   * \brief Constructs a Array instance with the given number of tuples
   *  and capacity.
   * \param [in] capacity the number of tuples to allocate space for.
   * \param [in] num_tuples the number of tuples accounted for in the Array.
   * \param [in] num_components the number of values per tuple.
   */
  Array( localIndex capacity, localIndex num_tuples,
         localIndex num_components );

  /*!
   * \brief Constructs a Array instance with external data and the given number
   *  of tuples.
   * \param [in] data the external data this Array will wrap.
   * \param [in] num_tuples the number of tuples accounted for in the Array.
   * \param [in] num_components the number of values per tuple.
   */
  Array( T* data, localIndex num_tuples, localIndex num_components );

//@}



//@{
//!  @name Array constructors for use with sidre.
#ifdef MINT_USE_SIDRE

  /*!
   * \brief Constructor for use with a sidre::View that already has data.
   * \param [in] view the sidre::View that holds this Array's data.
   * \param [in] num_tuples the number of tuples accounted for in the Array.
   */
  Array( sidre::View* view, localIndex num_tuples );

  /*!
   * \brief Constructor for use with an empty sidre::View with the given number
   *  of tuples and components.
   * \param [in] view the sidre::View that will hold this Array's data.
   * \param [in] num_tuples the number of tuples accounted for in the Array.
   * \param [in] num_components the number of values per tuple.
   */
  Array( sidre::View* view, localIndex num_tuples,
         localIndex num_components );

  /*!
   * \brief Constructor for use with an empty sidre::View with the given number
   *  of tuples and components and with the given capacity.
   * \param [in] view the sidre::View that will hold this Array's data.
   * \param [in] capacity the number of tuples to allocate space for.
   * \param [in] num_tuples the number of tuples accounted for in the Array.
   * \param [in] num_components the number of values per tuple.
   */
  Array( sidre::View* view, localIndex capacity, localIndex num_tuples,
         localIndex num_components );

#endif
//@}


  /*!
   * Destructor. Free's the associated buffer unless it's owned by sidre.
   */
  ~Array();


//@{
//!  @name Array methods to access the data.

  /*!
   * \brief Accessor, returns a reference to the given component of the
   *  specified tuple.
   * \param [in] pos the tuple to querry.
   * \param [in] component the component to return.
   * \return a reference to the given component of the specified tuple.
   */
  inline T & operator()( localIndex pos, localIndex component=0 )
  { return m_data[ pos * m_num_components + component ]; }

  /*!
   * \brief Constant accessor, returns a reference to the given component of the
   *  specified tuple.
   * \param [in] pos the tuple to querry.
   * \param [in] component the component to return.
   * \return a reference to the given component of the specified tuple.
   */
  inline T & operator()( localIndex pos, localIndex component=0 ) const
  { return m_data[ pos * m_num_components + component ]; }

  /*!
   * \brief Return a pointer to the array of data.
   * \return a pointer to the array of data.
   */
  inline T* getData()
  { return m_data; }

  /*!
   * \brief Return a constant pointer to the array of data.
   * \return a constant pointer to the array of data.
   */
  inline const T* getData() const
  { return m_data; }

//@}



//@{
//!  @name Array methods to modify the data.

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
  inline void append( const T* tuples, localIndex n );

  /*!
   * \brief Write tuples into existing spots in the array.
   * \param [in] tuples the new tuples to be written.
   * \param [in] n the number of tuples to write.
   * \param [in] pos the position at which to begin writing.
   * \pre pos + n <= m_num_tuples.
   */
  inline void set( const T* tuples, localIndex n, localIndex pos );

  /*!
   * \brief Insert tuples into the array at the given position.
   * \param [in] tuples the tuples to insert.
   * \param [in] n the number of tuples to insert.
   * \param [in] pos the position at which to begin the insertion.
   * \pre pos < m_num_tuples.
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  inline void insert( const T* tuples, localIndex n, localIndex pos );

  /*!
   * \brief Insert a value into the array at the given position.
   * \param [in] value the value to insert.
   * \param [in] pos the position at which to insert.
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  inline void insert( const T& value, localIndex pos );

//@}


//@{
//!  @name Array methods to query and set attributes


  /*!
   * \brief Return the number of tuples allocated for the data array.
   * \return the amount of space allocated for the data array.
   */
  inline localIndex getCapacity() const
  { return m_capacity; }

  /*!
   * \brief Increase the capacity. Does nothing if the new capacity is less
   *  than the current capacity.
   * \param [in] capacity the new number of tuples to allocate.
   */
  inline void reserve( localIndex capacity )
  { if ( capacity > m_capacity ) setCapacity( capacity ); }

  /*!
   * \brief Shrink the capacity to be equal to the size.
   */
  inline void shrink()
  { setCapacity(m_num_tuples); }

  /*!
   * \brief Returns true iff the Array stores no elements.
   * \return true iff the Array stores no elements.
   * \note If the Array is empty the capacity can still be greater than zero.
   */
  inline bool empty() const
  { return m_num_tuples == 0; }

  /*!
   * \brief Return the number of tuples stored in the data array.
   * \return the number of tuples stored in the data array.
   */
  inline localIndex size() const { return m_num_tuples; }

  /*!
   * \brief Update the number of tuples stored in the data array.
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  inline void resize( localIndex num_tuples );

  /*!
   * \brief Get the ratio by which the capacity increases upon dynamic resize.
   * \return The ration by which the capacity increases upon dynamic resize.
   */
  inline double getResizeRatio() const
  { return m_resize_ratio; }

  /*!
   * \breif Set the ratio by which the capacity increases upon dynamic resize.
   * \param [in] ratio the new resize ratio.
   */
  inline void setResizeRatio( double ratio )
  { m_resize_ratio = ratio; }

  /*!
   * \brief Get the chunk size of all allocations.
   * \return the chunk size of all allocations.
   */
  inline localIndex getNumComponents() const
  { return m_num_components; }


#ifdef MINT_USE_SIDRE
  /*!
   * \brief Return a pointer to the sidre::View that this Array wraps.
   * \return a const pointer to a sidre::View.
   */
  inline sidre::View* getView()
  { return m_view; }
#endif

//@}


private:

  /*!
   * \brief Set the number of tuples allocated for the data array.
   * \param [in] capacity the new number of tuples to allocate.
   */
  inline void setCapacity( localIndex capacity );

  /*!
   * \brief Reallocates the data array when the size exceeds the capacity.
   * \param [in] new_num_tuples the number of tuples which exceeds the current
   *  capacity.
   */
  inline void dynamicRealloc( localIndex new_num_tuples );

  /*!
   * \brief Make space for a subsequent insertion into the array.
   * \param [in] n the number of tuples to insert.
   * \param [in] pos the position at which to begin the insertion.
   * \return a pointer to the beginning of the insertion space.
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  inline T* reserveForInsert( localIndex n, localIndex pos );

#ifdef MINT_USE_SIDRE
  /*!
   * \brief Given a non-empty sidre::View of dimension 2, returns the length
   *  of the given dimension.
   * \param [in] view the sidre::View to examine.
   * \param [in] dim the dimension (0 or 1) to return the length of.
   * \return The length of dimension dim of view.
   */
  inline localIndex getViewShape( int dim ) const;

  /*!
   * \brief Allocates space within the Array's sidre::View.
   */
  inline void allocate_view_data();


  sidre::View* m_view;
#endif

  T* m_data;
  localIndex m_num_tuples;
  localIndex m_capacity;
  localIndex m_num_components;
  double m_resize_ratio;
  bool m_is_external;
};


//------------------------------------------------------------------------------
//                            Array IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template< typename T >
Array< T >::Array( localIndex num_tuples, localIndex num_components ) :
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

  setCapacity( m_num_tuples * m_resize_ratio );
}

//------------------------------------------------------------------------------
template< typename T >
Array< T >::Array( localIndex capacity, localIndex num_tuples,
                   localIndex num_components ) :
#ifdef MINT_USE_SIDRE
  m_view( AXOM_NULLPTR ),
#endif
  m_data( AXOM_NULLPTR ),
  m_num_tuples( num_tuples ),
  m_capacity(0),
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
  SLIC_ERROR_IF( m_num_tuples > capacity,
                 "Number of tuples (" << m_num_tuples << ") " <<
                 "cannot be greater than the tuple capacity " <<
                 "(" << capacity << ")." );

  setCapacity( capacity );
}

//------------------------------------------------------------------------------
template< typename T >
Array< T >::Array( T* data, localIndex num_tuples, localIndex num_components ) :
#ifdef MINT_USE_SIDRE
  m_view( AXOM_NULLPTR ),
#endif
  m_data( data ),
  m_num_tuples( num_tuples ),
  m_capacity( num_tuples ),
  m_num_components( num_components ),
  m_resize_ratio( 0.0 ),
  m_is_external( true )
{
  SLIC_ERROR_IF( m_num_tuples < 0,
                 "Number of tuples (" << m_num_tuples << ") " <<
                 "cannot be negative." );
  SLIC_ERROR_IF( m_num_components <= 0,
                 "Components per tuple (" << m_num_components << ") " <<
                 "must be greater than zero." );
}

#ifdef MINT_USE_SIDRE
//------------------------------------------------------------------------------
template< typename T >
Array< T >::Array( sidre::View* view, localIndex num_tuples ) :
  m_view( view ),
  m_data( AXOM_NULLPTR ),
  m_num_tuples( num_tuples ),
  m_capacity(0),
  m_num_components(0),
  m_resize_ratio( DEFAULT_RESIZE_RATIO ),
  m_is_external( false )
{
  SLIC_ERROR_IF( m_view == AXOM_NULLPTR, "Provided View cannot be null." );
  SLIC_ERROR_IF( m_view->isEmpty(), "View is empty." );

  m_capacity = getViewShape(0);
  m_num_components = getViewShape(1);

  SLIC_ERROR_IF( m_num_tuples < 0,
                 "Number of tuples (" << m_num_tuples << ") " <<
                 "cannot be negative." );
  SLIC_ERROR_IF( m_num_tuples > m_capacity,
                 "Number of tuples (" << m_num_tuples << ") " <<
                 "cannot be greater than the tuple capacity " <<
                 "(" << m_capacity << ")." );

  sidre::TypeID view_type = m_view->getTypeID();
  sidre::TypeID T_type = sidre::detail::SidreTT< T >::id;
  SLIC_ERROR_IF( view_type != T_type, "View data type (" << view_type << ")"
                                                         << "differs from this Array type (" << T_type <<
      ")." );

  m_data = static_cast< T* >( m_view->getVoidPtr() );
}

//------------------------------------------------------------------------------
template< typename T >
Array< T >::Array( sidre::View* view, localIndex num_tuples,
                   localIndex num_components ) :
  m_view( view ),
  m_data( AXOM_NULLPTR ),
  m_num_tuples( num_tuples ),
  m_capacity( num_tuples * DEFAULT_RESIZE_RATIO ),
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
                 "must be greater than zero." );

  sidre::TypeID T_type = sidre::detail::SidreTT< T >::id;
  sidre::SidreLength dims[ 2 ];
  dims[0] = m_capacity;
  dims[1] = m_num_components;
  m_view->allocate( T_type, m_capacity * m_num_components );
  m_view->apply( T_type, 2, dims );
  m_data = m_view->getVoidPtr();
}

//------------------------------------------------------------------------------
template< typename T >
Array< T >::Array( sidre::View* view, localIndex capacity,
                   localIndex num_tuples, localIndex num_components ) :
  m_view( view ),
  m_data( AXOM_NULLPTR ),
  m_num_tuples( num_tuples ),
  m_capacity( capacity ),
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
                 "must be greater than zero." );
  SLIC_ERROR_IF( m_num_tuples > m_capacity,
                 "Number of tuples (" << m_num_tuples << ") " <<
                 "cannot be greater than the tuple capacity " <<
                 "(" << m_capacity << ")." );

  sidre::TypeID T_type = sidre::detail::SidreTT< T >::id;
  sidre::SidreLength dims[ 2 ];
  dims[0] = m_capacity;
  dims[1] = m_num_components;
  m_view->allocate( T_type, m_capacity * m_num_components );
  m_view->apply( T_type, 2, dims );
  m_data = static_cast< T* >( m_view->getVoidPtr() );
}
#endif

//------------------------------------------------------------------------------
template< typename T >
Array< T >::~Array()
{
#ifdef MINT_USE_SIDRE
  if ( m_view == AXOM_NULLPTR && m_data != AXOM_NULLPTR )
  {
    utilities::free( m_data );
  }
#else
  if ( m_data != AXOM_NULLPTR )
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

  localIndex new_size = m_num_tuples + 1;
  if ( new_size > m_capacity )
  {
    dynamicRealloc( new_size );
  }

  m_data[ m_num_tuples ] = value;
  m_num_tuples = new_size;
}

//------------------------------------------------------------------------------
template< typename T >
inline void Array< T >::append( const T* tuples, localIndex n )
{
  localIndex new_size = m_num_tuples + n;
  if ( new_size > m_capacity )
  {
    dynamicRealloc( new_size );
  }

  T* cur_end = m_data + m_num_tuples * m_num_components;
  std::memcpy( cur_end, tuples, n * m_num_components * sizeof(T) );
  m_num_tuples = new_size;
}

//------------------------------------------------------------------------------
template< typename T >
inline void Array< T >::set( const T* tuples, localIndex n, localIndex pos )
{
  SLIC_ASSERT( pos >= 0 );
  SLIC_ASSERT( pos + n <= m_num_tuples );
  T* set_position = m_data + pos * m_num_tuples;
  localIndex byte_size = n * m_num_tuples * sizeof(T);
  std::memcpy( set_position, tuples, byte_size );
}

//------------------------------------------------------------------------------
template< typename T >
inline void Array< T >::insert( const T* tuples, localIndex n, localIndex pos )
{
  T* insert_pos = reserveForInsert( n, pos );
  std::memcpy( insert_pos, tuples, n * m_num_components * sizeof(T) );
}

//------------------------------------------------------------------------------
template< typename T >
inline void Array< T >::insert( const T& value, localIndex pos )
{
  SLIC_ASSERT_MSG( m_num_components != 1, "Number of components must be 1." );
  reserveForInsert( 1, pos );
  m_data[ pos ] = value;
}

//------------------------------------------------------------------------------
template< typename T >
inline void Array< T >::resize( localIndex num_tuples )
{
  SLIC_ASSERT( num_tuples >= 0 );
  SLIC_ERROR_IF( m_is_external, "Cannot change the capacity of external data.");

  m_num_tuples = num_tuples;
  if ( m_num_tuples > m_capacity )
  {
    setCapacity( m_resize_ratio * m_num_tuples );
  }
}

//------------------------------------------------------------------------------
template< typename T >
inline void Array< T >::setCapacity( localIndex capacity )
{
  SLIC_ASSERT( capacity >= 0 );
  SLIC_ERROR_IF( m_is_external, "Cannot change the capacity of external data.");

  m_capacity = capacity;
  if ( m_capacity < m_num_tuples )
  {
    m_num_tuples = m_capacity;
  }

#ifdef MINT_USE_SIDRE
  if ( m_view != AXOM_NULLPTR )
  {
    return allocate_view_data();
  }
#endif

  m_data = utilities::realloc( m_data, m_capacity * m_num_components );
  SLIC_ERROR_IF( m_data == AXOM_NULLPTR && m_capacity > 0,
                 "Array reallocation failed." );
}

//------------------------------------------------------------------------------
template< typename T >
inline void Array< T >::dynamicRealloc( localIndex new_num_tuples )
{
  SLIC_ERROR_IF( m_is_external, "Cannot change the capacity of external data.");
  SLIC_ERROR_IF( m_resize_ratio < 1.0, "Resize ratio of " << m_resize_ratio <<
                 " doesn't support dynamic resizing");
  m_capacity = new_num_tuples * m_resize_ratio + 0.5;

#ifdef MINT_USE_SIDRE
  if ( m_view != AXOM_NULLPTR )
  {
    return allocate_view_data();
  }
#endif

  m_data = utilities::realloc( m_data, m_capacity * m_num_components );
  SLIC_ERROR_IF( m_data == AXOM_NULLPTR && m_capacity > 0,
                 "Array reallocation failed." );
}

//------------------------------------------------------------------------------
template< typename T >
inline T* Array< T >::reserveForInsert( localIndex n, localIndex pos )
{
  SLIC_ASSERT( n > 0 );
  SLIC_ASSERT( pos >= 0 );
  SLIC_ASSERT( pos <= m_num_tuples );

  localIndex new_size = m_num_tuples + n;
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

  m_num_tuples = new_size;
  return insert_pos;
}

#ifdef MINT_USE_SIDRE
//------------------------------------------------------------------------------
template< typename T >
inline localIndex Array< T >::getViewShape( int dim ) const
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
inline void Array< T >::allocate_view_data()
{
  m_view->reallocate( m_capacity * m_num_components );

  sidre::TypeID T_type = sidre::detail::SidreTT< T >::id;
  sidre::SidreLength dims[ 2 ];
  dims[0] = m_capacity;
  dims[1] = m_num_components;
  m_view->apply(T_type, 2, dims );
  m_data = static_cast< T* >( m_view->getVoidPtr() );

  SLIC_ERROR_IF( m_data == AXOM_NULLPTR && m_capacity > 0,
                 "Array reallocation failed." );
}
#endif


} /* namespace mint */
} /* namespace axom */

#endif /* MINT_UTILS_ARRAY_HXX_ */

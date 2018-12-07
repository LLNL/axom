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

#ifndef SIDRE_ARRAY_HPP_
#define SIDRE_ARRAY_HPP_

#include "axom/core/Macros.hpp"      // for disable copy/assignment macro
#include "axom/core/utilities/Utilities.hpp"  // for memory allocation functions
#include "axom/core/utilities/Array.hpp"  // to inherit
#include "axom/core/Types.hpp"

#include "axom/slic/interface/slic.hpp"            // for slic logging macros

#include "View.hpp"          // for View definition
#include "Buffer.hpp"        // for Buffer definition

// C/C++ includes
#include <cstring>                  // for std::memcpy

namespace axom
{
namespace sidre
{

/* Provided so that 0 doesn't convert to nullptr and lead to ambiguous
 * constructor calls. */
namespace internal
{
constexpr axom::IndexType ZERO = 0;
}

/*!
 * \class Array
 *
 * \brief Provides a generic multi-component array, contained in Sidre.
 *
 *  The Array class provides a generic multi-component array container with
 *  dynamic re-allocation and insertion. Each element in the array is a tuple
 *  consisting of 1 or more components, which are stored contiguously.
 *
 *  Depending on which constructor is used, the Array object can have three
 *  different underlying storage types:
 *
 *  * <b> Sidre </b> <br />
 *
 *    An Array object may also be constructed from a View. From the
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
 * \see Group
 * \see View
 */
template< typename T >
class Array : public axom::utilities::Array<T>
{

public:

  /*!
   * \brief Default constructor. Disabled.
   */
  Array() = delete;

/// \name Sidre Array constructors
/// @{

  /*!
   * \brief Creates an Array instance from a View that already has data.
   *
   * \param [in] view the View that holds this Array's data.
   *
   * \note The Sidre view shape has two dimensions. The first dimension
   *  corresponds to the max capacity of the array and the second corresponds to
   *  the number of components per tuple.
   *
   * \pre view != nullptr
   * \pre view->isEmpty() == false.
   * \pre view->getNumDimensions() == 2
   *
   * \post capacity() == view->getDimension(0)
   * \post numComponents() == view->getDimension(1)
   * \post getResizeRatio() == DEFAULT_RESIZE_RATIO
   */
  Array( View* view );

  /*!
   * \brief Creates an Array instance of `num_tuples` size, where each
   *  tuple consists of `num_components` values and populates the associated
   *  View.
   *
   * \param [in] view the View that will hold this Array's data.
   * \param [in] num_tuples the number of tuples accounted for in the Array.
   * \param [in] num_components the number of values per tuple. If not
   *  specified defaults to 1.
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
   * \pre view != nullptr
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
  Array( View* view, axom::IndexType num_tuples, axom::IndexType num_components=1,
         axom::IndexType capacity=USE_DEFAULT );

/// @}

  /*!
   * Destructor. Free's the associated buffer unless the memory is external
   * or owned by Sidre.
   */
  virtual ~Array();

/// \name Array methods to query and set attributes
/// @{

  /*!
   * \brief Return true iff the external buffer constructor was called.
   */
  virtual bool isExternal() const { return false; }

  /*!
   * \brief Return true iff a sidre constructor was called.
   */
  virtual bool isInSidre() const
  {
    return m_view != nullptr;
  }

  /*!
   * \brief Return a pointer to the View that this Array wraps.
   */
  const View* getView() const { return m_view; }

/// @}

protected:

  /*!
   * \brief Update the number of tuples.
   *
   * \param [in] new_num_tuples the new number of tuples.
   */
  virtual void updateNumTuples( axom::IndexType new_num_tuples );

  /*!
   * \brief Set the number of tuples allocated for the data array.
   *
   * \param [in] capacity the new number of tuples to allocate.
   */
  virtual void setCapacity( axom::IndexType new_capacity );

  /*!
   * \brief Reallocates the data array when the size exceeds the capacity.
   *
   * \param [in] new_num_tuples the number of tuples which exceeds the current
   *  capacity.
   */
  virtual void dynamicRealloc( axom::IndexType new_num_tuples );


  /*!
   * \brief Return the TypeID corresponding to T. This function
   *  handles when T is not an enum.
   */
  template < typename U = T >
  static constexpr
  typename std::enable_if< !std::is_enum< U >::value, TypeID >::type
  sidreTypeId()
  { return detail::SidreTT< U >::id; }

  /*!
   * \brief Return the TypeID corresponding to T. This function handles
   *  when T is an enum.
   */
  template < typename U = T >
  static constexpr
  typename std::enable_if< std::is_enum< U >::value, TypeID >::type
  sidreTypeId()
  {
    return detail::SidreTT<
      typename std::underlying_type< U >::type >::id;
  }

  /*!
   * \brief Describes m_view as having dimensions
   * (m_num_tuples, m_num_components).
   */
  void describeView();

  /*!
   * \brief Given a non-empty View of dimension 2, returns the length
   *  of the given dimension.
   *
   * \param [in] view the View to examine.
   * \param [in] dim the dimension (0 or 1) to return the length of.
   *
   * \pre 0 <= dim <= 1
   */
  axom::IndexType getViewShape( int dim ) const;

  /*!
   * \brief Allocates space within the Array's View.
   */
  void reallocViewData();

  View* m_view;

  DISABLE_COPY_AND_ASSIGNMENT( Array );
  DISABLE_MOVE_AND_ASSIGNMENT( Array );
};


//------------------------------------------------------------------------------
//                            Array IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template< typename T >
Array< T >::Array( View* view ) :
  utilities::Array<T>(0, 0, 0),
  m_view( view )
{
  SLIC_ERROR_IF( m_view == nullptr, "Provided View cannot be null." );
  SLIC_ERROR_IF( m_view->isEmpty(), "Provided View cannot be empty." );

  this->m_num_tuples       = getViewShape(0);
  this->m_num_components = getViewShape(1);

  axom::IndexType buffer_size = m_view->getBuffer()->getNumElements();
  SLIC_ERROR_IF( buffer_size % this->m_num_components != 0,
                 "The buffer size (" << buffer_size << ") " <<
                 "is not a multiple of the number of components " <<
                 "(" << this->m_num_components << ")." );
  this->m_capacity = buffer_size / this->m_num_components;

  SLIC_ERROR_IF( this->m_num_tuples < 0,
                 "Number of tuples (" << this->m_num_tuples << ") " <<
                 "cannot be negative." );

  SLIC_ERROR_IF( this->m_num_components <= 0,
                 "Number of components (" << this->m_num_components << ") " <<
                 "must be greater than 0." );

  SLIC_ERROR_IF( this->m_num_tuples > this->m_capacity,
                 "Number of tuples (" << this->m_num_tuples << ") " <<
                 "cannot be greater than the tuple capacity " <<
                 "(" << this->m_capacity << ")." );

  TypeID view_type = m_view->getTypeID();
  TypeID T_type = sidreTypeId();
  SLIC_ERROR_IF( view_type != T_type,
                 "View data type (" << view_type << ")" <<
                 "differs from this Array type (" << T_type << ")." );

  this->m_data = static_cast< T* >( m_view->getVoidPtr() );
  SLIC_ERROR_IF( this->m_data == nullptr && this->m_capacity > 0,
                 "View returned a null pointer when the capacity " <<
                 "is greater than zero." );
}

//------------------------------------------------------------------------------
template< typename T >
Array< T >::Array( View* view, axom::IndexType num_tuples,
                   axom::IndexType num_components, axom::IndexType capacity ) :
  utilities::Array<T>(num_tuples, num_components, capacity),
  m_view( view )
{
  SLIC_ERROR_IF( m_view == nullptr, "Provided View cannot be null." );
  SLIC_ERROR_IF( !m_view->isEmpty(), "View must be empty." );
  SLIC_ERROR_IF( this->m_num_tuples < 0,
                 "Number of tuples (" << this->m_num_tuples << ") " <<
                 "cannot be negative." );
  SLIC_ERROR_IF( this->m_num_components <= 0,
                 "Components per tuple (" << this->m_num_components << ") " <<
                 "must be greater than 0." );

  if ( capacity == USE_DEFAULT )
  {
    capacity =
      ( this->m_num_tuples > utilities::Array<T>::MIN_DEFAULT_CAPACITY ) ?
      this->m_num_tuples : utilities::Array<T>::MIN_DEFAULT_CAPACITY;
  }
  SLIC_ERROR_IF( this->m_num_tuples > capacity,
                 "Number of tuples (" << this->m_num_tuples << ") " <<
                 "cannot be greater than the tuple capacity " <<
                 "(" << capacity << ")." );

  setCapacity( capacity );

  // sanity checks
  SLIC_ASSERT( capacity >= 0 );
  if ( capacity > 0 )
  {
    SLIC_ASSERT( this->m_data != nullptr );
  }
  SLIC_ASSERT( this->m_num_tuples >= 0 );
  SLIC_ASSERT( this->m_num_components >= 1 );
}

//------------------------------------------------------------------------------
template< typename T >
Array< T >::~Array()
{
  if ( m_view == nullptr )
  {
    utilities::free( this->m_data );
  }
  m_view = nullptr;
  this->m_data = nullptr;
}


//------------------------------------------------------------------------------
template< typename T >
inline void Array< T >::updateNumTuples( axom::IndexType new_num_tuples )
{
  SLIC_ASSERT( new_num_tuples >= 0 );
  SLIC_ASSERT( new_num_tuples <= this->m_capacity );
  this->m_num_tuples = new_num_tuples;

  if ( m_view != nullptr )
  {
    describeView();
  }
}

//------------------------------------------------------------------------------
template< typename T >
inline void Array< T >::setCapacity( axom::IndexType new_capacity )
{
  SLIC_ASSERT( new_capacity >= 0 );

  this->m_capacity = new_capacity;
  if ( this->m_capacity < this->m_num_tuples )
  {
    updateNumTuples( this->m_capacity );
  }

  if ( m_view != nullptr )
  {
    return reallocViewData();
  }

  this->m_data =
    utilities::realloc( this->m_data,
                        this->m_capacity * this->m_num_components );
  SLIC_ERROR_IF( this->m_data == nullptr && this->m_capacity > 0,
                 "Array reallocation failed." );
}

//------------------------------------------------------------------------------
template< typename T >
inline void Array< T >::dynamicRealloc( axom::IndexType new_num_tuples )
{
  SLIC_ERROR_IF( this->m_is_external,
                 "Cannot change the capacity of external data.");
  SLIC_ERROR_IF( this->m_resize_ratio < 1.0, "Resize ratio of " <<
                 this->m_resize_ratio << " doesn't support dynamic resizing");
  this->m_capacity = new_num_tuples * this->m_resize_ratio + 0.5;

  if ( m_view != nullptr )
  {
    return reallocViewData();
  }

  this->m_data =
    utilities::realloc( this->m_data,
                        this->m_capacity * this->m_num_components );
  SLIC_ERROR_IF( this->m_data == nullptr && this->m_capacity > 0,
                 "Array reallocation failed." );
}

//------------------------------------------------------------------------------
template< typename T >
inline void Array< T >::describeView()
{
  SLIC_ASSERT( m_view != nullptr );

  static constexpr TypeID T_type = sidreTypeId();
  IndexType dims[2];
  dims[0] = this->m_num_tuples;
  dims[1] = this->m_num_components;

  m_view->apply( T_type, 2, dims );
}

//------------------------------------------------------------------------------
template< typename T >
inline axom::IndexType Array< T >::getViewShape( int dim ) const
{
  SLIC_ERROR_IF( dim > 1, "Only two dimensional views supported." );
  SLIC_ERROR_IF( m_view->isEmpty(), "view cannot be empty." );
  SLIC_ERROR_IF( m_view->getNumDimensions() != 2,
                 "view must have dimension 2.");

  sidre::IndexType dims[ 2 ];
  m_view->getShape( 2, dims );
  return static_cast<axom::IndexType>(dims[ dim ]);
}

//------------------------------------------------------------------------------
template< typename T >
inline void Array< T >::reallocViewData()
{
  if ( m_view->isEmpty() )
  {
    constexpr sidre::TypeID T_type = sidreTypeId();
    m_view->allocate( T_type, this->m_capacity * this->m_num_components );
  }
  else
  {
    m_view->reallocate( this->m_capacity * this->m_num_components );
  }

  describeView();
  this->m_data = static_cast< T* >( m_view->getVoidPtr() );

  SLIC_ERROR_IF( this->m_data == nullptr && this->m_capacity > 0,
                 "Array reallocation failed." );
}

} /* namespace sidre */
} /* namespace axom */

#endif /* SIDRE_ARRAY_HPP_ */

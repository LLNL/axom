// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/Array.hpp"                /* for axom::Array */
#include "axom/core/memory_management.hpp"    /* for alloc() and free() */

#include "gtest/gtest.h"                      /* for TEST and EXPECT_* macros */

// C/C++ includes
#include <algorithm>                          /* for std::fill_n */

namespace axom
{

const char IGNORE_OUTPUT[] = ".*";

namespace internal
{

/*!
 * \brief Calculate the new capacity for an Array given an increase in the
 *  size.
 * \param [in] v, the Array in question.
 * \param [in] increase, the amount the size will increase by
 * \return the new capacity.
 */
template < typename T >
IndexType calc_new_capacity( Array< T > & v, IndexType increase )
{
  IndexType new_num_elements = v.size() + increase;
  if ( new_num_elements > v.capacity() )
  {
    return new_num_elements * v.getResizeRatio() + 0.5;
  }

  return v.capacity();
}

/*!
 * \brief Check if two Arrays are copies. Does not check the resize ratio.
 * \param [in] lhs, the first Array to compare.
 * \param [in] rhs, the second Array to compare.
 * \return the new capacity.
 */
template< typename T >
void check_copy( const Array< T >& lhs, const Array< T >& rhs )
{
  EXPECT_EQ( lhs.size(), rhs.size() );
  EXPECT_EQ( lhs.capacity(), rhs.capacity() );

  const T* lhs_data = lhs.data();
  const T* rhs_data = rhs.data();
  EXPECT_EQ( lhs_data, rhs_data );
}

/*!
 * \brief Check that the storage of an Array is working properly.
 * \param [in] v the Array to check.
 */
template < typename T >
void check_storage( Array< T >& v )
{
  EXPECT_TRUE( v.empty() );
  EXPECT_EQ( v.size(), 0 );

  IndexType capacity = v.capacity();
  const T* data_ptr = v.data();

  /* Push back up to half the capacity. */
  for ( T i = 0 ; i < capacity / 2 ; ++i )
  {
    v.push_back( i );
  }

  /* Check the array metadata. */
  EXPECT_TRUE( !v.empty() );
  EXPECT_EQ( v.size(), capacity / 2 );
  EXPECT_EQ( v.capacity(), capacity );
  EXPECT_EQ( v.data(), data_ptr );

  /* Push back up to the full capacity. */
  for ( T i = capacity / 2 ; i < capacity ; ++i )
  {
    v.push_back( i );
  }

  /* Check the array metadata. */
  EXPECT_TRUE( !v.empty() );
  EXPECT_EQ( v.size(), capacity );
  EXPECT_EQ( v.capacity(), capacity );
  EXPECT_EQ( v.data(), data_ptr );

  /* Check the array data using the [] operator and the raw pointer. */
  for ( IndexType i = 0 ; i < capacity ; ++i )
  {
    EXPECT_EQ( v[ i ], i );
    EXPECT_EQ( data_ptr[ i ], i );
  }

  /* Set the array data to new values using the [] operator. */
  for ( IndexType i = 0 ; i < capacity ; ++i )
  {
    v[ i ] = i - 5 * i + 7 ;
  }

  /* Check the array data using the [] operator and the raw pointer. */
  for ( IndexType i = 0 ; i < capacity ; ++i )
  {
    EXPECT_EQ( v[ i ], i - 5 * i + 7 );
    EXPECT_EQ( data_ptr[ i ], i - 5 * i + 7 );
  }

  /* Check the array metadata. */
  EXPECT_TRUE( !v.empty() );
  EXPECT_EQ( v.size(), capacity );
  EXPECT_EQ( v.capacity(), capacity );
  EXPECT_EQ( v.data(), data_ptr );
}

/*!
 * \brief Check that the fill method is working properly.
 * \param [in] v the Array to check.
 */
template < typename T >
void check_fill( Array< T >& v )
{
  constexpr T MAGIC_NUM_0 = 55;
  constexpr T MAGIC_NUM_1 = 6834;
  const IndexType capacity = v.capacity();
  const IndexType size = v.size();
  const double ratio = v.getResizeRatio();
  const T* const data_ptr = v.data();

  /* Fill the Array with MAGIC_NUM_0. */
  v.fill( MAGIC_NUM_0 );

  /* Check the meta data. */
  EXPECT_EQ( capacity, v.capacity() );
  EXPECT_EQ( size, v.size() );
  EXPECT_EQ( ratio, v.getResizeRatio() );
  EXPECT_EQ( data_ptr, v.data() );

  /* Check that the entries are all MAGIC_NUM_0. */
  for ( IndexType i = 0 ; i < size ; ++i )
  {
    EXPECT_EQ( v[ i ], MAGIC_NUM_0 );
  }

  /* Fill the Array with MAGIC_NUM_1. */
  v.fill( MAGIC_NUM_1 );

  /* Check the meta data. */
  EXPECT_EQ( capacity, v.capacity() );
  EXPECT_EQ( size, v.size() );
  EXPECT_EQ( ratio, v.getResizeRatio() );
  EXPECT_EQ( data_ptr, v.data() );

  /* Check that the entries are all MAGIC_NUM_1. */
  for ( IndexType i = 0 ; i < size ; ++i )
  {
    EXPECT_EQ( v[ i ], MAGIC_NUM_1 );
  }
}

/*!
 * \brief Check that the set method is working properly.
 * \param [in] v the Array to check.
 */
template < typename T >
void check_set( Array< T >& v )
{
  constexpr T ZERO = 0;
  const IndexType capacity = v.capacity();
  const IndexType size = v.size();
  const double ratio = v.getResizeRatio();
  const T* const data_ptr = v.data();

  /* Allocate a buffer half the size of the array. Fill it up with sequential
   * values. */
  const IndexType buffer_size = size / 2;
  T* buffer = allocate<T>( buffer_size );
  for ( IndexType i = 0 ; i < buffer_size ; ++i )
  {
    buffer[ i ] = i;
  }

  /* Set all the values in the array to zero. */
  v.fill( ZERO );

  /* Set the first half of the elements in the array to the sequential values in
   * buffer. */
  v.set( buffer, buffer_size, 0 );

  /* Check the array metadata. */
  EXPECT_EQ( capacity, v.capacity() );
  EXPECT_EQ( size, v.size() );
  EXPECT_EQ( ratio, v.getResizeRatio() );
  EXPECT_EQ( data_ptr, v.data() );

  /* Check that the first half of the elements in the array are equivalent to
   * those in buffer. */
  for ( IndexType i = 0 ; i < buffer_size ; ++i )
  {
    EXPECT_EQ( v[ i ], buffer[ i ] );
  }

  /* Check that the second half of the elements in the array are all zero. */
  for ( IndexType i = buffer_size ; i < size ; ++i )
  {
    EXPECT_EQ( v[ i ], ZERO );
  }

  /* Reset the values in buffer to the next sequential values. */
  for ( IndexType i = 0 ; i < buffer_size ; ++i )
  {
    buffer[ i ] = i + buffer_size;
  }

  /* Set the second half of the elements in the array to the new sequential
   * values in buffer. */
  v.set( buffer, buffer_size, buffer_size );

  /* Check the array metadata. */
  EXPECT_EQ( capacity, v.capacity() );
  EXPECT_EQ( size, v.size() );
  EXPECT_EQ( ratio, v.getResizeRatio() );
  EXPECT_EQ( data_ptr, v.data() );

  /* Check that all the elements in the array now hold sequential values. */
  for ( IndexType i = 0 ; i < 2 * buffer_size ; ++i )
  {
      EXPECT_EQ( v[ i ], i );
  }

  deallocate( buffer );
}

/*!
 * \brief Check that the resizing of an Array is working properly.
 * \param [in] v the Array to check.
 */
template < typename T >
void check_resize( Array< T >& v )
{
  /* Resize the array up to the capacity */
  IndexType capacity = v.capacity();
  v.resize( capacity );
  IndexType size = capacity;

  /* Check that the size equals the capacity. */
  EXPECT_EQ( v.size(), v.capacity() );

  /* Set the existing data in v */
  for ( IndexType i = 0 ; i < size ; ++i )
  {
    v[ i ] = static_cast<T>( i - 5 * i + 7 );
  }

  /* Push back a new element, should resize. */
  IndexType old_capacity = capacity;
  capacity = calc_new_capacity( v, 1 );
  v.push_back( size - 5 * size + 7 );
  size++;

  /* Check that it resized properly */
  EXPECT_GT( capacity, old_capacity );
  EXPECT_EQ( v.capacity(), capacity );
  EXPECT_EQ( v.size(), size );

  /* Check that the data is still intact. */
  for ( IndexType i = 0 ; i < size ; ++i )
  {
      EXPECT_EQ( v[ i ], i - 5 * i + 7 );
  }

  /* Push back 1000 elements */
  const IndexType n_elements = 1000;

  T* values = allocate< T >( n_elements );
  for ( IndexType i = 0 ; i < n_elements ; ++i )
  {
      IndexType i_real = i + size;
      values[ i ] = i_real - 5 * i_real + 7;
  }

  /* Push back the new elements. */
  capacity = calc_new_capacity( v, n_elements );
  v.insert(size, n_elements, values);
  size += n_elements;

  /* Check that size and capacity are as expected. */
  EXPECT_EQ( v.capacity(), capacity );
  EXPECT_EQ( v.size(), size );

  /* Check that the data is still intact. */
  for ( IndexType i = 0 ; i < size ; ++i )
  {
    EXPECT_EQ( v[ i ], i - 5 * i + 7 );
  }

  /* Reduce the size down to 500 elements */
  T* data_address = v.data();
  size = 500;
  v.resize(size);

  /* Check the metadata. */
  EXPECT_EQ( v.size(), size );
  EXPECT_EQ( v.capacity(), capacity );
  EXPECT_EQ( v.data(), data_address );

  /* Check the data. */
  for ( IndexType i = 0 ; i < size ; ++i )
  {
    EXPECT_EQ( v[ i ], i - 5 * i + 7 );
  }

  /* Shrink the vector */
  capacity = size;
  v.shrink();

  /* Check that the capacity and size are as expected. */
  EXPECT_EQ( v.capacity(), capacity );
  EXPECT_EQ( v.size(), size );

  /* Check that the data is intact. */
  for ( IndexType i = 0 ; i < size ; ++i )
  {
    EXPECT_EQ( v[ i ], i - 5 * i + 7 );
  }

  /* Push back a new element, should resize. */
  old_capacity = capacity;
  capacity = calc_new_capacity( v, 1 );
  v.push_back( size - 5 * size + 7 );
  size++;

  /* Check the new size and capacity. */
  EXPECT_GT( capacity, old_capacity );
  EXPECT_EQ( v.capacity(), capacity );
  EXPECT_EQ( v.size(), size );

  /* Check that the data is intact. */
  for ( IndexType i = 0 ; i < size ; ++i )
  {
    EXPECT_EQ( v[ i ], i - 5 * i + 7 );
  }

  /* Reset the data */
  T* data_ptr = v.data();
  for ( IndexType i = 0 ; i < size ; ++i )
  {
    data_ptr[ i ] = i;
  }

  /* Push back a bunch of elements to fill in up to the capacity. Resize should
   * not occur. */
  old_capacity = capacity;
  for ( IndexType i = size ; i < old_capacity ; ++i )
  {

    v.push_back( i );
    size++;
    EXPECT_EQ( v.capacity(), old_capacity );
    EXPECT_EQ( v.size(), size );
    EXPECT_EQ( v.data(), data_ptr );
  }

  EXPECT_EQ( v.size(), old_capacity );

  /* Push back a final element that should trigger a resize. */
  capacity = calc_new_capacity( v, old_capacity - size + 1 );
  v.push_back( size );
  size++;

  /* Check the new capacity and size. */
  EXPECT_GT( capacity, old_capacity );
  EXPECT_EQ( v.capacity(), capacity );
  EXPECT_EQ( v.size(), size );

  /* Check the data. */
  for ( IndexType i = 0 ; i < size ; ++i )
  {
    EXPECT_EQ( v[ i ], i );
  }

  deallocate( values );
  values = nullptr;
}

/*!
 * \brief Check that the insertion into an Array is working properly.
 * \param [in] v the Array to check.
 */
template < typename T >
void check_insert( Array< T >& v )
{
  /* Resize the array up to the capacity */
  IndexType capacity = v.capacity();
  v.resize( capacity );
  IndexType size = capacity;

  EXPECT_EQ( v.size(), v.capacity() );

  /* Set the existing data in v */
  for ( IndexType i = 0 ; i < size ; ++i )
  {
      v[ i ] = i - 5 * i + 7 ;
  }

  /* Insert a new element, should resize. */
  IndexType old_capacity = capacity;
  capacity = calc_new_capacity( v, 1 );
  v.insert( v.size(), 1, size - 5 * size + 7);
  size++;

  /* Check that it resized properly */
  EXPECT_GT( capacity, old_capacity );
  EXPECT_EQ( v.capacity(), capacity );
  EXPECT_EQ( v.size(), size );
  for ( IndexType i = 0 ; i < size ; ++i )
  {
    EXPECT_EQ( v[ i ], i - 5 * i + 7 );
  }

  /* Insert 1000 elements */
  const IndexType n_elements = 1000;
  T* values = allocate< T >( n_elements );
  for ( IndexType i = 0 ; i < n_elements ; ++i )
  {
      IndexType i_real = i + size;
      values[ i ] = i_real - 5 * i_real + 7 ;
  }

  capacity = calc_new_capacity( v, n_elements );
  v.insert( size, n_elements, values );
  size += n_elements;

  /* Check that it resizes properly */
  EXPECT_EQ( v.capacity(), capacity );
  EXPECT_EQ( v.size(), size );
  for ( IndexType i = 0 ; i < size ; ++i )
  {
    EXPECT_EQ( v[ i ], i - 5 * i + 7 );
  }

  capacity = size;
  v.shrink();
  IndexType n_insert_front = 100;

  /* Reset the data */
  T* data_ptr = v.data();
  for ( IndexType i = 0 ; i < size ; ++i )
  {
    data_ptr[ i ] = i + n_insert_front;
  }

  /* Insert into the front of the Array. */
  for ( IndexType i = n_insert_front - 1 ; i >= 0 ; i--)
  {
    capacity = calc_new_capacity( v, 1 );
    v.insert( 0, 1, i );
    size++;
  }

  /* Check that the insertion worked as expected */
  EXPECT_EQ( v.capacity(), capacity );
  EXPECT_EQ( v.size(), size );
  for ( IndexType i = 0 ; i < size ; ++i )
  {
      EXPECT_EQ( v[ i ], i );
  }

  deallocate( values );
  values = nullptr;
}

/*!
 * \brief Check that the insertion into an Array is working properly
 *        for iterators.
 * \param [in] v the Array to check.
 */
template < typename T >
void check_insert_iterator( Array< T >& v )
{
  /* Resize the array up to the capacity */
  IndexType capacity = v.capacity();
  v.resize( capacity );
  IndexType size = capacity;

  EXPECT_EQ( v.size(), v.capacity() );

  /* Set the existing data in v */
  for ( IndexType i = 0 ; i < size ; ++i )
  {
      v[ i ] = i - 5 * i + 7 ;
  }

  /* Insert a new element, should resize. */
  IndexType old_capacity = capacity;
  capacity = calc_new_capacity( v, 1 );
  typename axom::Array< T >::ArrayIterator ret = 
    v.insert( v.end(), 1, size - 5 * size + 7);
  size++;

  /* Check that it resized properly */
  EXPECT_GT( capacity, old_capacity );
  EXPECT_EQ( v.capacity(), capacity );
  EXPECT_EQ( v.size(), size );
  EXPECT_EQ( ret, v.end() - 1);
  for ( IndexType i = 0 ; i < size ; ++i )
  {
    EXPECT_EQ( v[ i ], i - 5 * i + 7 );
  }

  /* Insert 1000 elements */
  const IndexType n_elements = 1000;
  T* values = allocate< T >( n_elements );
  for ( IndexType i = 0 ; i < n_elements ; ++i )
  {
      IndexType i_real = i + size;
      values[ i ] = i_real - 5 * i_real + 7 ;
  }

  capacity = calc_new_capacity( v, n_elements );
  typename axom::Array< T >::ArrayIterator ret2 =
    v.insert( v.end(), n_elements, values );
  
  EXPECT_EQ( ret2, v.begin() + size);
  size += n_elements;

  /* Check that it resizes properly */
  EXPECT_EQ( v.capacity(), capacity );
  EXPECT_EQ( v.size(), size );
  for ( IndexType i = 0 ; i < size ; ++i )
  {
    EXPECT_EQ( v[ i ], i - 5 * i + 7 );
  }

  capacity = size;
  v.shrink();
  IndexType n_insert_front = 100;

  /* Reset the data */
  T* data_ptr = v.data();
  for ( IndexType i = 0 ; i < size ; ++i )
  {
    data_ptr[ i ] = i + n_insert_front;
  }

  /* Insert into the front of the Array. */
  for ( IndexType i = n_insert_front - 1 ; i >= 0 ; i--)
  {
    capacity = calc_new_capacity( v, 1 );
    typename axom::Array< T >::ArrayIterator ret3 =
      v.insert( v.begin(), 1, i );
    EXPECT_EQ(ret3, v.begin());
    size++;
  }

  /* Check that the insertion worked as expected */
  EXPECT_EQ( v.capacity(), capacity );
  EXPECT_EQ( v.size(), size );
  for ( IndexType i = 0 ; i < size ; ++i )
  {
      EXPECT_EQ( v[ i ], i );
  }

  deallocate( values );
  values = nullptr;
}

/*!
 * \brief Check that emplace() into an Array is working properly.
 * \param [in] v the Array to check.
 */
template < typename T >
void check_emplace( Array< T >& v )
{
/* Resize the array up to the capacity */
  IndexType capacity = v.capacity();
  v.resize( capacity );
  IndexType size = capacity;

  EXPECT_EQ( v.size(), v.capacity() );

  /* Set the existing data in v */
  for ( IndexType i = 0 ; i < size ; ++i )
  {
      v[ i ] = i - 5 * i + 7 ;
  }

  /* Emplace a new element, should resize. */
  IndexType old_capacity = capacity;
  capacity = calc_new_capacity( v, 1 );
  typename axom::Array< T >::ArrayIterator ret = 
    v.emplace( v.end(), size - 5 * size + 7);
  size++;

  /* Check that it resized properly */
  EXPECT_GT( capacity, old_capacity );
  EXPECT_EQ( v.size(), size );
  EXPECT_EQ( ret, v.end() - 1);
  for ( IndexType i = 0 ; i < size ; ++i )
  {
    EXPECT_EQ( v[ i ], i - 5 * i + 7 );
  }

  /* Emplace_back 1000 elements */
  const IndexType n_elements = 1000;
  for ( IndexType i = 0 ; i < n_elements ; ++i )
  {
      IndexType i_real = i + size;
      v.emplace_back( i_real - 5 * i_real + 7 );

  }

  size += n_elements;

  /* Check that it resizes properly */
  EXPECT_EQ( v.size(), size );
  for ( IndexType i = 0 ; i < size ; ++i )
  {
    EXPECT_EQ( v[ i ], i - 5 * i + 7 );
  }

  capacity = size;
  v.shrink();
  IndexType n_insert_front = 100;

  /* Reset the data */
  T* data_ptr = v.data();
  for ( IndexType i = 0 ; i < size ; ++i )
  {
    data_ptr[ i ] = i + n_insert_front;
  }

  /* Emplace into the front of the Array. */
  for ( IndexType i = n_insert_front - 1 ; i >= 0 ; i--)
  {
    capacity = calc_new_capacity( v, 1 );
    typename axom::Array< T >::ArrayIterator ret3 =
      v.emplace( v.begin(), i );
    EXPECT_EQ(ret3, v.begin());
    size++;
  }

  /* Check that the emplace worked as expected */
  EXPECT_EQ( v.capacity(), capacity );
  EXPECT_EQ( v.size(), size );
  for ( IndexType i = 0 ; i < size ; ++i )
  {
      EXPECT_EQ( v[ i ], i );
  }
}

template< typename T >
void check_swap( Array< T >& v)
{ 
  axom::Array< T > v_two ( v.size() );

  /* Push 0...size elements */
  for (int i = 0; i < v.size(); i++)
  {
    v[ i ] = i;
    v_two[ i ] = -i;
  }

  /* Create copies */
  axom::Array< T > v_copy (v);
  axom::Array< T > v_two_copy (v_two);

  EXPECT_EQ(v, v_copy);
  EXPECT_EQ(v_two, v_two_copy);

  /* Swap */
  v.swap(v_two);

  EXPECT_EQ(v_two, v_copy);
  EXPECT_EQ(v, v_two_copy);

  /* Swap back */
  v.swap(v_two);

  EXPECT_EQ(v, v_copy);
  EXPECT_EQ(v_two, v_two_copy);
}

/*!
 * \brief Check an external array for defects.
 * \param [in] v the external array to check.
 */
template< typename T >
void check_external( Array< T >& v )
{
  ASSERT_TRUE( v.isExternal() );

  /* Check that the array is full. */
  ASSERT_EQ( v.size(), v.capacity() );

  const IndexType size = v.size();
  const IndexType num_values = size;
  T* const data_ptr = v.data();

  /* Set the elements in the array. */
  for ( IndexType i = 0 ; i < size ; ++i )
  {
      v[ i ] = i;
  }

  /* Check the elements using the raw pointer. */
  for (IndexType i = 0 ; i < num_values ; ++i )
  {
    EXPECT_EQ( data_ptr[ i ], i );
  }

  /* Set the elements using the raw pointer. */
  for ( IndexType i = 0 ; i < size ; ++i )
  {
      data_ptr[ i ] = i * i;
  }

  /* Check the elements using the () operator. */
  for ( IndexType i = 0 ; i < size ; ++i )
  {
      EXPECT_EQ( v[ i ], i * i );
  }

  EXPECT_EQ( size, v.size() );
  EXPECT_EQ( size, v.capacity() );
  EXPECT_EQ( data_ptr, v.data() );

  /* Since the array is full all of the following calls should require a
   * reallocation and cause a fatal error. */
  T val = T();
  EXPECT_DEATH_IF_SUPPORTED( v.push_back( val ), IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( v.insert( 0,1, val ), IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( v.reserve( size + 1 ), IGNORE_OUTPUT );
}

}   /* end namespace internal */

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
TEST( core_array, checkStorage )
{
  constexpr IndexType ZERO = 0;

  for ( IndexType capacity = 2 ; capacity < 512 ; capacity *= 2 )
  {
    Array< int > v_int( ZERO, capacity );
    internal::check_storage( v_int );

    Array< double > v_double( ZERO, capacity );
    internal::check_storage( v_double );
  }
}

//------------------------------------------------------------------------------
TEST( core_array, checkFill )
{
  for ( IndexType capacity = 2 ; capacity < 512 ; capacity *= 2 )
  {
    IndexType size = capacity / 2;
      Array< int > v_int( size, capacity );
      internal::check_fill( v_int );

      Array< double > v_double( size, capacity );
      internal::check_fill( v_double );
  }
}

//------------------------------------------------------------------------------
TEST( core_array, checkSet )
{
  for ( IndexType size = 2 ; size < 512 ; size *= 2 )
  {
      Array< int > v_int( size );
      internal::check_set( v_int );

      Array< double > v_double( size );
      internal::check_set( v_double );
  }
}

//------------------------------------------------------------------------------
TEST( core_array, checkResize )
{
  constexpr IndexType ZERO = 0;

  for ( double ratio = 1.0 ; ratio <= 2.0 ; ratio += 0.5 )
  {
    for ( IndexType capacity = 2 ; capacity <= 512; capacity *= 2 )
    {
        Array< int > v_int( ZERO, capacity );
        v_int.setResizeRatio( ratio );
        internal::check_resize( v_int );

        Array< double > v_double( ZERO, capacity );
        v_double.setResizeRatio( ratio );
        internal::check_resize( v_double );
    }
  }
}

//------------------------------------------------------------------------------
TEST( core_array_DeathTest, checkResize )
{
  /* Resizing isn't allowed with a ratio less than 1.0. */
  Array< int > v_int( axom::internal::ZERO, 100 );
  v_int.setResizeRatio( 0.99 );
  EXPECT_DEATH_IF_SUPPORTED( internal::check_resize( v_int ), IGNORE_OUTPUT );
}

//------------------------------------------------------------------------------
TEST( core_array, checkInsert )
{
  constexpr IndexType ZERO = 0;

  for ( double ratio = 1.0 ; ratio <= 2.0 ; ratio += 0.5 )
  {
    for ( IndexType capacity = 2 ; capacity <= 512 ; capacity *= 2 )
    {
        Array< int > v_int( ZERO, capacity );
        v_int.setResizeRatio( ratio );
        internal::check_insert( v_int );

        Array< double > v_double( ZERO, capacity );
        v_double.setResizeRatio( ratio );
        internal::check_insert( v_double );
    }
  }
}

//------------------------------------------------------------------------------
TEST( core_array, checkInsertIterator )
{
  constexpr IndexType ZERO = 0;

  for ( double ratio = 1.0 ; ratio <= 2.0 ; ratio += 0.5 )
  {
    for ( IndexType capacity = 10 ; capacity <= 512 ; capacity *= 2 )
    {
        Array< int > v_int( ZERO, capacity );
        v_int.setResizeRatio( ratio );
        internal::check_insert_iterator( v_int );

        Array< double > v_double( ZERO, capacity );
        v_double.setResizeRatio( ratio );
        internal::check_insert_iterator( v_double );
    }
  }
}

//------------------------------------------------------------------------------
TEST( core_array, checkEmplace )
{
  constexpr IndexType ZERO = 0;

  for ( double ratio = 1.0 ; ratio <= 2.0 ; ratio += 0.5 )
  {
    for ( IndexType capacity = 10 ; capacity <= 512 ; capacity *= 2 )
    {
        Array< int > v_int( ZERO, capacity );
        v_int.setResizeRatio( ratio );
        internal::check_emplace( v_int );

        Array< double > v_double( ZERO, capacity );
        v_double.setResizeRatio( ratio );
        internal::check_emplace( v_double );
    }
  }
}

//------------------------------------------------------------------------------
TEST( core_array, checkSwap )
{
  for ( IndexType size = 10 ; size <= 512 ; size *= 2 )
  {
    Array< int > v_int( size );
    internal::check_swap( v_int );

    Array< double > v_double( size );
    internal::check_swap( v_double );
  }
}

//------------------------------------------------------------------------------
TEST( core_array_DeathTest, checkExternal )
{
  constexpr double MAGIC_NUM = 5683578.8;
  constexpr IndexType MAX_SIZE = 256;
  constexpr IndexType MAX_VALUES = MAX_SIZE;
  union DataBuffer
  {
    int ints[ MAX_SIZE ];
    double doubles[ MAX_SIZE ];
  };


  DataBuffer buffer;
  std::fill_n( buffer.doubles, MAX_VALUES, MAGIC_NUM );

  for ( IndexType size = 16 ; size <= MAX_SIZE ; size *= 2 )
  {
    Array< int > v_int( buffer.ints, size, size / 16 );
    EXPECT_EQ( v_int.data(), buffer.ints );
    internal::check_external( v_int );

    Array< double > v_double( buffer.doubles, size, size / 16 );
    EXPECT_EQ( v_double.data(), buffer.doubles );
    internal::check_external( v_double );

    /* Set v_double's data to MAGIC_NUM */
    v_double.fill( MAGIC_NUM );

    /* Check that the data still exists in the buffer */
    for ( IndexType i = 0 ; i < MAX_VALUES ; ++i )
    {
      EXPECT_EQ( buffer.doubles[ i ], MAGIC_NUM );
    }
  }
}

//------------------------------------------------------------------------------
TEST( core_array, checkIterator )
{
  constexpr int SIZE = 1000;
  axom::Array< int > v_int( SIZE );

  /* Push 0...999 elements */
  for (int i = 0; i < SIZE; i++)
  {
    v_int[ i ] = i;
  }

  EXPECT_EQ( *v_int.begin(), 0 );
  EXPECT_EQ( *(v_int.end() - 1), SIZE - 1 );
  EXPECT_EQ( v_int.size(), SIZE );

  /* Erase nothing */
  axom::Array<int>::ArrayIterator ret1 = 
    v_int.erase( v_int.begin() + SIZE / 2, v_int.begin() + SIZE / 2);

  EXPECT_EQ( ret1, v_int.begin() + SIZE / 2 );
  EXPECT_EQ( v_int.size(), SIZE );

  /* Erase half the elements */
  axom::Array<int>::ArrayIterator ret2 =
    v_int.erase (v_int.begin(), v_int.begin() + SIZE / 2 );

  EXPECT_EQ( ret2, v_int.begin() );
  EXPECT_EQ( *v_int.begin(), SIZE / 2 );
  EXPECT_EQ( *(v_int.end() - 1), SIZE - 1 );
  EXPECT_EQ( v_int.size(), SIZE / 2 );

  /* Erase first, last elements */
  axom::Array<int>::ArrayIterator ret3 = v_int.erase( v_int.begin() );

  EXPECT_EQ( ret3, v_int.begin() );
  EXPECT_EQ( *v_int.begin(), SIZE / 2 + 1 );

  axom::Array<int>::ArrayIterator ret4 = v_int.erase( v_int.end() - 1);

  EXPECT_EQ( ret4, v_int.end() );
  EXPECT_EQ( *(v_int.end() - 1), SIZE - 2 );

  /* Clear the rest of the array */
  v_int.clear();
  EXPECT_EQ( v_int.size(), 0);
}

//------------------------------------------------------------------------------ 
TEST( core_array, check_move_copy) 
{ 
  constexpr int MAGIC_INT = 255; 
  constexpr double MAGIC_DOUBLE = 5683578.8; 

  for ( IndexType capacity = 2 ; capacity < 512 ; capacity *= 2 ) 
  { 
    IndexType size = capacity; 

    /* Check copy and move semantics for Array of ints */ 
    Array< int > v_int( size, capacity ); 
    v_int.fill(MAGIC_INT); 

    std::vector <int> ints( size , MAGIC_INT ); 
    Array< int > v_int_external( ints.data(), size, capacity ); 

    Array< int > v_int_copy_ctor( v_int ); 
    Array< int > v_int_copy_assign( 0, 0 ); 
    v_int_copy_assign = v_int; 
    EXPECT_EQ( v_int, v_int_copy_ctor ); 
    EXPECT_EQ( v_int, v_int_copy_assign ); 

    Array< int > v_int_external_copy_ctor( v_int_external ); 
    Array< int > v_int_external_copy_assign( 0, 0 ); 
    v_int_external_copy_assign = v_int_external; 
    EXPECT_EQ( v_int_external, v_int_external_copy_ctor ); 
    EXPECT_EQ( v_int_external, v_int_external_copy_assign ); 

    Array< int > v_int_move_assign( 0, 0 ); 
    v_int_move_assign = std::move( v_int_copy_assign ); 
    Array< int > v_int_move_ctor = std::move( v_int_copy_ctor ); 
    EXPECT_EQ( v_int, v_int_move_assign ); 
    EXPECT_EQ( v_int, v_int_move_ctor ); 
    EXPECT_EQ( v_int_copy_assign.data(), nullptr ); 
    EXPECT_EQ( v_int_copy_ctor.data(), nullptr ); 

    /* Check copy and move semantics for Array of doubles */ 
    Array< double > v_double( size, capacity ); 
    v_double.fill( MAGIC_DOUBLE ); 

    std::vector <double> doubles( size , MAGIC_DOUBLE ); 
    Array< double > v_double_external( doubles.data(), size, capacity ); 

    Array< double > v_double_copy_ctor( v_double ); 
    Array< double > v_double_copy_assign( 0, 0 ); 
    v_double_copy_assign = v_double; 
    EXPECT_EQ( v_double, v_double_copy_ctor ); 
    EXPECT_EQ( v_double, v_double_copy_assign ); 

    Array< double > v_double_external_copy_ctor( v_double_external ); 
    Array< double > v_double_external_copy_assign( 0, 0 ); 
    v_double_external_copy_assign = v_double_external; 
    EXPECT_EQ( v_double_external, v_double_external_copy_ctor ); 
    EXPECT_EQ( v_double_external, v_double_external_copy_assign ); 

    Array< double > v_double_move_assign( 0, 0 ); 
    v_double_move_assign = std::move( v_double_copy_assign );       
    Array< double > v_double_move_ctor = std::move( v_double_copy_ctor ); 
    EXPECT_EQ( v_double, v_double_move_assign ); 
    EXPECT_EQ( v_double, v_double_move_ctor ); 
    EXPECT_EQ( v_double_copy_assign.data(), nullptr ); 
    EXPECT_EQ( v_double_copy_ctor.data(), nullptr ); 
  } 
}

} /* end namespace axom */

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}

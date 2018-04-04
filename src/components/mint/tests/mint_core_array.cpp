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
#include "mint/config.hpp"              /* for mint defintions */
#include "mint/Array.hpp"               /* for mint::Array */

#include "axom_utils/Utilities.hpp"     /* for utilities::max */

#include "slic/UnitTestLogger.hpp"      /* for UnitTestLogger */
#include "slic/slic.hpp"                /* for slic macros */

#include "gtest/gtest.h"                /* for TEST and EXPECT_* macros */

#ifdef MINT_USE_SIDRE
#include "sidre/sidre.hpp"
#endif

// C/C++ includes
#include <algorithm>                    /* for std::fill_n */

namespace axom
{
namespace mint
{

const char IGNORE_OUTPUT[] = ".*";

namespace internal
{

/*!
 * \brief Calculate the new capacity for and Array given an increase in the
 *  size.
 * \param [in] v, the Array in question.
 * \param [in] increase, the ammount the size will increase by
 * \return the new capacity.
 */
template < typename T >
IndexType calc_new_capacity( Array< T > & v, IndexType increase )
{
  IndexType new_num_tuples = v.size() + increase;
  if ( new_num_tuples > v.capacity() )
  {
    return new_num_tuples * v.getResizeRatio() + 0.5;
  }

  return v.capacity();
}

/*!
 * \brief Check if two Arrays are equivalent. Does not check the resize ratio.
 * \param [in] lhs, the first Array to compare.
 * \param [in] rhs, the second Array to compare.
 * \return the new capacity.
 */
template< typename T >
void check_equality( const Array< T >& lhs, const Array< T >& rhs )
{
  EXPECT_EQ( lhs.size(), rhs.size() );
  EXPECT_EQ( lhs.numComponents(), rhs.numComponents() );
  EXPECT_EQ( lhs.capacity(), rhs.capacity() );

  const T* lhs_data = lhs.getData();
  const T* rhs_data = rhs.getData();
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
  IndexType num_components = v.numComponents();
  const T* data_ptr = v.getData();

  if ( num_components == 1 )
  {
    for ( T i = 0 ; i < capacity / 2 ; ++i )
    {
      v.append( i );
    }
  }
  else
  {
    T tuple[ num_components ];
    for ( IndexType i = 0 ; i < capacity / 2 ; ++i )
    {
      for ( IndexType j = 0 ; j < num_components ; ++j )
      {
        tuple[ j ] = i * num_components + j;
      }
      v.append( tuple, 1 );
    }
  }

  EXPECT_TRUE( !v.empty() );
  EXPECT_EQ( v.size(), capacity / 2 );
  EXPECT_EQ( v.capacity(), capacity );
  EXPECT_EQ( v.getData(), data_ptr );

  if ( num_components == 1 )
  {
    for ( T i = capacity / 2 ; i < capacity ; ++i )
    {
      v.append( i );
    }
  }
  else
  {
    T tuple[ num_components ];
    for ( IndexType i = capacity / 2 ; i < capacity ; ++i )
    {
      for ( IndexType j = 0 ; j < num_components ; ++j )
      {
        tuple[ j ] = i * num_components + j;
      }
      v.append( tuple, 1 );
    }
  }

  EXPECT_TRUE( !v.empty() );
  EXPECT_EQ( v.size(), capacity );
  EXPECT_EQ( v.capacity(), capacity );
  EXPECT_EQ( v.getData(), data_ptr );

  for ( IndexType i = 0 ; i < capacity ; ++i )
  {
    for ( IndexType j = 0 ; j < num_components ; ++j )
    {
      EXPECT_EQ( v( i, j ), i * num_components + j );
      EXPECT_EQ( data_ptr[ i * num_components + j ], i * num_components + j );
    }
  }

  for ( IndexType i = 0 ; i < capacity ; ++i )
  {
    for ( IndexType j = 0 ; j < num_components ; ++j )
    {
      v( i, j ) = i * j - 5 * i + 7 * j;
    }
  }

  for ( IndexType i = 0 ; i < capacity ; ++i )
  {
    for ( IndexType j = 0 ; j < num_components ; ++j )
    {
      EXPECT_EQ( v( i, j ), i * j - 5 * i + 7 * j );
      EXPECT_EQ( data_ptr[ i * num_components + j ], i * j - 5 * i + 7 * j );
    }
  }

  EXPECT_TRUE( !v.empty() );
  EXPECT_EQ( v.size(), capacity );
  EXPECT_EQ( v.capacity(), capacity );
  EXPECT_EQ( v.getData(), data_ptr );
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
  const IndexType num_components = v.numComponents();
  const double ratio = v.getResizeRatio();
  const T* const data_ptr = v.getData();

  v.fill( MAGIC_NUM_0 );

  EXPECT_EQ( capacity, v.capacity() );
  EXPECT_EQ( size, v.size() );
  EXPECT_EQ( num_components, v.numComponents() );
  EXPECT_EQ( ratio, v.getResizeRatio() );
  EXPECT_EQ( data_ptr, v.getData() );
  for ( IndexType i = 0; i < size; ++i )
  {
    for ( IndexType j = 0; j < num_components; ++j )
    {
      EXPECT_EQ( v( i, j ), MAGIC_NUM_0 );
    }
  }

  v.fill( MAGIC_NUM_1 );

  EXPECT_EQ( capacity, v.capacity() );
  EXPECT_EQ( size, v.size() );
  EXPECT_EQ( num_components, v.numComponents() );
  EXPECT_EQ( ratio, v.getResizeRatio() );
  EXPECT_EQ( data_ptr, v.getData() );
  for ( IndexType i = 0; i < size; ++i )
  {
    for ( IndexType j = 0; j < num_components; ++j )
    {
      EXPECT_EQ( v( i, j ), MAGIC_NUM_1 );
    }
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
  const IndexType num_components = v.numComponents();
  const double ratio = v.getResizeRatio();
  const T* const data_ptr = v.getData();

  const IndexType buffer_size = size / 2;
  T * buffer = new T[ buffer_size * num_components ];

  for ( IndexType i = 0; i < buffer_size * num_components; ++i )
  {
    buffer[ i ] = i; 
  }

  v.fill( ZERO );
  v.set( buffer, buffer_size, 0 );

  EXPECT_EQ( capacity, v.capacity() );
  EXPECT_EQ( size, v.size() );
  EXPECT_EQ( num_components, v.numComponents() );
  EXPECT_EQ( ratio, v.getResizeRatio() );
  EXPECT_EQ( data_ptr, v.getData() );
  for ( IndexType i = 0; i < buffer_size; ++i )
  {
    for ( IndexType j = 0; j < num_components; ++j )
    {
      EXPECT_EQ( v( i, j ), i * num_components + j );
    }
  }
  for ( IndexType i = buffer_size; i < size; ++i )
  {
    for ( IndexType j = 0; j < num_components; ++j )
    {
      EXPECT_EQ( v( i, j ), ZERO );
    }
  }

  for ( IndexType i = 0; i < buffer_size * num_components; ++i )
  {
    buffer[ i ] = i + buffer_size * num_components; 
  }

  v.set( buffer, buffer_size, buffer_size );

  EXPECT_EQ( capacity, v.capacity() );
  EXPECT_EQ( size, v.size() );
  EXPECT_EQ( num_components, v.numComponents() );
  EXPECT_EQ( ratio, v.getResizeRatio() );
  EXPECT_EQ( data_ptr, v.getData() );
  for ( IndexType i = 0; i < 2 * buffer_size; ++i )
  {
    for ( IndexType j = 0; j < num_components; ++j )
    {
      EXPECT_EQ( v( i, j ), i * num_components + j );
    }
  }
  for ( IndexType i = 2 * buffer_size; i < size; ++i )
  {
    for ( IndexType j = 0; j < num_components; ++j )
    {
      EXPECT_EQ( v( i, j ), ZERO );
    }
  }
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
  IndexType num_components = v.numComponents();

  EXPECT_EQ( v.size(), v.capacity() );

  /* Set the existing data in v */
  for ( IndexType i = 0 ; i < size ; ++i )
  {
    for ( IndexType j = 0 ; j < num_components ; ++j )
    {
      v( i, j ) = i * j - 5 * i + 7 * j;
    }
  }

  /* Append a new tuple, should resize. */
  IndexType old_capacity = capacity;
  capacity = calc_new_capacity( v, 1 );
  T tuple[ num_components ];
  for ( IndexType j = 0 ; j < num_components ; ++j )
  {
    tuple[ j ] = size * j - 5 * size + 7 * j;
  }
  v.append( tuple, 1 );
  size++;

  /* Check that it resized properly */
  EXPECT_GT( capacity, old_capacity );
  EXPECT_EQ( v.capacity(), capacity );
  EXPECT_EQ( v.size(), size );
  for ( IndexType i = 0 ; i < size ; ++i )
  {
    for ( IndexType j = 0 ; j < num_components ; ++j )
    {
      EXPECT_EQ( v( i, j ), i * j - 5 * i + 7 * j );
    }
  }

  /* Append 1000 tuples */
  const IndexType n_tuples = 1000;
  T values[ n_tuples * num_components ];
  for ( IndexType i = 0 ; i < n_tuples ; ++i )
  {
    for ( IndexType j = 0 ; j < num_components ; ++j )
    {
      IndexType i_real = i + size;
      values[ i * num_components + j ] = i_real * j - 5 * i_real + 7 * j;
    }
  }

  capacity = calc_new_capacity( v, n_tuples );
  v.append( values, n_tuples );
  size += n_tuples;

  /* Check that it resize properly */
  EXPECT_EQ( v.capacity(), capacity );
  EXPECT_EQ( v.size(), size );
  for ( IndexType i = 0 ; i < size ; ++i )
  {
    for ( IndexType j = 0 ; j < num_components ; ++j )
    {
      EXPECT_EQ( v( i, j ), i * j - 5 * i + 7 * j );
    }
  }

  /* Reduce the size down to 500 tuples */
  T* data_address = v.getData();
  size = 500;
  v.resize(size);
  EXPECT_EQ( v.size(), size );
  EXPECT_EQ( v.capacity(), capacity );
  EXPECT_EQ( v.getData(), data_address );
  for ( IndexType i = 0 ; i < size ; ++i )
  {
    for ( IndexType j = 0 ; j < num_components ; ++j )
    {
      EXPECT_EQ( v( i, j ), i * j - 5 * i + 7 * j );
    }
  }

  /* Shrink the vector */
  capacity = size;
  v.shrink();

  EXPECT_EQ( v.capacity(), capacity );
  EXPECT_EQ( v.size(), size );
  for ( IndexType i = 0 ; i < size ; ++i )
  {
    for ( IndexType j = 0 ; j < num_components ; ++j )
    {
      EXPECT_EQ( v( i, j ), i * j - 5 * i + 7 * j );
    }
  }

  /* Append a new tuple, should resize. */
  old_capacity = capacity;
  capacity = calc_new_capacity( v, 1 );
  for ( IndexType j = 0 ; j < num_components ; ++j )
  {
    tuple[ j ] = size * j - 5 * size + 7 * j;
  }
  v.append( tuple, 1 );
  size++;

  /* Check that it resized properly */
  EXPECT_GT( capacity, old_capacity );
  EXPECT_EQ( v.capacity(), capacity );
  EXPECT_EQ( v.size(), size );
  for ( IndexType i = 0 ; i < size ; ++i )
  {
    for ( IndexType j = 0 ; j < num_components ; ++j )
    {
      EXPECT_EQ( v( i, j ), i * j - 5 * i + 7 * j );
    }
  }

  /* Reset the data */
  T* data_ptr = v.getData();
  for ( IndexType i = 0 ; i < size * num_components ; ++i )
  {
    data_ptr[ i ] = i;
  }

  /* Append a bunch of tuples to fill in up to the capacity. Resize should
     not occur. */
  old_capacity = capacity;
  capacity = calc_new_capacity( v, old_capacity - size + 1 );
  for ( IndexType i = size ; i < old_capacity ; ++i )
  {
    for ( IndexType j = 0 ; j < num_components ; ++j )
    {
      tuple[ j ] = i * num_components + j;
    }

    v.append( tuple, 1 );
    size++;
    EXPECT_EQ( v.capacity(), old_capacity );
    EXPECT_EQ( v.size(), size );
    EXPECT_EQ( v.getData(), data_ptr );
  }

  EXPECT_EQ( v.size(), old_capacity );

  /* Append a final tuple that should trigger a resize. */
  for ( IndexType j = 0 ; j < num_components ; ++j )
  {
    tuple[ j ] = size * num_components + j;
  }

  v.append( tuple, 1 );
  size++;
  EXPECT_GT( capacity, old_capacity );
  EXPECT_EQ( v.capacity(), capacity );
  EXPECT_EQ( v.size(), size );

  for ( IndexType i = 0 ; i < size ; ++i )
  {
    for ( IndexType j = 0 ; j < num_components ; ++j )
    {
      EXPECT_EQ( v( i, j ), i * num_components + j );
    }
  }

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
  IndexType num_components = v.numComponents();

  EXPECT_EQ( v.size(), v.capacity() );

  /* Set the existing data in v */
  for ( IndexType i = 0 ; i < size ; ++i )
  {
    for ( IndexType j = 0 ; j < num_components ; ++j )
    {
      v( i, j ) = i * j - 5 * i + 7 * j;
    }
  }

  /* Append a new tuple, should resize. */
  IndexType old_capacity = capacity;
  capacity = calc_new_capacity( v, 1 );
  T tuple[ num_components ];
  for ( IndexType j = 0 ; j < num_components ; ++j )
  {
    tuple[ j ] = size * j - 5 * size + 7 * j;
  }
  v.insert( tuple, 1, v.size() );
  size++;

  /* Check that it resized properly */
  EXPECT_GT( capacity, old_capacity );
  EXPECT_EQ( v.capacity(), capacity );
  EXPECT_EQ( v.size(), size );
  for ( IndexType i = 0 ; i < size ; ++i )
  {
    for ( IndexType j = 0 ; j < num_components ; ++j )
    {
      EXPECT_EQ( v( i, j ), i * j - 5 * i + 7 * j );
    }
  }

  /* Append 1000 tuples */
  const IndexType n_tuples = 1000;
  T values[ n_tuples * num_components ];
  for ( IndexType i = 0 ; i < n_tuples ; ++i )
  {
    for ( IndexType j = 0 ; j < num_components ; ++j )
    {
      IndexType i_real = i + size;
      values[ i * num_components + j ] = i_real * j - 5 * i_real + 7 * j;
    }
  }

  capacity = calc_new_capacity( v, n_tuples );
  v.insert( values, n_tuples, size );
  size += n_tuples;

  /* Check that it resizes properly */
  EXPECT_EQ( v.capacity(), capacity );
  EXPECT_EQ( v.size(), size );
  for ( IndexType i = 0 ; i < size ; ++i )
  {
    for ( IndexType j = 0 ; j < num_components ; ++j )
    {
      EXPECT_EQ( v( i, j ), i * j - 5 * i + 7 * j );
    }
  }

  capacity = size;
  v.shrink();
  IndexType n_insert_front = 100;

  /* Reset the data */
  T* data_ptr = v.getData();
  for ( IndexType i = 0 ; i < size * num_components ; ++i )
  {
    data_ptr[ i ] = i + num_components * n_insert_front;
  }

  /* Insert into the front of the Array. */
  for ( IndexType i = n_insert_front - 1 ; i >= 0 ; i--)
  {
    for ( IndexType j = 0 ; j < num_components ; ++j )
    {
      tuple[ j ] = i * num_components + j;
    }
    capacity = calc_new_capacity( v, 1 );
    v.insert( tuple, 1, 0 );
    size++;
  }

  /* Check that the insertion worked as expected */
  EXPECT_EQ( v.capacity(), capacity );
  EXPECT_EQ( v.size(), size );
  for ( IndexType i = 0 ; i < size ; ++i )
  {
    for ( IndexType j = 0 ; j < num_components ; ++j )
    {
      EXPECT_EQ( v( i, j ), i * num_components + j );
    }
  }

}

/*!
 * \brief Make a copy of an Array through sidre and check it for defects.
 * \param [in] v the Array to copy.
 */
template< typename T >
void check_sidre( Array< T >& v )
{
#ifdef MINT_USE_SIDRE
  Array< T > cpy( const_cast< sidre::View* >( v.getView() ) );
  cpy.setResizeRatio( v.getResizeRatio() );

  check_equality( v, cpy );
  cpy.resize(0);
  check_storage( cpy );
  check_insert( cpy );
#endif
}

/*!
 * \brief Check an external array for defects.
 * \param [in] v the external array to check.
 */
template< typename T >
void check_external( Array< T >& v )
{
  EXPECT_TRUE( v.isExternal() );
  EXPECT_EQ( v.size(), v.capacity() );

  const IndexType size = v.size();
  const IndexType num_components = v.numComponents();
  const IndexType num_values = size * num_components;
  T* const data_ptr = v.getData();

  for ( IndexType i = 0 ; i < size ; ++i )
  {
    for ( IndexType j = 0 ; j < num_components ; ++j )
    {
      v( i, j ) = i * num_components + j;
    }
  }

  for (IndexType i = 0 ; i < num_values ; ++i )
  {
    EXPECT_EQ( data_ptr[ i ], i );
  }

  for ( IndexType i = 0 ; i < size ; ++i )
  {
    for ( IndexType j = 0 ; j < num_components ; ++j )
    {
      data_ptr[ i * num_components + j ] = i * j - i - j;
    }
  }

  for ( IndexType i = 0 ; i < size ; ++i )
  {
    for ( IndexType j = 0 ; j < num_components ; ++j )
    {
      EXPECT_EQ( v(i, j), i * j - i - j );
    }
  }

  EXPECT_EQ( size, v.size() );
  EXPECT_EQ( data_ptr, v.getData() );

  T tuple[ num_components ];
  EXPECT_DEATH_IF_SUPPORTED( v.append( tuple, 1 ), IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( v.insert( tuple, 1, 0 ), IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( v.reserve( size + 1 ), IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( v.shrink(), IGNORE_OUTPUT );
}

}   /* end namespace internal */

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
TEST( mint_core_array, checkStorage )
{
#ifdef MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();
#endif

  constexpr IndexType ZERO = 0;

  for ( IndexType capacity = 2 ; capacity < 1024 ; capacity *= 2 )
  {
    for ( IndexType n_components = 1 ; n_components <= 4 ; n_components++ )
    {
      Array< int > v_int( ZERO, n_components, capacity );
      internal::check_storage( v_int );

      Array< double > v_double( ZERO, n_components, capacity );
      internal::check_storage( v_double );

#ifdef MINT_USE_SIDRE
      Array< int > v_int_sidre( root->createView( "int" ), ZERO, n_components,
                                capacity );
      internal::check_storage( v_int_sidre );

      Array< double > v_double_sidre( root->createView( "double" ), ZERO,
                                      n_components, capacity );
      internal::check_storage( v_double_sidre );

      root->destroyViewsAndData();
#endif
    }
  }
}

//------------------------------------------------------------------------------
TEST( mint_core_array, checkFill )
{
#ifdef MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();
#endif

  for ( IndexType capacity = 2 ; capacity < 1024 ; capacity *= 2 )
  {
    IndexType size = capacity / 2;  
    for ( IndexType n_components = 1 ; n_components <= 4 ; n_components++ )
    {
      Array< int > v_int( size, n_components, capacity );
      internal::check_fill( v_int );

      Array< double > v_double( size, n_components, capacity );
      internal::check_fill( v_double );

#ifdef MINT_USE_SIDRE
      Array< int > v_int_sidre( root->createView( "int" ), size, n_components,
                                capacity );
      internal::check_fill( v_int_sidre );

      Array< double > v_double_sidre( root->createView( "double" ), size,
                                      n_components, capacity );
      internal::check_fill( v_double_sidre );

      root->destroyViewsAndData();
#endif
    }
  }
}

//------------------------------------------------------------------------------
TEST( mint_core_array, checkSet )
{
#ifdef MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();
#endif

  for ( IndexType capacity = 2 ; capacity < 1024 ; capacity *= 2 )
  {
    IndexType size = capacity / 2;  
    for ( IndexType n_components = 1 ; n_components <= 4 ; n_components++ )
    {
      Array< int > v_int( size, n_components, capacity );
      internal::check_set( v_int );

      Array< double > v_double( size, n_components, capacity );
      internal::check_set( v_double );

#ifdef MINT_USE_SIDRE
      Array< int > v_int_sidre( root->createView( "int" ), size, n_components,
                                capacity );
      internal::check_set( v_int_sidre );

      Array< double > v_double_sidre( root->createView( "double" ), size,
                                      n_components, capacity );
      internal::check_set( v_double_sidre );

      root->destroyViewsAndData();
#endif
    }
  }
}

//------------------------------------------------------------------------------
TEST( mint_core_array, checkResize )
{
#ifdef MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();
#endif

  constexpr IndexType ZERO = 0;

  /* Resizing isn't allowed with a ratio less than 1.0. */
  Array< int > v_int( ZERO, 1, 100 );
  v_int.setResizeRatio( 0.99 );
  EXPECT_DEATH_IF_SUPPORTED( internal::check_resize( v_int ), IGNORE_OUTPUT );

  for ( double ratio = 1.0; ratio <= 3.0; ratio += 0.5 )
  {
    for ( IndexType capacity = 2 ; capacity <= 1024 ; capacity *= 2 )
    {
      for ( IndexType n_components = 1 ; n_components <= 4 ; n_components++ )
      {
        Array< int > v_int( ZERO, n_components, capacity );
        v_int.setResizeRatio( ratio );
        internal::check_resize( v_int );

        Array< double > v_double( ZERO, n_components, capacity );
        v_double.setResizeRatio( ratio );
        internal::check_resize( v_double );

#ifdef MINT_USE_SIDRE
        Array< int > v_int_sidre( root->createView( "int" ), ZERO, n_components,
                                  capacity );
        v_int_sidre.setResizeRatio( ratio );
        internal::check_resize( v_int_sidre );

        Array< double > v_double_sidre( root->createView( "double" ), ZERO,
                                        n_components, capacity );
        v_double_sidre.setResizeRatio( ratio );
        internal::check_resize( v_double_sidre );

        root->destroyViewsAndData();
#endif
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST( mint_core_array_DeathTest, checkResize )
{
  /* Resizing isn't allowed with a ratio less than 1.0. */
  Array< int > v_int( 0, 1, 100 );
  v_int.setResizeRatio( 0.99 );
  EXPECT_DEATH_IF_SUPPORTED( internal::check_resize( v_int ), IGNORE_OUTPUT );
}

//------------------------------------------------------------------------------
TEST( mint_core_array, checkInsert )
{
#ifdef MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();
#endif

  constexpr IndexType ZERO = 0;

  for ( double ratio = 1.0 ; ratio <= 3.0 ; ratio += 0.5 )
  {
    for ( IndexType capacity = 2 ; capacity <= 1024 ; capacity *= 2 )
    {
      for ( IndexType n_components = 1 ; n_components <= 4 ; n_components++ )
      {
        Array< int > v_int( ZERO, n_components, capacity );
        v_int.setResizeRatio( ratio );
        internal::check_insert( v_int );

        Array< double > v_double( ZERO, n_components, capacity );
        v_double.setResizeRatio( ratio );
        internal::check_insert( v_double );

#ifdef MINT_USE_SIDRE
        Array< int > v_int_sidre( root->createView("int"), ZERO, n_components,
                                  capacity );
        v_int_sidre.setResizeRatio( ratio );
        internal::check_insert( v_int_sidre );

        Array< double > v_double_sidre( root->createView("double"), ZERO,
                                        n_components, capacity );
        v_double_sidre.setResizeRatio( ratio );
        internal::check_insert( v_double_sidre );

        root->destroyViewsAndData();
#endif
      }
    }
  }
}

/* Sidre specific tests */
#ifdef MINT_USE_SIDRE

//------------------------------------------------------------------------------
TEST( mint_core_array, checkSidre )
{
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  constexpr IndexType ZERO = 0;

  for ( double ratio = 1.0 ; ratio <= 3.0 ; ratio += 0.5 )
  {
    for ( IndexType capacity = 2 ; capacity <= 1024 ; capacity *= 2 )
    {
      for ( IndexType n_components = 1 ; n_components <= 4 ; n_components++ )
      {
        Array< int > v_int( root->createView("int"), ZERO, n_components,
                            capacity );
        v_int.setResizeRatio( ratio );
        internal::check_storage( v_int );
        internal::check_sidre( v_int );

        Array< double > v_double( root->createView("double"), ZERO, n_components,
                                  capacity );
        v_double.setResizeRatio( ratio );
        internal::check_storage( v_double );
        internal::check_sidre( v_double );

        root->destroyViewsAndData();
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST( mint_core_array, checkSidrePermanence)
{
  constexpr double MAGIC_NUM = 5683578.8;

  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  constexpr IndexType ZERO = 0;

  for ( double ratio = 1.0 ; ratio <= 3.0 ; ratio += 0.5 )
  {
    for ( IndexType capacity = 2 ; capacity <= 1024 ; capacity *= 2 )
    {
      for ( IndexType n_components = 1 ; n_components <= 4 ; n_components++ )
      {
        const double * array_data_ptr;
        IndexType num_values;
        { /* Begin scope */
          Array< double > v( root->createView("double"), ZERO, n_components,
                                    capacity );
          array_data_ptr = v.getData();
          num_values = v.size() * v.numComponents();
          v.setResizeRatio( ratio );
          internal::check_storage( v );

          /* Set v's data to MAGIC_NUM */
          for ( IndexType i = 0; i < v.size(); ++i )
          {
            for ( IndexType j = 0; j < v.numComponents(); ++j )
            {
              v( i, j ) = MAGIC_NUM;
            }
          }
        } /* End scope, v has been deallocated. */

        /* Check that the data still exists in sidre */
        sidre::View * view = root->getView( "double" );
        const double * view_data_ptr = static_cast< const double * >(
                                                           view->getVoidPtr() );
        EXPECT_EQ( view_data_ptr, array_data_ptr );
        EXPECT_EQ( view->getNumDimensions(), 2 );

        sidre::SidreLength dims[2];
        view->getShape( 2, dims );
        EXPECT_EQ( dims[0], capacity );
        EXPECT_EQ( dims[1], n_components );

        for ( IndexType i = 0; i < num_values; ++i ) {
          EXPECT_EQ( view_data_ptr[ i ], MAGIC_NUM );
        }

        root->destroyViewsAndData();
      }
    }
  }
}

#endif

//------------------------------------------------------------------------------
TEST( mint_core_array_DeathTest, checkExternal )
{
  constexpr double MAGIC_NUM = 5683578.8;
  constexpr IndexType MAX_SIZE = 1024;
  constexpr IndexType MAX_COMPONENTS = 4;
  constexpr IndexType MAX_VALUES = MAX_SIZE * MAX_COMPONENTS;
  union DataBuffer
  {
    int ints[ MAX_SIZE * MAX_COMPONENTS];
    double doubles[ MAX_SIZE * MAX_COMPONENTS];
  };


  DataBuffer buffer;
  std::fill_n( buffer.doubles, MAX_VALUES, MAGIC_NUM );

  for ( IndexType size = 2 ; size <= MAX_SIZE ; size *= 2 )
  {
    for ( IndexType n_comp = 1 ; n_comp <= MAX_COMPONENTS ; n_comp++ )
    {

      Array< int > v_int( buffer.ints, size, n_comp );
      EXPECT_EQ( v_int.getData(), buffer.ints );
      internal::check_external( v_int );

      Array< double > v_double( buffer.doubles, size, n_comp );
      EXPECT_EQ( v_double.getData(), buffer.doubles );
      internal::check_external( v_double );

      /* Set v_double's data to MAGIC_NUM */
      v_double.fill( MAGIC_NUM );
    }

    /* Check that the data still exists in the buffer */
    for ( IndexType i = 0; i < MAX_VALUES; ++i )
    {
      EXPECT_EQ( buffer.doubles[ i ], MAGIC_NUM );
    }
  }
}

} /* end namespace mint */
} /* end namespace axom */

//------------------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}

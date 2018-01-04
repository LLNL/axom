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

namespace axom
{
namespace mint
{

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

  IndexType num_values = lhs.size() * lhs.numComponents();
  const T* lhs_data = lhs.getData();
  const T* rhs_data = rhs.getData();

  for ( IndexType i = 0 ; i < num_values ; ++i )
  {
    EXPECT_EQ( lhs_data[ i ], rhs_data[ i ] );
  }
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
  Array< T > cpy( const_cast< sidre::View* >( v.getView() ), v.size() );
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

  /* To Do? Check stderr output agains regular expression. */
  T tuple[ num_components ];
  EXPECT_DEATH_IF_SUPPORTED( v.append( tuple, 1 ), ".*" );
  EXPECT_DEATH_IF_SUPPORTED( v.insert( tuple, 1, 0 ), ".*" );
  EXPECT_DEATH_IF_SUPPORTED( v.resize( size ), ".*" );
  EXPECT_DEATH_IF_SUPPORTED( v.reserve( size + 1 ), ".*" );
  EXPECT_DEATH_IF_SUPPORTED( v.shrink(), ".*" );
}

}   /* end namespace internal */

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
TEST( mint_array, checkStorage )
{
#ifdef MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();
#endif

  for ( IndexType capacity = 2 ; capacity < 1024 ; capacity *= 2 )
  {
    for ( int num_components = 1 ; num_components <= 4 ; num_components++ )
    {
      Array< int > v_int( 0, num_components, capacity );
      internal::check_storage( v_int );

      Array< double > v_double( 0, num_components, capacity );
      internal::check_storage( v_double );

#ifdef MINT_USE_SIDRE
      Array< int > v_int_sidre( root->createView( "int" ), 0, num_components,
                                capacity );
      internal::check_storage( v_int_sidre );

      Array< double > v_double_sidre( root->createView( "double" ), 0,
                                      num_components, capacity );
      internal::check_storage( v_double_sidre );

      root->destroyViewsAndData();
#endif
    }
  }
}

//------------------------------------------------------------------------------
TEST( mint_array, checkResize )
{
#ifdef MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();
#endif

  for ( double ratio = 1.0 ; ratio <= 3.0 ; ratio += 0.5 )
  {
    for ( IndexType capacity = 2 ; capacity < 1024 ; capacity *= 2 )
    {
      for ( int num_components = 1 ; num_components <= 4 ; num_components++ )
      {
        Array< int > v_int( 0, num_components, capacity );
        v_int.setResizeRatio( ratio );
        internal::check_resize( v_int );

        Array< double > v_double( 0, num_components, capacity );
        v_double.setResizeRatio( ratio );
        internal::check_resize( v_double );

#ifdef MINT_USE_SIDRE
        Array< int > v_int_sidre( root->createView( "int" ), 0, num_components,
                                  capacity );
        v_int_sidre.setResizeRatio( ratio );
        internal::check_insert( v_int_sidre );

        Array< double > v_double_sidre( root->createView( "double" ), 0,
                                        num_components, capacity );
        v_double_sidre.setResizeRatio( ratio );
        internal::check_insert( v_double_sidre );

        root->destroyViewsAndData();
#endif
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST( mint_array, checkInsert )
{
#ifdef MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();
#endif

  for ( double ratio = 1.0 ; ratio <= 3.0 ; ratio += 0.5 )
  {
    for ( IndexType capacity = 2 ; capacity <= 1024 ; capacity *= 2 )
    {
      for ( int num_components = 1 ; num_components <= 4 ; num_components++ )
      {
        Array< int > v_int( 0, num_components, capacity );
        v_int.setResizeRatio( ratio );
        internal::check_insert( v_int );

        Array< double > v_double( 0, num_components, capacity );
        v_double.setResizeRatio( ratio );
        internal::check_insert( v_double );

#ifdef MINT_USE_SIDRE
        Array< int > v_int_sidre( root->createView("int"), 0, num_components,
                                  capacity );
        v_int_sidre.setResizeRatio( ratio );
        internal::check_insert( v_int_sidre );

        Array< double > v_double_sidre( root->createView("double"), 0,
                                        num_components, capacity );
        v_double_sidre.setResizeRatio( ratio );
        internal::check_insert( v_double_sidre );

        root->destroyViewsAndData();
#endif
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST( mint_array, checkSidre )
{
#ifdef MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  for ( double ratio = 1.0 ; ratio <= 3.0 ; ratio += 0.5 )
  {
    for ( IndexType capacity = 2 ; capacity <= 1024 ; capacity *= 2 )
    {
      for ( int num_components = 1 ; num_components <= 4 ; num_components++ )
      {
        Array< int > v_int( root->createView("int"), 0, num_components,
                            capacity );
        v_int.setResizeRatio( ratio );
        internal::check_storage( v_int );
        internal::check_sidre( v_int );

        Array< double > v_double( root->createView("double"), 0, num_components,
                                  capacity );
        v_double.setResizeRatio( ratio );
        internal::check_storage( v_double );
        internal::check_sidre( v_double );

        root->destroyViewsAndData();
      }
    }
  }

#else
  EXPECT_TRUE( true );
#endif
}

//------------------------------------------------------------------------------
TEST( mint_array_DeathTest, checkExternal )
{
  const int MAX_SIZE = 1024;
  const int MAX_COMPONENTS = 4;
  union DataBuffer
  {
    int ints[ MAX_SIZE * MAX_COMPONENTS];
    double doubles[ MAX_SIZE * MAX_COMPONENTS];
  };

  DataBuffer buffer;
  for ( IndexType size = 2 ; size <= MAX_SIZE ; size *= 2 )
  {
    for ( int n_comp = 1 ; n_comp <= MAX_COMPONENTS ; n_comp++ )
    {
      Array< int > v_int( buffer.ints, size, n_comp );
      internal::check_external( v_int );

      Array< double > v_double( buffer.doubles, size, n_comp );
      internal::check_external( v_double );
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

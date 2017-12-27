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

#include "mint/Array.hpp"

#include "axom_utils/Utilities.hpp"     /* for utilities::max */
#include "mint/DataTypes.hpp"           /* for localIndex */
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
localIndex calc_new_capacity( Array< T > & v, localIndex increase )
{
  localIndex new_num_tuples = v.size() + increase;
  if ( new_num_tuples > v.getCapacity() )
  {
    return new_num_tuples * v.getResizeRatio() + 0.5;
  }

  return v.getCapacity();
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
  EXPECT_EQ( lhs.getNumComponents(), rhs.getNumComponents() );
  EXPECT_EQ( lhs.getCapacity(), rhs.getCapacity() );

  localIndex num_values = lhs.size() * lhs.getNumComponents();
  const T* lhs_data = lhs.getData();
  const T* rhs_data = rhs.getData();

  for ( localIndex i = 0 ; i < num_values ; ++i )
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

  localIndex capacity = v.getCapacity();
  localIndex num_components = v.getNumComponents();
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
    for ( localIndex i = 0 ; i < capacity / 2 ; ++i )
    {
      for ( localIndex j = 0 ; j < num_components ; ++j )
      {
        tuple[ j ] = i * num_components + j;
      }
      v.append( tuple, 1 );
    }
  }

  EXPECT_TRUE( !v.empty() );
  EXPECT_EQ( v.size(), capacity / 2 );
  EXPECT_EQ( v.getCapacity(), capacity );
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
    for ( localIndex i = capacity / 2 ; i < capacity ; ++i )
    {
      for ( localIndex j = 0 ; j < num_components ; ++j )
      {
        tuple[ j ] = i * num_components + j;
      }
      v.append( tuple, 1 );
    }
  }

  EXPECT_TRUE( !v.empty() );
  EXPECT_EQ( v.size(), capacity );
  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.getData(), data_ptr );

  for ( localIndex i = 0 ; i < capacity ; ++i )
  {
    for ( localIndex j = 0 ; j < num_components ; ++j )
    {
      EXPECT_EQ( v( i, j ), i * num_components + j );
      EXPECT_EQ( data_ptr[ i * num_components + j ], i * num_components + j );
    }
  }

  for ( localIndex i = 0 ; i < capacity ; ++i )
  {
    for ( localIndex j = 0 ; j < num_components ; ++j )
    {
      v( i, j ) = i * j - 5 * i + 7 * j;
    }
  }

  for ( localIndex i = 0 ; i < capacity ; ++i )
  {
    for ( localIndex j = 0 ; j < num_components ; ++j )
    {
      EXPECT_EQ( v( i, j ), i * j - 5 * i + 7 * j );
      EXPECT_EQ( data_ptr[ i * num_components + j ], i * j - 5 * i + 7 * j );
    }
  }

  EXPECT_TRUE( !v.empty() );
  EXPECT_EQ( v.size(), capacity );
  EXPECT_EQ( v.getCapacity(), capacity );
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
  localIndex capacity = v.getCapacity();
  v.resize( capacity );
  localIndex size = capacity;
  localIndex num_components = v.getNumComponents();

  EXPECT_EQ( v.size(), v.getCapacity() );

  /* Set the existing data in v */
  for ( localIndex i = 0 ; i < size ; ++i )
  {
    for ( localIndex j = 0 ; j < num_components ; ++j )
    {
      v( i, j ) = i * j - 5 * i + 7 * j;
    }
  }

  /* Append a new tuple, should resize. */
  localIndex old_capacity = capacity;
  capacity = calc_new_capacity( v, 1 );
  T tuple[ num_components ];
  for ( localIndex j = 0 ; j < num_components ; ++j )
  {
    tuple[ j ] = size * j - 5 * size + 7 * j;
  }
  v.append( tuple, 1 );
  size++;

  /* Check that it resized properly */
  EXPECT_GT( capacity, old_capacity );
  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.size(), size );
  for ( localIndex i = 0 ; i < size ; ++i )
  {
    for ( localIndex j = 0 ; j < num_components ; ++j )
    {
      EXPECT_EQ( v( i, j ), i * j - 5 * i + 7 * j );
    }
  }

  /* Append 1000 tuples */
  const localIndex n_tuples = 1000;
  T values[ n_tuples * num_components ];
  for ( localIndex i = 0 ; i < n_tuples ; ++i )
  {
    for ( localIndex j = 0 ; j < num_components ; ++j )
    {
      localIndex i_real = i + size;
      values[ i * num_components + j ] = i_real * j - 5 * i_real + 7 * j;
    }
  }

  capacity = calc_new_capacity( v, n_tuples );
  v.append( values, n_tuples );
  size += n_tuples;

  /* Check that it resize properly */
  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.size(), size );
  for ( localIndex i = 0 ; i < size ; ++i )
  {
    for ( localIndex j = 0 ; j < num_components ; ++j )
    {
      EXPECT_EQ( v( i, j ), i * j - 5 * i + 7 * j );
    }
  }

  /* Reduce the size down to 500 tuples */
  T* data_address = v.getData();
  size = 500;
  v.resize(size);
  EXPECT_EQ( v.size(), size );
  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.getData(), data_address );
  for ( localIndex i = 0 ; i < size ; ++i )
  {
    for ( localIndex j = 0 ; j < num_components ; ++j )
    {
      EXPECT_EQ( v( i, j ), i * j - 5 * i + 7 * j );
    }
  }

  /* Shrink the vector */
  capacity = size;
  v.shrink();

  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.size(), size );
  for ( localIndex i = 0 ; i < size ; ++i )
  {
    for ( localIndex j = 0 ; j < num_components ; ++j )
    {
      EXPECT_EQ( v( i, j ), i * j - 5 * i + 7 * j );
    }
  }

  /* Append a new tuple, should resize. */
  old_capacity = capacity;
  capacity = calc_new_capacity( v, 1 );
  for ( localIndex j = 0 ; j < num_components ; ++j )
  {
    tuple[ j ] = size * j - 5 * size + 7 * j;
  }
  v.append( tuple, 1 );
  size++;

  /* Check that it resized properly */
  EXPECT_GT( capacity, old_capacity );
  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.size(), size );
  for ( localIndex i = 0 ; i < size ; ++i )
  {
    for ( localIndex j = 0 ; j < num_components ; ++j )
    {
      EXPECT_EQ( v( i, j ), i * j - 5 * i + 7 * j );
    }
  }

  /* Reset the data */
  T* data_ptr = v.getData();
  for ( localIndex i = 0 ; i < size * num_components ; ++i )
  {
    data_ptr[ i ] = i;
  }

  /* Append a bunch of tuples to fill in up to the capacity. Resize should
     not occur. */
  old_capacity = capacity;
  capacity = calc_new_capacity( v, old_capacity - size + 1 );
  for ( localIndex i = size ; i < old_capacity ; ++i )
  {
    for ( localIndex j = 0 ; j < num_components ; ++j )
    {
      tuple[ j ] = i * num_components + j;
    }

    v.append( tuple, 1 );
    size++;
    EXPECT_EQ( v.getCapacity(), old_capacity );
    EXPECT_EQ( v.size(), size );
    EXPECT_EQ( v.getData(), data_ptr );
  }

  EXPECT_EQ( v.size(), old_capacity );

  /* Append a final tuple that should trigger a resize. */
  for ( localIndex j = 0 ; j < num_components ; ++j )
  {
    tuple[ j ] = size * num_components + j;
  }

  v.append( tuple, 1 );
  size++;
  EXPECT_GT( capacity, old_capacity );
  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.size(), size );

  for ( localIndex i = 0 ; i < size ; ++i )
  {
    for ( localIndex j = 0 ; j < num_components ; ++j )
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
  localIndex capacity = v.getCapacity();
  v.resize( capacity );
  localIndex size = capacity;
  localIndex num_components = v.getNumComponents();

  EXPECT_EQ( v.size(), v.getCapacity() );

  /* Set the existing data in v */
  for ( localIndex i = 0 ; i < size ; ++i )
  {
    for ( localIndex j = 0 ; j < num_components ; ++j )
    {
      v( i, j ) = i * j - 5 * i + 7 * j;
    }
  }

  /* Append a new tuple, should resize. */
  localIndex old_capacity = capacity;
  capacity = calc_new_capacity( v, 1 );
  T tuple[ num_components ];
  for ( localIndex j = 0 ; j < num_components ; ++j )
  {
    tuple[ j ] = size * j - 5 * size + 7 * j;
  }
  v.insert( tuple, 1, v.size() );
  size++;

  /* Check that it resized properly */
  EXPECT_GT( capacity, old_capacity );
  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.size(), size );
  for ( localIndex i = 0 ; i < size ; ++i )
  {
    for ( localIndex j = 0 ; j < num_components ; ++j )
    {
      EXPECT_EQ( v( i, j ), i * j - 5 * i + 7 * j );
    }
  }

  /* Append 1000 tuples */
  const localIndex n_tuples = 1000;
  T values[ n_tuples * num_components ];
  for ( localIndex i = 0 ; i < n_tuples ; ++i )
  {
    for ( localIndex j = 0 ; j < num_components ; ++j )
    {
      localIndex i_real = i + size;
      values[ i * num_components + j ] = i_real * j - 5 * i_real + 7 * j;
    }
  }

  capacity = calc_new_capacity( v, n_tuples );
  v.insert( values, n_tuples, size );
  size += n_tuples;

  /* Check that it resizes properly */
  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.size(), size );
  for ( localIndex i = 0 ; i < size ; ++i )
  {
    for ( localIndex j = 0 ; j < num_components ; ++j )
    {
      EXPECT_EQ( v( i, j ), i * j - 5 * i + 7 * j );
    }
  }

  capacity = size;
  v.shrink();
  localIndex n_insert_front = 100;

  /* Reset the data */
  T* data_ptr = v.getData();
  for ( localIndex i = 0 ; i < size * num_components ; ++i )
  {
    data_ptr[ i ] = i + num_components * n_insert_front;
  }

  /* Insert into the front of the Array. */
  for ( localIndex i = n_insert_front - 1 ; i >= 0 ; i--)
  {
    for ( localIndex j = 0 ; j < num_components ; ++j )
    {
      tuple[ j ] = i * num_components + j;
    }
    capacity = calc_new_capacity( v, 1 );
    v.insert( tuple, 1, 0 );
    size++;
  }

  /* Check that the insertion worked as expected */
  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.size(), size );
  for ( localIndex i = 0 ; i < size ; ++i )
  {
    for ( localIndex j = 0 ; j < num_components ; ++j )
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
void load_and_check( Array< T >& v )
{
#ifdef MINT_USE_SIDRE
  Array< T > cpy = Array< T >(v.getView(), v.size() );
  cpy.setResizeRatio( v.getResizeRatio() );

  check_equality( v, cpy );
  cpy.resize(0);
  check_storage( cpy );
  check_resize( cpy );
#endif
}

}   /* end namespace internal */


//------------------------------------------------------------------------------
TEST( mint_array, checkStorage )
{
#ifdef MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();
#endif

  localIndex capacity = 2;
  for ( int i = 1 ; i < 10 ; ++i )
  {
    for ( int num_components = 1 ; num_components < 5 ; num_components++ )
    {
      Array< int > v_int = Array< int >( capacity, 0, num_components );
      internal::check_storage( v_int );

      Array< long > v_long = Array< long >( capacity, 0, num_components );
      internal::check_storage( v_long );

      Array< float > v_float = Array< float >( capacity, 0, num_components );
      internal::check_storage( v_float );

      Array< double > v_double = Array< double >( capacity, 0, num_components );
      internal::check_storage( v_double );

#ifdef MINT_USE_SIDRE
      v_int = Array< int >( root->createView("int"), capacity, 0,
                            num_components);
      internal::check_storage( v_int );
      internal::load_and_check( v_int );

      v_long  = Array< long >( root->createView("long"), capacity, 0,
                               num_components);
      internal::check_storage( v_long );
      internal::load_and_check( v_long );

      v_float  = Array< float >( root->createView("float"), capacity, 0,
                                 num_components);
      internal::check_storage( v_float );
      internal::load_and_check( v_float );

      v_double  = Array< double >( root->createView("double"), capacity, 0,
                                   num_components);
      internal::check_storage( v_double );
      internal::load_and_check( v_double );


      root->destroyViewsAndData();
#endif
    }
    capacity *= 2;
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
    localIndex capacity = 2;
    for ( int i = 1 ; i < 10 ; ++i )
    {
      for ( int num_components = 1 ; num_components < 5 ; num_components++ )
      {
        Array< int > v_int = Array< int >( capacity, 0, num_components );
        v_int.setResizeRatio( ratio );
        internal::check_resize( v_int );

        Array< long > v_long = Array< long >( capacity, 0, num_components );
        v_long.setResizeRatio( ratio );
        internal::check_resize( v_long );

        Array< float > v_float = Array< float >( capacity, 0, num_components );
        v_float.setResizeRatio( ratio );
        internal::check_resize( v_float );

        Array< double > v_double =
          Array< double >( capacity, 0, num_components );
        v_double.setResizeRatio( ratio );
        internal::check_resize( v_double );

#ifdef MINT_USE_SIDRE
        v_int = Array< int >( root->createView("int"), capacity, 0,
                              num_components);
        v_int.setResizeRatio( ratio );
        internal::check_resize( v_int );
        internal::load_and_check( v_int );

        v_long  = Array< long >( root->createView("long"), capacity, 0,
                                 num_components);
        v_long.setResizeRatio( ratio );
        internal::check_resize( v_long );
        internal::load_and_check( v_long );

        v_float  = Array< float >( root->createView("float"), capacity, 0,
                                   num_components);
        v_float.setResizeRatio( ratio );
        internal::check_resize( v_float );
        internal::load_and_check( v_float );

        v_double  = Array< double >( root->createView("double"), capacity, 0,
                                     num_components);
        v_double.setResizeRatio( ratio );
        internal::check_resize( v_double );
        internal::load_and_check( v_double );

        root->destroyViewsAndData();
#endif
      }
      capacity *= 2;
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
    localIndex capacity = 2;
    for ( int i = 1 ; i < 10 ; ++i )
    {
      for ( int num_components = 1 ; num_components < 5 ; num_components++ )
      {
        Array< int > v_int = Array< int >( capacity, 0, num_components );
        v_int.setResizeRatio( ratio );
        internal::check_insert( v_int );

        Array< long > v_long = Array< long >( capacity, 0, num_components );
        v_long.setResizeRatio( ratio );
        internal::check_insert( v_long );

        Array< float > v_float = Array< float >( capacity, 0, num_components );
        v_float.setResizeRatio( ratio );
        internal::check_insert( v_float );

        Array< double > v_double =
          Array< double >( capacity, 0, num_components );
        v_double.setResizeRatio( ratio );
        internal::check_insert( v_double );

#ifdef MINT_USE_SIDRE
        v_int = Array< int >( root->createView("int"), capacity, 0,
                              num_components);
        v_int.setResizeRatio( ratio );
        internal::check_insert( v_int );
        internal::load_and_check( v_int );

        v_long  = Array< long >( root->createView("long"), capacity, 0,
                                 num_components);
        v_long.setResizeRatio( ratio );
        internal::check_insert( v_long );
        internal::load_and_check( v_long );

        v_float  = Array< float >( root->createView("float"), capacity, 0,
                                   num_components);
        v_float.setResizeRatio( ratio );
        internal::check_insert( v_float );
        internal::load_and_check( v_float );

        v_double  = Array< double >( root->createView("double"), capacity, 0,
                                     num_components);
        v_double.setResizeRatio( ratio );
        internal::check_insert( v_double );
        internal::load_and_check( v_double );

        root->destroyViewsAndData();
#endif
      }
      capacity *= 2;
    }
  }
}


} /* end namespace mint */
} /* end namespace axom */

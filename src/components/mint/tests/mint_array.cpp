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

//------------------------------------------------------------------------------
template< typename T >
void check_equality( const Array< T >& lhs, const Array< T >& rhs )
{
  EXPECT_EQ( lhs.size(), rhs.size() );
  EXPECT_EQ( lhs.getNumComponents(), rhs.getNumComponents() );
  EXPECT_EQ( lhs.getCapacity(), rhs.getCapacity() );

  localIndex num_values = lhs.size() * lhs.getNumComponents();
  const T* lhs_data = lhs.getData();
  const T* rhs_data = rhs.getData();

  for ( localIndex i = 0; i < num_values; ++i ) {
    EXPECT_EQ( lhs_data[ i ], rhs_data[ i ] );
  }
}

//------------------------------------------------------------------------------
template< typename T >
void load_and_check( Array< T >& v )
{
#ifdef MINT_USE_SIDRE
  Array< T > cpy = Array< T >(v.getView(), v.size() );
  check_equality( v, cpy );
#endif
}

//------------------------------------------------------------------------------
template < typename T >
localIndex calc_new_capacity( Array< T > & v, localIndex increase )
{
  localIndex new_num_tuples = v.size() + increase;
  if ( new_num_tuples > v.getCapacity() )
  { return new_num_tuples * v.getResizeRatio() + 0.5; }

  return v.getCapacity();
}

//------------------------------------------------------------------------------
template < typename T >
void check_storage( Array< T >& v )
{
  EXPECT_TRUE( v.empty() );
  EXPECT_EQ( v.size(), 0 );

  localIndex capacity = v.getCapacity();
  localIndex num_components = v.getNumComponents();
  const T* data_ptr = v.getData();

  if ( num_components == 1 ) {
    for ( T i = 0 ; i < capacity / 2 ; ++i ) { 
      v.append( i );
    }
  }
  else {
    T tuple[ num_components ];
    for ( localIndex i = 0 ; i < capacity / 2 ; ++i ) {
      for ( localIndex j = 0; j < num_components; ++j ) {
        tuple[ j ] = i * num_components + j;
      }
      v.append( tuple, 1 );
    }
  }

  EXPECT_TRUE( !v.empty() );
  EXPECT_EQ( v.size(), capacity / 2 );
  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.getData(), data_ptr );

  if ( num_components == 1 ) {
    for ( T i = capacity / 2; i < capacity; ++i ) { 
      v.append( i );
    }
  }
  else {
    T tuple[ num_components ];
    for ( localIndex i = capacity / 2; i < capacity; ++i ) {
      for ( localIndex j = 0; j < num_components; ++j ) {
        tuple[ j ] = i * num_components + j;
      }
      v.append( tuple, 1 );
    }
  }

  EXPECT_TRUE( !v.empty() );
  EXPECT_EQ( v.size(), capacity );
  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.getData(), data_ptr );

  for ( localIndex i = 0 ; i < capacity ; ++i ) {
    for ( localIndex j = 0; j < num_components; ++j ) {
      EXPECT_EQ( v( i, j ), i * num_components + j );
      EXPECT_EQ( data_ptr[ i * num_components + j ], i * num_components + j );
    }
  }

  for ( localIndex i = 0 ; i < capacity ; ++i ) {
    for ( localIndex j = 0; j < num_components; ++j ) {
      v( i, j ) = i * j - 5 * i + 7 * j;
    }
  }

  for ( localIndex i = 0 ; i < capacity ; ++i ) {
    for ( localIndex j = 0; j < num_components; ++j ) {
      EXPECT_EQ( v( i, j ), i * j - 5 * i + 7 * j );
      EXPECT_EQ( data_ptr[ i * num_components + j ], i * j - 5 * i + 7 * j );
    }
  }

  EXPECT_TRUE( !v.empty() );
  EXPECT_EQ( v.size(), capacity );
  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.getData(), data_ptr );
}

#if 0

//------------------------------------------------------------------------------
template < typename T >
void check_resize( int num_components, double resize_ratio )
{
  localIndex capacity = 0;
  localIndex size = 0;
  Array< T > v( capacity, num_components );

  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.size(), size );
  EXPECT_EQ( v.getNumComponents(), num_components );
  EXPECT_EQ( v.getResizeRatio(), resize_ratio );

  capacity = calc_new_capacity( v, 1 );
  v.add( 0 );
  size += 1;

  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.size(), size );
  EXPECT_EQ( v[0], 0 );

  const localIndex n_vals = 1000 * num_components;
  T values[ n_vals ];
  for ( localIndex i = 0 ; i < n_vals ; ++i )
  {
    values[ i ] = i + 1;
  }

  capacity = calc_new_capacity( v, n_vals );
  v.add( values, n_vals );
  size += n_vals;

  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.size(), size );
  for ( localIndex i = 0 ; i < size ; ++i )
  {
    EXPECT_EQ( v[ i ], i );
  }

  capacity = calc_new_capacity( v, 1 );
  v.add( n_vals + 1 );
  size += 1;

  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.size(), size );
  EXPECT_EQ( v[ size - 1 ], n_vals + 1 );

  T* data_address = v.getData();
  size = 500 * num_components;
// TODO: ???
//  v.setSize(size);
  EXPECT_EQ( v.size(), size );
  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.getData(), data_address );
  for ( localIndex i = 0 ; i < size ; ++i )
  {
    EXPECT_EQ( v[ i ], i );
  }

  capacity = 250 * num_components;
  size = capacity;

// TODO: ???
//  v.setCapacity( capacity );

  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.size(), size );
  for ( localIndex i = 0 ; i < size ; ++i )
  {
    EXPECT_EQ( v[ i ], i );
  }

  size += n_vals;
  capacity = std::ceil( size / num_components ) * num_components;

// TODO: ???
//  v.setSize( size );

  EXPECT_EQ( v.size(), size );
  EXPECT_EQ( v.getCapacity(), capacity );

  v.set( values, n_vals, size - n_vals );
  EXPECT_EQ( v.size(), size );
  EXPECT_EQ( v.getCapacity(), capacity );

  for ( localIndex i = 0 ; i < size - n_vals ; ++i )
  {
    EXPECT_EQ( v[ i ], i );
  }

  for ( localIndex i = size - n_vals ; i < size ; ++i )
  {
    EXPECT_EQ( v[ i ], values[i - size + n_vals] );
  }

  capacity = 1000 * num_components;
  size = capacity;
// TODO: ???
//  v.setCapacity( capacity );
//  v.setSize( size );
  data_address = v.getData();

  for ( localIndex i = 0 ; i < size ; ++i )
  {
    data_address[ i ] = i * i;
  }

  EXPECT_EQ( v.size(), size );
  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.getData(), data_address );

  for ( localIndex i = 0 ; i < size ; ++i )
  {
    EXPECT_EQ( v[ i ], i * i );
  }
}

//------------------------------------------------------------------------------
template < typename T >
void check_insert( double resize_ratio )
{
  localIndex capacity = 0;
  localIndex size = 0;
  Array< T > v( capacity, resize_ratio );

  localIndex num_vals = 1000;

  for ( localIndex i = 0 ; i < num_vals ; ++i )
  {
    capacity = calc_new_capacity( v, 1 );
    v.insert( 2 * i, size );
    size++;
  }

  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.size(), size );
  EXPECT_EQ( v.size(), num_vals );

  for ( localIndex i = 0 ; i < num_vals ; ++i )
  {
    EXPECT_EQ( v[ i ], 2 * i );
  }

  for ( localIndex i = 0 ; i < num_vals ; ++i )
  {
    capacity = calc_new_capacity( v, 1 );
    v.insert( 2 * i + 1, 2 * i + 1 );
    size++;
  }

  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.size(), size );
  EXPECT_EQ( v.size(), 2 * num_vals );

  for ( localIndex i = 0 ; i < 2 * num_vals ; ++i )
  {
    EXPECT_EQ( v[ i ], i );
  }

  T values[ num_vals ];
  for ( localIndex i = 0 ; i < num_vals ; ++i )
  {
    values[ i ] = i * ( i - 4 ) + 5;
  }

  capacity = calc_new_capacity( v, num_vals );
  size += num_vals;
  v.insert( values, num_vals, 0 );

  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.size(), size );
  EXPECT_EQ( v.size(), 3 * num_vals );

  for ( localIndex i = 0 ; i < num_vals ; ++i )
  {
    EXPECT_EQ( v[ i ], values[ i ] );
  }

  for ( localIndex i = num_vals ; i < 3 * num_vals ; ++i )
  {
    EXPECT_EQ( v[ i ], i - num_vals );
  }
}

#endif

}   /* end namespace internal */


//------------------------------------------------------------------------------
TEST( mint_array, checkStorage )
{
#ifdef MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();
#endif  

  localIndex capacity = 2;
  for ( int i = 1 ; i < 10 ; ++i ) {
    for ( int num_components = 1; num_components < 5; num_components++ ) {
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

#if 0
//------------------------------------------------------------------------------
TEST( mint_array, checkResize )
{
  for ( int num_components = 1 ; num_components < 5 ; ++num_components )
  {
    for ( double resize_ratio = 1.0 ; resize_ratio < 3 ; resize_ratio += 0.3 )
    {
      internal::check_resize< int > ( num_components, resize_ratio );
      internal::check_resize< long int >( num_components, resize_ratio );
      internal::check_resize< float >( num_components, resize_ratio );
      internal::check_resize< double >( num_components, resize_ratio );
    }
  }
}

//------------------------------------------------------------------------------
TEST( mint_array, checkInsert )
{
  for ( double resize_ratio = 1.0 ; resize_ratio < 3 ; resize_ratio += 0.3 )
  {
    internal::check_insert< int > ( resize_ratio );
    internal::check_insert< long int >( resize_ratio );
    internal::check_insert< float >( resize_ratio );
    internal::check_insert< double >( resize_ratio );
  }
}

/* test copy / move */

#endif

} /* end namespace mint */
} /* end namespace axom */


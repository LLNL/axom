/*
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

#include "mint/Vector.hpp"
#include "mint/DataTypes.hpp"           /* for localIndex */
#include "axom_utils/Utilities.hpp"     /* for utilities::max */
#include "gtest/gtest.h"                /* for TEST and EXPECT_* macros */
#include "slic/slic.hpp"                /* for slic macros */
#include "slic/UnitTestLogger.hpp"      /* for UnitTestLogger */

#include <cmath>                        /* for std::ceil */


namespace axom
{
namespace mint
{

namespace internal
{

//------------------------------------------------------------------------------
template < typename T >
localIndex calc_new_capacity( Vector< T > & v, localIndex increase )
{
  localIndex newSize = v.getNumTuples() + increase;
  if ( newSize > v.getCapacity() )
  {
    double n_tuples = newSize * v.getResizeRatio() / v.getNumComponents();
    return std::ceil( n_tuples ) * v.getNumComponents();
  }

  return v.getCapacity();
}

//------------------------------------------------------------------------------
template < typename T >
void check_storage( localIndex capacity )
{
  Vector< T > v( capacity );

  EXPECT_TRUE( v.empty() );
  EXPECT_EQ( v.getNumTuples(), 0 );
  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.getResizeRatio(), 2.0 );
  EXPECT_EQ( v.getNumComponents(), 1 );

  const T * data_address = v.getData();

  for ( localIndex i = 0 ; i < capacity / 2 ; ++i )
  {
    v.add( i );
  }

  EXPECT_TRUE( !v.empty() );
  EXPECT_EQ( v.getNumTuples(), capacity / 2 );
  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.getData(), data_address );

  for ( localIndex i =  capacity / 2 ; i < capacity ; ++i )
  {
    v.add( i );
  }

  EXPECT_TRUE( !v.empty() );
  EXPECT_EQ( v.getNumTuples(), capacity );
  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.getData(), data_address );

  for ( localIndex i = 0 ; i < capacity ; ++i )
  {
    EXPECT_EQ( v[ i ], i );
  }

  for ( localIndex i = 0 ; i < capacity ; ++i )
  {
    v[ i ] = capacity - i;
  }

  for ( localIndex i = 0 ; i < capacity ; ++i )
  {
    EXPECT_EQ( v[ i ], capacity - i );
  }
}


//------------------------------------------------------------------------------
template < typename T >
void check_resize( int num_components, double resize_ratio )
{
  localIndex capacity = 0;
  localIndex size = 0;
  Vector< T > v( capacity, num_components );

  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.getNumTuples(), size );
  EXPECT_EQ( v.getNumComponents(), num_components );
  EXPECT_EQ( v.getResizeRatio(), resize_ratio );

  capacity = calc_new_capacity( v, 1 );
  v.add( 0 );
  size += 1;

  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.getNumTuples(), size );
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
  EXPECT_EQ( v.getNumTuples(), size );
  for ( localIndex i = 0 ; i < size ; ++i )
  {
    EXPECT_EQ( v[ i ], i );
  }

  capacity = calc_new_capacity( v, 1 );
  v.add( n_vals + 1 );
  size += 1;

  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.getNumTuples(), size );
  EXPECT_EQ( v[ size - 1 ], n_vals + 1 );

  T * data_address = v.getData();
  size = 500 * num_components;
// TODO: ???
//  v.setSize(size);
  EXPECT_EQ( v.getNumTuples(), size );
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
  EXPECT_EQ( v.getNumTuples(), size );
  for ( localIndex i = 0 ; i < size ; ++i )
  {
    EXPECT_EQ( v[ i ], i );
  }

  size += n_vals;
  capacity = std::ceil( size / num_components ) * num_components;

// TODO: ???
//  v.setSize( size );

  EXPECT_EQ( v.getNumTuples(), size );
  EXPECT_EQ( v.getCapacity(), capacity );

  v.set( values, n_vals, size - n_vals );
  EXPECT_EQ( v.getNumTuples(), size );
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

  EXPECT_EQ( v.getNumTuples(), size );
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
  Vector< T > v( capacity, resize_ratio );

  localIndex num_vals = 1000;

  for ( localIndex i = 0 ; i < num_vals ; ++i )
  {
    capacity = calc_new_capacity( v, 1 );
    v.insert( 2 * i, size );
    size++;
  }

  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.getNumTuples(), size );
  EXPECT_EQ( v.getNumTuples(), num_vals );

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
  EXPECT_EQ( v.getNumTuples(), size );
  EXPECT_EQ( v.getNumTuples(), 2 * num_vals );

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
  EXPECT_EQ( v.getNumTuples(), size );
  EXPECT_EQ( v.getNumTuples(), 3 * num_vals );

  for ( localIndex i = 0 ; i < num_vals ; ++i )
  {
    EXPECT_EQ( v[ i ], values[ i ] );
  }

  for ( localIndex i = num_vals ; i < 3 * num_vals ; ++i )
  {
    EXPECT_EQ( v[ i ], i - num_vals );
  }
}

}   /* end namespace internal */

//------------------------------------------------------------------------------
TEST( mint_vector, checkStorage )
{
  localIndex capacity = 2;
  for ( int i = 1 ; i < 17 ; ++i )
  {
    internal::check_storage< int > ( capacity );
    internal::check_storage< long int >( capacity );
    internal::check_storage< float >( capacity );
    internal::check_storage< double >( capacity );
    capacity *= 2;
  }
}

//------------------------------------------------------------------------------
TEST( mint_vector, checkResize )
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
TEST( mint_vector, checkInsert )
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

} /* end namespace mint */
} /* end namespace axom */

//------------------------------------------------------------------------------
using axom::slic::UnitTestLogger;

int main( int argc, char * argv[] )
{
  int result = 0;
  ::testing::InitGoogleTest( &argc, argv );
  UnitTestLogger logger;
  result = RUN_ALL_TESTS();
  return result;
}

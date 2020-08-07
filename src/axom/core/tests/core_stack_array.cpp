// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/StackArray.hpp"           /* for axom::StackArray */
#include "gtest/gtest.h"                      /* for TEST and EXPECT_* macros */
#include <string>

namespace axom
{

namespace internal
{


template< typename T, int N, typename LAMBDA >
void check( const StackArray< T, N > & arr, LAMBDA && getValue )
{
  for ( int i = 0 ; i < N ; ++i )
  {
    EXPECT_EQ( arr[ i ], getValue( i ) );
  }

  const T* ptr = arr;
  EXPECT_EQ( ptr, &arr[0] );

  for ( int i = 0 ; i < N ; ++i )
  {
    EXPECT_EQ( ptr[ i ], getValue( i ) );
  }
}


template< typename T, int N, typename LAMBDA >
StackArray< T, N > test_storage( LAMBDA && getValue )
{
  StackArray< T, N > arr;

  for ( int i = 0 ; i < N ; ++i )
  {
    arr[ i ] = getValue( i );
  }

  check( arr, std::forward< LAMBDA >( getValue ) );
  return arr;
}


template< typename T, int N, typename LAMBDA >
void test_copy( LAMBDA && getValue )
{
  StackArray< T, N > arr;

  for ( int i = 0 ; i < N ; ++i )
  {
    arr[ i ] = getValue( i );
  }

  StackArray< T, N > copy(arr);
  EXPECT_NE( &arr[0], &copy[0] );

  for ( int i = 0 ; i < N ; ++i )
  {
    EXPECT_EQ( arr[ i ], copy[ i ] );
  }

  check( copy, std::forward< LAMBDA >( getValue ) );
}


template< typename T, typename LAMBDA >
void test_list_initilization( LAMBDA && getValue )
{
  const StackArray< T, 5 > arr = {
    { getValue(0), getValue(1), getValue(2),
      getValue(3), getValue(4) }
  };

  check( arr, std::forward< LAMBDA >( getValue ) );
}

} /* namespace internal */


struct Tensor
{
  double x, y, z;

  Tensor() = default;
  Tensor(const Tensor& other) = default;


  explicit Tensor( double val ) :
    x( 3 * val ), y( 3 * val + 1 ), z( 3 * val + 2 )
  {}

  Tensor& operator=( const Tensor& other )
  {
    x = other.x;
    y = other.y;
    z = other.z;
    return *this;
  }

  bool operator==( const Tensor& other ) const
  {
#ifndef WIN32
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#endif

    return x == other.x && y == other.y && z == other.z;

#ifndef WIN32
#pragma GCC diagnostic pop
#endif
  }
};

std::ostream& operator<<(std::ostream& stream, const Tensor & t )
{
  stream << "(" << t.x << ", " << t.y << ", " << t.z << ")";
  return stream;
}

TEST( StackArray, storage )
{
  constexpr int N = 100;     /* Number of values to store */
  internal::test_storage< int, N >(
    []( int i ) {
      return i;
    } );

  internal::test_storage< Tensor, N >(
    []( int i ) {
      return Tensor( i );
    } );

  internal::test_storage< std::string, N >(
    []( int i ) {
      return std::to_string( i );
    } );
}

TEST( StackArray, copy )
{
  constexpr int N = 100;     /* Number of values to store */
  internal::test_copy< int, N >(
    []( int i ) {
      return i;
    } );

  internal::test_copy< Tensor, N >(
    []( int i ) {
      return Tensor( i );
    } );

  internal::test_copy< std::string, N >(
    []( int i ) {
      return std::to_string( i );
    } );
}

TEST( StackArray, list_initilization )
{
  internal::test_list_initilization< int >(
    []( int i ) {
      return i;
    } );

  internal::test_list_initilization< Tensor >(
    []( int i ) {
      return Tensor( i );
    } );

  internal::test_list_initilization< std::string >(
    []( int i ) {
      return std::to_string( i );
    } );
}

} /* namespace axom */

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}

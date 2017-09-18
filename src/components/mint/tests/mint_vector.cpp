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
#include "axom_utils/Utilities.hpp"     // for utilities::max
#include "gtest/gtest.h"                // for TEST and EXPECT_* macros
#include "slic/slic.hpp"                // for slic macros
#include "slic/UnitTestLogger.hpp"      // for UnitTestLogger


#include <cmath>                        // for std::ceil


namespace axom {
namespace mint {


/*!
 * \brief Creates a 1D ParticleMesh and writes it out to disk using
 *  mint::write_vtk and then reads the file back in to check for correctness.
 */
TEST( mint_vector, checkStorage ) {
  const int capacity = 1000;

  Vector< int, int > vii( capacity );
  Vector< double, int> vdi( capacity );
  Vector< int, long long uint> viu( capacity );
  Vector< double, long long uint> vdu( capacity );

  EXPECT_TRUE( vii.empty() );
  EXPECT_TRUE( vdi.empty() );
  EXPECT_TRUE( viu.empty() );
  EXPECT_TRUE( vdi.empty() );


  EXPECT_EQ( vii.getSize(), 0 );
  EXPECT_EQ( vdi.getSize(), 0 );
  EXPECT_EQ( viu.getSize(), 0 );
  EXPECT_EQ( vdu.getSize(), 0 );
  
  EXPECT_EQ( vii.getCapacity(), capacity );
  EXPECT_EQ( vdi.getCapacity(), capacity );
  EXPECT_EQ( viu.getCapacity(), capacity );
  EXPECT_EQ( vdu.getCapacity(), capacity );
  
  EXPECT_EQ( vii.getResizeRatio(), 2.0 );
  EXPECT_EQ( vdi.getResizeRatio(), 2.0 ); 
  EXPECT_EQ( viu.getResizeRatio(), 2.0 ); 
  EXPECT_EQ( vdu.getResizeRatio(), 2.0 );

  EXPECT_EQ( vii.getChunkSize(), 1 );
  EXPECT_EQ( vdi.getChunkSize(), 1 ); 
  EXPECT_EQ( viu.getChunkSize(), 1 ); 
  EXPECT_EQ( vdu.getChunkSize(), 1 ); 

  for ( int i = 0; i < capacity; ++i ) {
    vii.add( i );
    vdi.add( i );
    viu.add( i );
    vdu.add( i );
  }

  EXPECT_TRUE( !vii.empty() );
  EXPECT_TRUE( !vdi.empty() );
  EXPECT_TRUE( !viu.empty() );
  EXPECT_TRUE( !vdi.empty() );

  EXPECT_EQ( vii.getSize(), capacity );
  EXPECT_EQ( vdi.getSize(), capacity );
  EXPECT_EQ( viu.getSize(), capacity );
  EXPECT_EQ( vdu.getSize(), capacity );
  
  EXPECT_EQ( vii.getCapacity(), capacity );
  EXPECT_EQ( vdi.getCapacity(), capacity );
  EXPECT_EQ( viu.getCapacity(), capacity );
  EXPECT_EQ( vdu.getCapacity(), capacity );

  for ( int i = 0; i < capacity; ++i ) {
    EXPECT_EQ( vii[ i ], i );
    EXPECT_EQ( vdi[ i ], i );
    EXPECT_EQ( viu[ i ], i );
    EXPECT_EQ( vdu[ i ], i );
  }

  for ( int i = 0; i < capacity; ++i ) {
    vii[ i ] = capacity - i;
    vdi[ i ] = capacity - i;
    viu[ i ] = capacity - i;
    vdu[ i ] = capacity - i;
  }

  for ( int i = 0; i < capacity; ++i ) {
    EXPECT_EQ( vii[ i ], capacity - i );
    EXPECT_EQ( vdi[ i ], capacity - i );
    EXPECT_EQ( viu[ i ], capacity - i );
    EXPECT_EQ( vdu[ i ], capacity - i );
  }
}

TEST( mint_vector, checkResize ) {
  int capacity = 0;
  int size = 0;
  int chunk_size = 1;
  int resize_ratio = 2.0;
  Vector< double, int> v( capacity, chunk_size, resize_ratio );

  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.getSize(), size );
  EXPECT_EQ( v.getChunkSize(), chunk_size );
  EXPECT_EQ( v.getResizeRatio(), resize_ratio );

  v.add( 0.0 );
  size += 1;
  capacity = resize_ratio * size;

  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.getSize(), size );
  EXPECT_EQ( v[0], 0.0 );

  const int numValues = 1000;
  double values[ numValues ];
  for ( int i = 0; i < numValues; ++i ) {
    values[ i ] = i + 1; 
  }

  v.add( values, numValues );
  size += numValues;
  capacity = resize_ratio * size;

  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.getSize(), size );
  for ( int i = 0; i < size; ++i ) {
    EXPECT_EQ( v[ i ], i );
  }

  v.add( numValues + 1 );
  size += 1;
  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.getSize(), size );
  EXPECT_EQ( v[ size - 1 ], numValues + 1 );

  size = 500;
  v.setSize(500);
  EXPECT_EQ( v.getSize(), size );
  EXPECT_EQ( v.getCapacity(), capacity );
  for ( int i = 0; i < size; ++i ) {
    EXPECT_EQ( v[ i ], i );
  }

  capacity = 250;
  size = capacity;
  v.setCapacity( capacity );
  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.getSize(), size );
  for ( int i = 0; i < size; ++i ) {
    EXPECT_EQ( v[ i ], i );
  }

  size += numValues;
  capacity = size;
  v.setSize( size );
  EXPECT_EQ( v.getSize(), size );
  EXPECT_EQ( v.getCapacity(), capacity );

  v.set( values, numValues, size - numValues );
  EXPECT_EQ( v.getSize(), size );
  EXPECT_EQ( v.getCapacity(), capacity );

  for ( int i = 0 ; i < 250; ++i ) {
    EXPECT_EQ( v[ i ], i );
  }

  for ( int i = 250; i < size; ++i ) {
    EXPECT_EQ( v[ i ], i - 250 + 1 );
  }
}

/* Do the chunk stuff. */
TEST( mint_vector, checkResizeChunk ) {
  int capacity = 0;
  int size = 0;
  int chunk_size = 7;
  int resize_ratio = 2.0;
  Vector< double, int> v( capacity, chunk_size, resize_ratio );

  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.getSize(), size );
  EXPECT_EQ( v.getChunkSize(), chunk_size );
  EXPECT_EQ( v.getResizeRatio(), resize_ratio );

  v.add( 0.0 );
  size += 1;
  capacity = std::ceil( double( size ) / chunk_size ) * chunk_size;


  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.getSize(), size );
  EXPECT_EQ( v[0], 0.0 );

  const int numValues = 1000;
  double values[ numValues ];
  for ( int i = 0; i < numValues; ++i ) {
    values[ i ] = i + 1; 
  }

  v.add( values, numValues );
  size += numValues;
  capacity = std::ceil( double( size * resize_ratio ) / chunk_size );
  capacity *= chunk_size;

  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.getSize(), size );
  for ( int i = 0; i < size; ++i ) {
    EXPECT_EQ( v[ i ], i );
  }

  v.add( numValues + 1 );
  size += 1;
  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.getSize(), size );
  EXPECT_EQ( v[ size - 1 ], numValues + 1 );

  size = 500;
  v.setSize(500);
  EXPECT_EQ( v.getSize(), size );
  EXPECT_EQ( v.getCapacity(), capacity );
  for ( int i = 0; i < size; ++i ) {
    EXPECT_EQ( v[ i ], i );
  }

  capacity = 250;
  v.setCapacity( capacity );
  capacity = std::ceil( double( capacity ) / chunk_size ) * chunk_size;
  size = capacity;

  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.getSize(), size );
  for ( int i = 0; i < size; ++i ) {
    EXPECT_EQ( v[ i ], i );
  }

  size += numValues;
  capacity = std::ceil( double( size ) / chunk_size ) * chunk_size;
  v.setSize( size );
  EXPECT_EQ( v.getSize(), size );
  EXPECT_EQ( v.getCapacity(), capacity );

  v.set( values, numValues, size - numValues );
  EXPECT_EQ( v.getSize(), size );
  EXPECT_EQ( v.getCapacity(), capacity );

  for ( int i = 0 ; i < size - numValues; ++i ) {
    EXPECT_EQ( v[ i ], i );
  }

  for ( int i = size - numValues; i < size; ++i ) {
    EXPECT_EQ( v[ i ], i + numValues - size + 1 );
  }
}

/* test reserve */
/* test insert */


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

} /* end namespace mint */
} /* end namespace axom */

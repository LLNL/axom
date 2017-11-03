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
#include "mint/DataTypes.hpp"
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
  const localIndex capacity = 1000;

  Vector< int > vi( capacity );
  Vector< double > vd( capacity );

  EXPECT_TRUE( vi.empty() );
  EXPECT_TRUE( vd.empty() );

  EXPECT_EQ( vi.getSize(), 0 );
  EXPECT_EQ( vd.getSize(), 0 );

  EXPECT_EQ( vi.getCapacity(), capacity );
  EXPECT_EQ( vd.getCapacity(), capacity );

  EXPECT_EQ( vi.getResizeRatio(), 2.0 );
  EXPECT_EQ( vd.getResizeRatio(), 2.0 ); 

  EXPECT_EQ( vi.getNumComponents(), 1 );
  EXPECT_EQ( vd.getNumComponents(), 1 ); 

  for ( localIndex i = 0; i < capacity; ++i ) {
    vi.add( i );
    vd.add( i );
  }

  EXPECT_TRUE( !vi.empty() );
  EXPECT_TRUE( !vd.empty() );

  EXPECT_EQ( vi.getSize(), capacity );
  EXPECT_EQ( vd.getSize(), capacity );

  EXPECT_EQ( vi.getCapacity(), capacity );
  EXPECT_EQ( vd.getCapacity(), capacity );


  for ( localIndex i = 0; i < capacity; ++i ) {
    EXPECT_EQ( vi[ i ], i );
    EXPECT_EQ( vd[ i ], i );
  }

  for ( localIndex i = 0; i < capacity; ++i ) {
    vi[ i ] = capacity - i;
    vd[ i ] = capacity - i;
  }

  for ( localIndex i = 0; i < capacity; ++i ) {
    EXPECT_EQ( vi[ i ], capacity - i );
    EXPECT_EQ( vd[ i ], capacity - i );
  }
}

TEST( mint_vector, checkResize ) {
  localIndex capacity = 0;
  localIndex size = 0;
  int num_components = 1;
  double resize_ratio = 2.0;
  Vector< double > v( capacity, num_components, resize_ratio );

  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.getSize(), size );
  EXPECT_EQ( v.getNumComponents(), num_components );
  EXPECT_EQ( v.getResizeRatio(), resize_ratio );

  v.add( 0.0 );
  size += 1;
  capacity = resize_ratio * size;

  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.getSize(), size );
  EXPECT_EQ( v[0], 0.0 );

  const localIndex numValues = 1000;
  double values[ numValues ];
  for ( localIndex i = 0; i < numValues; ++i ) {
    values[ i ] = i + 1; 
  }

  v.add( values, numValues );
  size += numValues;
  capacity = resize_ratio * size;

  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.getSize(), size );
  for ( localIndex i = 0; i < size; ++i ) {
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
  for ( localIndex i = 0; i < size; ++i ) {
    EXPECT_EQ( v[ i ], i );
  }

  capacity = 250;
  size = capacity;
  v.setCapacity( capacity );
  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.getSize(), size );
  for ( localIndex i = 0; i < size; ++i ) {
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

  for ( localIndex i = 0 ; i < 250; ++i ) {
    EXPECT_EQ( v[ i ], i );
  }

  for ( localIndex i = 250; i < size; ++i ) {
    EXPECT_EQ( v[ i ], i - 250 + 1 );
  }
}

/* Do the chunk stuff. */
TEST( mint_vector, checkResizeMultiComponent ) {
  localIndex capacity = 0;
  localIndex size = 0;
  int num_components = 7;
  double resize_ratio = 2.0;
  Vector< double > v( capacity, num_components, resize_ratio );

  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.getSize(), size );
  EXPECT_EQ( v.getNumComponents(), num_components );
  EXPECT_EQ( v.getResizeRatio(), resize_ratio );

  v.add( 0.0 );
  size += 1;
  capacity = std::ceil( double( size ) / num_components ) * num_components;


  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.getSize(), size );
  EXPECT_EQ( v[0], 0.0 );

  const localIndex numValues = 1000;
  double values[ numValues ];
  for ( localIndex i = 0; i < numValues; ++i ) {
    values[ i ] = i + 1; 
  }

  v.add( values, numValues );
  size += numValues;
  capacity = std::ceil( double( size * resize_ratio ) / num_components );
  capacity *= num_components;

  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.getSize(), size );
  for ( localIndex i = 0; i < size; ++i ) {
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
  for ( localIndex i = 0; i < size; ++i ) {
    EXPECT_EQ( v[ i ], i );
  }

  capacity = 250;
  v.setCapacity( capacity );
  capacity = std::ceil( double( capacity ) / num_components ) * num_components;
  size = capacity;

  EXPECT_EQ( v.getCapacity(), capacity );
  EXPECT_EQ( v.getSize(), size );
  for ( localIndex i = 0; i < size; ++i ) {
    EXPECT_EQ( v[ i ], i );
  }

  size += numValues;
  capacity = std::ceil( double( size ) / num_components ) * num_components;
  v.setSize( size );
  EXPECT_EQ( v.getSize(), size );
  EXPECT_EQ( v.getCapacity(), capacity );

  v.set( values, numValues, size - numValues );
  EXPECT_EQ( v.getSize(), size );
  EXPECT_EQ( v.getCapacity(), capacity );

  for ( localIndex i = 0 ; i < size - numValues; ++i ) {
    EXPECT_EQ( v[ i ], i );
  }

  for ( localIndex i = size - numValues; i < size; ++i ) {
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

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

#include "primal/Plane.hpp"
#include "axom_utils/vector_utilities.hpp"

#include "gtest/gtest.h"

// C/C++ includes
#include <cmath>

namespace primal   = axom::primal;
namespace numerics = axom::numerics;

//------------------------------------------------------------------------------
// INTERNAL HELPER METHODS
//------------------------------------------------------------------------------
namespace
{

void ensure_unit_norm( const double* v, int n )
{
  const double norm = std::sqrt( numerics::dot_product(v,v,n) );
  EXPECT_DOUBLE_EQ( norm, 1.0 );
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
// UNIT TEST
//------------------------------------------------------------------------------

TEST( primal_plane, construct_from_normal_and_point )
{
  double normal[ 3 ] = { 0.0, 0.0, 10.0 };
  double x[ 3 ]      = { 0.0, 0.0, 2.0 };

  // test 3D
  primal::Plane< double,3 > P( normal, x );
  EXPECT_DOUBLE_EQ( P.getOffset(), 2.0 );
  EXPECT_DOUBLE_EQ( P.getNormal()[0], 0.0 );
  EXPECT_DOUBLE_EQ( P.getNormal()[1], 0.0 );
  EXPECT_DOUBLE_EQ( P.getNormal()[2], 1.0 );
  EXPECT_EQ( P.getDimension(), 3 );

  // test 2D
  normal[ 0 ] = 1.0; normal[ 1 ] = 2.0;
  x[ 0 ]      = 1.0; x[ 1 ]      = 2.0;

  primal::Plane< double, 2 > P2( normal, x );
  ensure_unit_norm( P2.getNormal(), 2 );
  EXPECT_DOUBLE_EQ( P2.getOffset(), std::sqrt( 5.0 ) );
  EXPECT_EQ( P2.getDimension(), 2 );
}

//------------------------------------------------------------------------------
TEST( primal_plane, construct_from_normal_and_offset )
{
  double normal[ 3 ] = { 0.0, 0.0, 1.0 };
  double offset      = 2.0;

  // test 3D
  primal::Plane< double, 3 > P( normal, offset, true );
  EXPECT_DOUBLE_EQ( P.getNormal()[0], 0.0 );
  EXPECT_DOUBLE_EQ( P.getNormal()[1], 0.0 );
  EXPECT_DOUBLE_EQ( P.getNormal()[2], 1.0 );
  EXPECT_DOUBLE_EQ( P.getOffset(), offset );

  // test 2D
  normal[ 0 ] = 1.0; normal[ 1 ] = 2.0;
  offset = std::sqrt( 5.0 );
  primal::Plane< double,2 > P2( normal, offset );
  ensure_unit_norm( P2.getNormal(), 2 );
  EXPECT_DOUBLE_EQ( P2.getOffset(), offset );
}

//------------------------------------------------------------------------------
TEST( primal_plane, construct_from_points )
{
  double x1[ 3 ] = { 1.0, 1.0, 3.0 };
  double x2[ 3 ] = { 2.0, 2.0, 3.0 };
  double x3[ 3 ] = { 1.0, 3.0, 3.0 };

  // test 3D
  primal::Plane< double, 3 > P( x1, x2, x3 );
  ensure_unit_norm( P.getNormal(), 3 );
  EXPECT_DOUBLE_EQ( P.getOffset(), 3.0 );

  // test 2D
  double a[ 2 ] = { 2.0, -1.0 };
  double b[ 2 ] = { 2.0,  2.0 };
  primal::Plane< double, 2 > P2( a, b, AXOM_NULLPTR );
  ensure_unit_norm( P2.getNormal(), 2 );
  EXPECT_DOUBLE_EQ( P2.getOffset(), -2.0 );
}

//------------------------------------------------------------------------------
TEST( primal_plane, signed_distance_and_orientation )
{
  double x1[ 3 ] = { 1.0, 1.0, 3.0 };
  double x2[ 3 ] = { 2.0, 2.0, 3.0 };
  double x3[ 3 ] = { 1.0, 3.0, 3.0 };

  double signed_distance = 0.0;       // stores computed signed distance.
  double q[ 3 ]  = { 0.0, 0.0, 0.0 }; // test query point

  // STEP 0: test 3D
  primal::Plane< double, 3 > P( x1, x2, x3 );

  // (a) test point below plane
  signed_distance = P.computeSignedDistance( q );
  EXPECT_DOUBLE_EQ( signed_distance, -3.0 );
  EXPECT_EQ( P.getOrientedSide( q ), primal::ON_NEGATIVE_SIDE );

  // (b) test point above plane
  q[ 2 ] = 6.0;
  signed_distance = P.computeSignedDistance( q );
  EXPECT_DOUBLE_EQ( signed_distance, 3.0 );
  EXPECT_EQ( P.getOrientedSide( q ), primal::ON_POSITIVE_SIDE );

  // (c) test point on plane
  q[ 2 ] = 3.0;
  signed_distance = P.computeSignedDistance( q );
  EXPECT_DOUBLE_EQ( signed_distance, 0.0 );
  EXPECT_EQ( P.getOrientedSide( q ), primal::ON_BOUNDARY );

  // STEP 1: test 2D
  double a[ 2 ] = { 2.0, -1.0 };
  double b[ 2 ] = { 2.0,  2.0 };
  primal::Plane< double, 2 > P2( a, b, AXOM_NULLPTR );

  // (a) test point above plane
  signed_distance = P2.computeSignedDistance( q );
  EXPECT_DOUBLE_EQ( signed_distance, 2.0 );
  EXPECT_EQ( P2.getOrientedSide( q ), primal::ON_POSITIVE_SIDE );

  // (b) test point below plane
  q[ 0 ] = 4.0;
  signed_distance = P2.computeSignedDistance( q );
  EXPECT_DOUBLE_EQ( signed_distance, -2.0 );
  EXPECT_EQ( P2.getOrientedSide( q ), primal::ON_NEGATIVE_SIDE );

  // (c) test point on plane
  q[ 0 ] = 2.0;
  signed_distance = P2.computeSignedDistance( q );
  EXPECT_DOUBLE_EQ( signed_distance, 0.0 );
  EXPECT_EQ( P2.getOrientedSide( q ), primal::ON_BOUNDARY );
}

//------------------------------------------------------------------------------
TEST( primal_plane, project_point)
{
  double x1[ 3 ] = { 1.0, 1.0, 3.0 };
  double x2[ 3 ] = { 2.0, 2.0, 3.0 };
  double x3[ 3 ] = { 1.0, 3.0, 3.0 };

  double q[ 3 ]     = { 0.0, 0.0, 0.0 };
  double qproj[ 3 ] = { 0.0, 0.0, 0.0 };

  // STEP 0: test 3D
  primal::Plane< double, 3 > P( x1, x2, x3 );

  // (a) test project point below plane
  P.projectPoint( q, qproj );
  EXPECT_EQ( P.getOrientedSide( qproj ), primal::ON_BOUNDARY );
  EXPECT_DOUBLE_EQ( qproj[ 0 ], 0.0 );
  EXPECT_DOUBLE_EQ( qproj[ 1 ], 0.0 );
  EXPECT_DOUBLE_EQ( qproj[ 2 ], 3.0 );

  // (b) test project point above plane
  q[ 2 ] = 6.0;
  qproj[ 0 ] = qproj[ 1 ] = qproj[ 2 ] = 0.0;
  P.projectPoint( q, qproj );
  EXPECT_EQ( P.getOrientedSide( qproj ), primal::ON_BOUNDARY );
  EXPECT_DOUBLE_EQ( qproj[ 0 ], 0.0 );
  EXPECT_DOUBLE_EQ( qproj[ 1 ], 0.0 );
  EXPECT_DOUBLE_EQ( qproj[ 2 ], 3.0 );

  // (c) test project point (already) on plane
  q[ 2 ] = 3.0;
  qproj[ 0 ] = qproj[ 1 ] = qproj[ 2 ] = 0.0;
  P.projectPoint( q, qproj );
  EXPECT_EQ( P.getOrientedSide( qproj ), primal::ON_BOUNDARY );
  EXPECT_DOUBLE_EQ( qproj[ 0 ], q[ 0 ] );
  EXPECT_DOUBLE_EQ( qproj[ 1 ], q[ 1 ] );
  EXPECT_DOUBLE_EQ( qproj[ 2 ], q[ 2 ] );

  // STEP 1: test 2D
  double a[ 2 ] = { 2.0, -1.0 };
  double b[ 2 ] = { 2.0,  2.0 };
  primal::Plane< double, 2 > P2( a, b, AXOM_NULLPTR );

  // (a) test project point below plane
  q[ 0 ] = 4.0;
  qproj[ 0 ] = qproj[ 1 ] = qproj[ 2 ] = 0.0;
  P2.projectPoint( q, qproj );
  EXPECT_EQ( P2.getOrientedSide( qproj ), primal::ON_BOUNDARY );
  EXPECT_DOUBLE_EQ( qproj[ 0 ], 2.0 );
  EXPECT_DOUBLE_EQ( qproj[ 1 ], 0.0 );

  // (b) test project point above plane
  q[ 0 ] = 0.0;
  qproj[ 0 ] = qproj[ 1 ] = qproj[ 2 ] = 0.0;
  P2.projectPoint( q, qproj );
  EXPECT_EQ( P2.getOrientedSide( qproj ), primal::ON_BOUNDARY );
  EXPECT_DOUBLE_EQ( qproj[ 0 ], 2.0 );
  EXPECT_DOUBLE_EQ( qproj[ 1 ], 0.0 );

  // (c) test project point (already) on plane
  q[ 0 ] = 2.0;
  qproj[ 0 ] = qproj[ 1 ] = qproj[ 2 ] = 0.0;
  P2.projectPoint( q, qproj );
  EXPECT_EQ( P2.getOrientedSide( qproj ), primal::ON_BOUNDARY );
  EXPECT_DOUBLE_EQ( qproj[ 0 ], q[ 0 ] );
  EXPECT_DOUBLE_EQ( qproj[ 1 ], q[ 1 ] );
}

//------------------------------------------------------------------------------
TEST( primal_plane, flip )
{
  double x1[ 3 ] = { 1.0, 1.0, 3.0 };
  double x2[ 3 ] = { 2.0, 2.0, 3.0 };
  double x3[ 3 ] = { 1.0, 3.0, 3.0 };

  double q[ 3 ] = { 0.0, 0.0, 0.0 };

  // STEP 0: test 3D
  primal::Plane< double, 3 > P( x1, x2, x3 );
  EXPECT_DOUBLE_EQ( P.getOrientedSide( q ), primal::ON_NEGATIVE_SIDE );
  P.flip( );
  EXPECT_DOUBLE_EQ( P.getOrientedSide( q ), primal::ON_POSITIVE_SIDE );

  // STEP 1: test 2D
  double a[ 2 ] = { 2.0, -1.0 };
  double b[ 2 ] = { 2.0,  2.0 };
  primal::Plane< double, 2 > P2( a, b, AXOM_NULLPTR );
  EXPECT_DOUBLE_EQ( P2.getOrientedSide( q ), primal::ON_POSITIVE_SIDE );
  P2.flip( );
  EXPECT_DOUBLE_EQ( P2.getOrientedSide( q ), primal::ON_NEGATIVE_SIDE );
}

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

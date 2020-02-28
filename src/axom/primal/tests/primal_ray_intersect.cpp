// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/primal/operators/detail/intersect_ray_impl.hpp"

#include "axom/core/numerics/Matrix.hpp"     // for numerics::Matrix
#include "axom/core/numerics/matvecops.hpp"  // for matrix-vector operators

#include "gtest/gtest.h"

// C/C++ includes
#include <cmath>

// namespace aliases
namespace primal   = axom::primal;
namespace numerics = axom::numerics;

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST( primal_ray_intersect, ray_aabb_intersection_2D )
{
  constexpr int N     = 3;
  constexpr double LO = 0.0f;
  constexpr double HI = 1.0f;

// Defines the test box:
#define TEST_BOX2D LO,HI,LO,HI

  double x[ N ];
  numerics::linspace( LO, HI, x, N );

  // test BOTTOM (-y)
  {
    const double n0[] = { 0.0,  1.0 };
    const double n1[] = { 0.0, -1.0 };
    const double y0   = -1.0f;
    for ( int i=0; i < N; ++i )
    {
      double tmin = 0.0;
      double tmax = 0.0;

      double s[] = { x[i], y0 };
      EXPECT_TRUE( primal::detail::intersect_ray(s,n0,TEST_BOX2D,tmin,tmax) );
      EXPECT_FALSE( primal::detail::intersect_ray(s,n1,TEST_BOX2D,tmin,tmax) );
    }
  }

  // test RIGHT (+x)
  {
    const double n0[] = { -1.0, 0.0 };
    const double n1[] = {  1.0, 0.0 };
    const double x0   = 2.0f;
    for ( int i=0; i < N; ++i )
    {
      double tmin = 0.0;
      double tmax = 0.0;

      double s[] = { x0, x[ i ] };
      EXPECT_TRUE( primal::detail::intersect_ray(s,n0,TEST_BOX2D,tmin,tmax) );
      EXPECT_FALSE( primal::detail::intersect_ray(s,n1,TEST_BOX2D,tmin,tmax) );
    }
  }

  // test TOP (+y)
  {
    const double n0[] = { 0.0, -1.0 };
    const double n1[] = { 0.0,  1.0 };
    const double y0   = 2.0f;
    for ( int i=0; i < N; ++i )
    {
      double tmin = 0.0;
      double tmax = 0.0;

      double s[] = { x[ i ], y0 };
      EXPECT_TRUE( primal::detail::intersect_ray(s,n0,TEST_BOX2D,tmin,tmax) );
      EXPECT_FALSE( primal::detail::intersect_ray(s,n1,TEST_BOX2D,tmin,tmax) );
    }
  }

  // test LEFT (-x)
  {
    const double n0[] = {  1.0, 0.0 };
    const double n1[] = { -1.0, 0.0 };
    const double x0   = -1.0f;
    for ( int i=0; i < N; ++i )
    {
      double tmin = 0.0;
      double tmax = 0.0;

      double s[] = { x0, x[ i ] };
      EXPECT_TRUE( primal::detail::intersect_ray(s,n0,TEST_BOX2D,tmin,tmax) );
      EXPECT_FALSE( primal::detail::intersect_ray(s,n1,TEST_BOX2D,tmin,tmax) );
    }
  }

  // Test a bunch of rays emitted from the box center
  constexpr int NUM_ANGLES = 20;
  double angles[ NUM_ANGLES ];
  numerics::linspace( 0.0, 360.0, angles, NUM_ANGLES );

  numerics::Matrix< double > rotation_matrix( 2, 2 );

  constexpr double PI_OVER_180 = M_PI / 180.0;
  const double xc[] = { 0.5, 0.5 };
  const double e1[] = { 1.0, 0.0 };
  for( int i=0; i < NUM_ANGLES; ++i )
  {
    const double t    = angles[ i ] * PI_OVER_180;
    const double cost = cos( t );
    const double sint = sin( t );

    rotation_matrix( 0,0 ) = cost; rotation_matrix( 0,1 ) = -sint;
    rotation_matrix( 1,0 ) = sint; rotation_matrix( 1,1 ) = cost;

    double n[ 2 ];
    numerics::matrix_vector_multiply( rotation_matrix, e1, n );

    double inv_n[ ] = { -n[ 0 ], -n[ 1 ] };

    double tmin = 0.0;
    double tmax = 0.0;
    EXPECT_TRUE( primal::detail::intersect_ray(xc,n,TEST_BOX2D,tmin,tmax) );
    EXPECT_TRUE( primal::detail::intersect_ray(xc,inv_n,TEST_BOX2D,tmin,tmax) );
  }

#undef TEST_BOX2D
}

//------------------------------------------------------------------------------
TEST( primal_ray_intersect, ray_aabb_intersection_3D )
{
  constexpr int N     = 3;
  constexpr double LO = 0.0f;
  constexpr double HI = 1.0f;

// Defines the test box:
#define TEST_BOX3D LO,HI,LO,HI,LO,HI

  double x[ N ];
  numerics::linspace( LO, HI, x, N );

  // test BOTTOM (-z)
  {

    const double n0[] = { 0.0, 0.0,  1.0 };
    const double n1[] = { 0.0, 0.0, -1.0 };
    const double z0   = -1.0f;

    for ( int i=0; i < N; ++i )
    {
      for ( int j=0; j < N; ++j )
      {
        double tmin = 0.0;
        double tmax = 0.0;

        double s[] = { x[i], x[j], z0 };
        EXPECT_TRUE( primal::detail::intersect_ray(s,n0,TEST_BOX3D,tmin,tmax));
        EXPECT_FALSE( primal::detail::intersect_ray(s,n1,TEST_BOX3D,tmin,tmax));
      } // END for all j
    } // END for all i

  }

  // test TOP (+z)
  {

    const double n0[] = { 0.0, 0.0, -1.0 };
    const double n1[] = { 0.0, 0.0, 1.0 };
    const double z0   = 2.0f;

    for ( int i=0; i < N; ++i )
    {
      for ( int j=0; j < N; ++j )
      {
        double tmin = 0.0;
        double tmax = 0.0;

        double s[] = { x[i], x[j], z0 };
        EXPECT_TRUE( primal::detail::intersect_ray(s,n0,TEST_BOX3D,tmin,tmax));
        EXPECT_FALSE( primal::detail::intersect_ray(s,n1,TEST_BOX3D,tmin,tmax));
      } // END for all j
    } // END for all i

  }

  // test LEFT (-y)
  {

    const double n0[] = { 0.0,  1.0, 0.0 };
    const double n1[] = { 0.0, -1.0, 0.0 };
    const double y0   = -1.0f;

    for ( int i=0; i < N; ++i )
    {
      for ( int j=0; j < N; ++j )
      {
        double tmin = 0.0;
        double tmax = 0.0;

        double s[] = { x[ i ], y0, x[ j ] };
        EXPECT_TRUE( primal::detail::intersect_ray(s,n0,TEST_BOX3D,tmin,tmax));
        EXPECT_FALSE( primal::detail::intersect_ray(s,n1,TEST_BOX3D,tmin,tmax));
      } // END for all j
    } // END for all i

  }

  // test RIGHT (+y)
  {

    const double n0[] = { 0.0, -1.0, 0.0 };
    const double n1[] = { 0.0,  1.0, 0.0 };
    const double y0   = 2.0f;

    for ( int i=0; i < N; ++i )
    {
      for ( int j=0; j < N; ++j )
      {
        double tmin = 0.0;
        double tmax = 0.0;

        double s[] = { x[ i ], y0, x[ j ] };
        EXPECT_TRUE( primal::detail::intersect_ray(s,n0,TEST_BOX3D,tmin,tmax));
        EXPECT_FALSE( primal::detail::intersect_ray(s,n1,TEST_BOX3D,tmin,tmax));
      } // END for all j
    } // END for all i

  }

  // test BACK (-x)
  {

    const double n0[] = {  1.0, 0.0, 0.0 };
    const double n1[] = { -1.0, 0.0, 0.0 };
    const double x0   = -1.0f;

    for ( int i=0; i < N; ++i )
    {
      for ( int j=0; j < N; ++j )
      {
        double tmin = 0.0;
        double tmax = 0.0;

        double s[] = { x0, x[ i ], x[ j ] };
        EXPECT_TRUE( primal::detail::intersect_ray(s,n0,TEST_BOX3D,tmin,tmax));
        EXPECT_FALSE( primal::detail::intersect_ray(s,n1,TEST_BOX3D,tmin,tmax));
      } // END for all j
    } // END for all i

  }

  // test FRONT (+x)
  {

    const double n0[] = { -1.0, 0.0, 0.0 };
    const double n1[] = {  1.0, 0.0, 0.0 };
    const double x0   = 2.0f;

    for ( int i=0; i < N; ++i )
    {
      for ( int j=0; j < N; ++j )
      {
        double tmin = 0.0;
        double tmax = 0.0;

        double s[] = { x0, x[ i ], x[ j ] };
        EXPECT_TRUE( primal::detail::intersect_ray(s,n0,TEST_BOX3D,tmin,tmax));
        EXPECT_FALSE( primal::detail::intersect_ray(s,n1,TEST_BOX3D,tmin,tmax));
      } // END for all j
    } // END for all i

  }

  // Test a bunch of rays emitted from the box center
  constexpr int NUM_ANGLES = 20;
  double angles[ NUM_ANGLES ];
  numerics::linspace( 0.0, 360.0, angles, NUM_ANGLES );

  numerics::Matrix< double > Rx(3,3); Rx(0,0) = 1.0;
  numerics::Matrix< double > Ry(3,3); Ry(1,1) = 1.0;
  numerics::Matrix< double > Rz(3,3); Rz(2,2) = 1.0;

  constexpr double PI_OVER_180 = M_PI / 180.0;
  const double xc[] = { 0.5, 0.5, 0.5 };
  const double e1[] = { 1.0, 0.0, 0.0 };
  const double e2[] = { 0.0, 1.0, 0.0 };
  const double e3[] = { 0.0, 0.0, 1.0 };
  for( int i=0; i < NUM_ANGLES; ++i )
  {
    const double t    = angles[ i ] * PI_OVER_180;
    const double cost = cos( t );
    const double sint = sin( t );

    double tmin = 0.0;
    double tmax = 0.0;

    // Update Rx
    Rx(1,1)=cost; Rx(1,2)= -sint;
    Rx(2,1)=sint; Rx(2,2)= cost;

    double nx[ 3 ];
    numerics::matrix_vector_multiply( Rx, e1, nx );
    EXPECT_TRUE( primal::detail::intersect_ray(xc,nx,TEST_BOX3D,tmin,tmax) );

    // Update Ry
    Ry(0,0)=cost;  Ry(0,2)=sint;
    Ry(2,0)=-sint; Ry(2,2)=cost;

    double ny[ 3 ];
    numerics::matrix_vector_multiply( Ry, e2, ny );
    EXPECT_TRUE( primal::detail::intersect_ray(xc,ny,TEST_BOX3D,tmin,tmax) );

    // Update Rz
    Rz(0,0)=cost; Rz(0,1)=-sint;
    Rz(1,0)=sint; Rz(1,1)=cost;

    double nz[ 3 ];
    numerics::matrix_vector_multiply( Rz, e3, nz );
    EXPECT_TRUE( primal::detail::intersect_ray(xc,nz,TEST_BOX3D,tmin,tmax) );
  }



#undef TEST_BOX3D
}

//------------------------------------------------------------------------------
#include "axom/slic/core/UnitTestLogger.hpp"
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




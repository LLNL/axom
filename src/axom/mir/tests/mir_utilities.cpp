// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MIR_UTILITIES_TEST_H_
#define MIR_UTILITIES_TEST_H_

#include "gtest/gtest.h"

#include "axom/slic.hpp"
#include "axom/mir.hpp"

using namespace axom;

TEST(mir_interpolation, float_linear_interpolation)
{
  axom::float64 f0 = 50.0;
  axom::float64 f1 = 100.0;

  axom::float64 t = 0.0;
  axom::float64 interpolant = mir::utilities::lerpFloat(f0, f1, t);
  EXPECT_DOUBLE_EQ(  interpolant,  50.0  );

  t = 1.0;
  interpolant = mir::utilities::lerpFloat(f0, f1, t);
  EXPECT_DOUBLE_EQ(  interpolant,  100.0  );
  
  t = 0.5;
  interpolant = mir::utilities::lerpFloat(f0, f1, t);
  EXPECT_DOUBLE_EQ(  interpolant,  75.0  );

  t = 0.66;
  interpolant = mir::utilities::lerpFloat(f0, f1, t);
  EXPECT_DOUBLE_EQ(  interpolant,  83.0  );
}

//----------------------------------------------------------------------

TEST(mir_interpolation, point_linear_interpolation)
{
  mir::Point2 p0(0.25, 0.0);
  mir::Point2 p1(1.25, 1.0);

  axom::float64 t = 0.0;
  mir::Point2 interpolant = mir::utilities::interpolateVertexPosition(p0, p1, t);
  EXPECT_DOUBLE_EQ(  interpolant.m_x, 0.25);
  EXPECT_DOUBLE_EQ(  interpolant.m_y, 0.0);

  t = 1.0;
  interpolant = mir::utilities::interpolateVertexPosition(p0, p1, t);
  EXPECT_DOUBLE_EQ(  interpolant.m_x, 1.25);
  EXPECT_DOUBLE_EQ(  interpolant.m_y, 1.0);

  t = 0.5;
  interpolant = mir::utilities::interpolateVertexPosition(p0, p1, t);
  EXPECT_DOUBLE_EQ(  interpolant.m_x, 0.75);
  EXPECT_DOUBLE_EQ(  interpolant.m_y, 0.5);

  t = 0.66;
  interpolant = mir::utilities::interpolateVertexPosition(p0, p1, t);
  EXPECT_NEAR(  interpolant.m_x, 0.91, 0.00001);
  EXPECT_NEAR(  interpolant.m_y, 0.66, 0.00001);
}

//----------------------------------------------------------------------

TEST(mir_distance_utility, compute_distance)
{
  mir::Point2 p0( 0.25, 0.0 );
  mir::Point2 p1( 1.25, 0.0 );
  axom::float64 dist = mir::utilities::distance(p0, p1);
  EXPECT_DOUBLE_EQ( dist, 1.0 ); 

  mir::Point2 p2( 0.25, 0.0 );
  mir::Point2 p3( 1.25, 1.0 );
  dist = mir::utilities::distance(p2, p3);
  EXPECT_NEAR( dist, 1.4142, 0.001 ); 

  mir::Point2 p4( 1.0, 1.0 );
  mir::Point2 p5( 1.0, 1.0 );
  dist = mir::utilities::distance(p4, p5);
  EXPECT_DOUBLE_EQ( dist, 0.0 ); 

}

//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::UnitTestLogger logger;  // create & initialize test logger,

  result = RUN_ALL_TESTS();
  return result;
}


#endif //  MIR_UTILITIES_TEST_H_

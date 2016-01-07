/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

/*!
 *******************************************************************************
 * \file
 *
 * \date Jan 5, 2016
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */

#include "gtest/gtest.h"

#include "quest/Intersection.hpp"
#include "quest/Point.hpp"
#include "quest/Ray.hpp"
#include "quest/Segment.hpp"
#include "quest/Vector.hpp"


TEST( quest_intersection, ray_segment_intersection )
{
  typedef quest::Point< double,2 >   PointType;
  typedef quest::Segment< double,2 > SegmentType;
  typedef quest::Vector< double,2 >  VectorType;
  typedef quest::Ray< double,2 >     RayType;

  // STEP 0: construct segment
  PointType A(0.0);
  PointType B(1.0,1);
  SegmentType S( A, B );

  // STEP 1: construct ray
  PointType origin = PointType::make_point( 0.5,-0.5 );
  VectorType direction;
  direction[0] = 0.0;
  direction[1] = 0.5;
  RayType R( origin,direction.unitVector() );

  // STEP 2: compute intersection
  PointType ip;
  bool intersects = quest::intersect( R, S, ip );
  EXPECT_TRUE( intersects );
  EXPECT_DOUBLE_EQ(0.5,ip[0]);
  EXPECT_DOUBLE_EQ(0.0,ip[1]);

  // STEP 3: construct non-intersecting ray
  origin[1] = 0.5; // shift R up
  RayType R2( origin, direction.unitVector() );
  bool intersects2 = quest::intersect( R2, S, ip );
  EXPECT_FALSE( intersects2 );
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using asctoolkit::slic::UnitTestLogger;

int main(int argc, char * argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}

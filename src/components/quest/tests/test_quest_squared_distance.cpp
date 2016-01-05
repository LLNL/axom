/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


/*
 * $Id$
 */

/*!
 *******************************************************************************
 * \file
 *
 * \date Dec 9, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */

#include "gtest/gtest.h"

#include "quest/Point.hpp"
#include "quest/SquaredDistance.hpp"
#include "quest/Triangle.hpp"
#include "quest/Segment.hpp"


//------------------------------------------------------------------------------
TEST( quest_squared_distance, point_to_point )
{
  quest::Point< double,2 > A = quest::Point< double,2 >::make_point( 0.0, 0.0 );
  quest::Point< double,2 > B = quest::Point< double,2 >::make_point( 1.5, 1.5 );

  double dist = quest::squared_distance( A, B );
  EXPECT_EQ(4.5, dist);

  double dist2 = quest::squared_distance( B, A );
  EXPECT_EQ( 4.5, dist2 );
}

//------------------------------------------------------------------------------
TEST( quest_squared_distance, point_to_triangle )
{

  // STEP 0: Setup triangle ABC in 3D
  quest::Point< double,3 > A =
          quest::Point< double,3 >::make_point( 0.0, 0.0, 0.0 );
  quest::Point< double,3 > B =
          quest::Point< double,3 >::make_point( 1.5, 1.5, 0.0 );
  quest::Point< double,3 > C =
          quest::Point< double,3 >::make_point( 2.5, 0.0, 0.0 );
  quest::Triangle< double,3 > tri(A,B,C);

  double dist = 0.0;

  // STEP 1: Check distance of triangle nodes to the triangle
  dist = quest::squared_distance( A, tri );
  EXPECT_DOUBLE_EQ( 0.0f, dist );

  dist = quest::squared_distance( B, tri );
  EXPECT_DOUBLE_EQ( 0.0f, dist );

  dist = quest::squared_distance( C, tri );
  EXPECT_DOUBLE_EQ( 0.0f, dist );

  // STEP 2: Check distance of edge midpoints
  quest::Point< double,3 > mid_of_AB = quest::Point< double,3 >::midpoint(A,B);
  dist = quest::squared_distance( mid_of_AB, tri );
  EXPECT_DOUBLE_EQ( 0.0f, dist );

  quest::Point< double,3 > mid_of_BC = quest::Point< double,3 >::midpoint(B,C);
  dist = quest::squared_distance( mid_of_BC, tri );
  EXPECT_DOUBLE_EQ( 0.0f, dist );

  quest::Point< double,3 > mid_of_AC = quest::Point< double,3 >::midpoint(A,C);
  dist = quest::squared_distance( mid_of_AC, tri );
  EXPECT_DOUBLE_EQ( 0.0f, dist );

  // STEP 3: Check distance of point q, inside ABC and orthogonally
  // projected along z+
  quest::Point< double,3 > q=quest::Point< double,3 >::make_point(1.5,0.5,0.5);
  dist = quest::squared_distance( q, tri );
  EXPECT_DOUBLE_EQ( 0.25f, dist );

  // STEP 4: Check distance of point xA, within the voronoi region of A
  quest::Point< double,3 > xA =
          quest::Point< double,3 >::make_point(-1.0,-1.0,0.0);
  dist = quest::squared_distance( xA, tri );
  EXPECT_DOUBLE_EQ( quest::squared_distance(xA,A), dist );

  // STEP 5: Check distance of point xB, within the voronoi region B
  quest::Point< double,3 > xB =
          quest::Point< double,3 >::make_point( 1.5, 2.0, 0.0 );
  dist = quest::squared_distance( xB, tri );
  EXPECT_DOUBLE_EQ( quest::squared_distance(xB,B), dist );

  // STEP 6: Check distance of point xC within the voronoi region C
  quest::Point< double,3 > xC =
          quest::Point< double,3 >::make_point( 3.0, -1.0, 0.5 );
  dist = quest::squared_distance( xC, tri );
  EXPECT_DOUBLE_EQ( quest::squared_distance(xC,C), dist );

}

//------------------------------------------------------------------------------
TEST( quest_squared_distance, point_to_segment )
{
   quest::Point< double,2 > A(0.0);
   quest::Point< double,2 > B(1.0,1);
   quest::Segment< double,2 > S( A,B );

   quest::Point< double,2 > q0 = quest::Point< double,2 >::make_point( 0.5,10 );
   double dist = quest::squared_distance( q0, S );
   EXPECT_DOUBLE_EQ( 100.0f, dist );

  // STEP 0: check source point
  dist = quest::squared_distance( A, S );
  EXPECT_DOUBLE_EQ( 0.0f, dist );

  // STEP 1: check target point
  dist = quest::squared_distance( B, S );
  EXPECT_DOUBLE_EQ( 0.0f, dist );

  // STEP 2: check midpoint
  quest::Point< double,2 > Q = quest::Point< double,2 >::midpoint( A, B );
  dist = quest::squared_distance( Q, S );
  EXPECT_DOUBLE_EQ( 0.0f, dist );

  // STEP 3: check projection to target point
  quest::Point< double,2 > q1(2.0);
  dist = quest::squared_distance( q1,S );
  EXPECT_DOUBLE_EQ( quest::squared_distance(q1,B), dist );

  // STEP 4: check projection source point
  quest::Point< double,2 > q2(-1.0);
  dist = quest::squared_distance( q2,S );
  EXPECT_DOUBLE_EQ( quest::squared_distance(q2,A), dist );

}


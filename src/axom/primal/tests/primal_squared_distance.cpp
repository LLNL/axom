// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/operators/squared_distance.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/Segment.hpp"

using namespace axom;

//------------------------------------------------------------------------------
TEST(primal_squared_distance, array_to_array)
{
  double A[5];
  double B[5];
  int N;
  double d, expect_d;

  {
    SCOPED_TRACE("1D");
    A[0] = 3.25;
    B[0] = 6.6;
    N = 1;
    expect_d = 11.2225;
    d = primal::squared_distance(A, B, N);
    EXPECT_DOUBLE_EQ(d, expect_d);
  }

  {
    SCOPED_TRACE("2D");
    A[0] = 0;
    A[1] = 0;
    B[0] = 1.5;
    B[1] = 1.5;
    N = 2;
    expect_d = 4.5;
    d = primal::squared_distance(A, B, N);
    EXPECT_DOUBLE_EQ(d, expect_d);
  }

  {
    SCOPED_TRACE("3D");
    A[0] = 6;
    A[1] = 2.3;
    A[2] = -1;
    B[0] = 8;
    B[1] = 1.5;
    B[2] = 0;
    N = 3;
    expect_d = 5.64;
    d = primal::squared_distance(A, B, N);
    EXPECT_DOUBLE_EQ(d, expect_d);
  }
}

//------------------------------------------------------------------------------
TEST(primal_squared_distance, point_to_point)
{
  primal::Point<double, 2> A = primal::Point<double, 2>::make_point(0.0, 0.0);
  primal::Point<double, 2> B = primal::Point<double, 2>::make_point(1.5, 1.5);

  double dist = primal::squared_distance(A, B);
  EXPECT_EQ(4.5, dist);

  double dist2 = primal::squared_distance(B, A);
  EXPECT_EQ(4.5, dist2);
}

//------------------------------------------------------------------------------
TEST(primal_squared_distance, point_to_triangle)
{
  // STEP 0: Setup triangle ABC in 3D
  primal::Point<double, 3> A =
    primal::Point<double, 3>::make_point(0.0, 0.0, 0.0);
  primal::Point<double, 3> B =
    primal::Point<double, 3>::make_point(1.5, 1.5, 0.0);
  primal::Point<double, 3> C =
    primal::Point<double, 3>::make_point(2.5, 0.0, 0.0);
  primal::Triangle<double, 3> tri(A, B, C);

  double dist = 0.0;

  // STEP 1: Check distance of triangle nodes to the triangle
  dist = primal::squared_distance(A, tri);
  EXPECT_DOUBLE_EQ(0.0f, dist);

  dist = primal::squared_distance(B, tri);
  EXPECT_DOUBLE_EQ(0.0f, dist);

  dist = primal::squared_distance(C, tri);
  EXPECT_DOUBLE_EQ(0.0f, dist);

  // STEP 2: Check distance of edge midpoints
  primal::Point<double, 3> mid_of_AB = primal::Point<double, 3>::midpoint(A, B);
  dist = primal::squared_distance(mid_of_AB, tri);
  EXPECT_DOUBLE_EQ(0.0f, dist);

  primal::Point<double, 3> mid_of_BC = primal::Point<double, 3>::midpoint(B, C);
  dist = primal::squared_distance(mid_of_BC, tri);
  EXPECT_DOUBLE_EQ(0.0f, dist);

  primal::Point<double, 3> mid_of_AC = primal::Point<double, 3>::midpoint(A, C);
  dist = primal::squared_distance(mid_of_AC, tri);
  EXPECT_DOUBLE_EQ(0.0f, dist);

  // STEP 3: Check distance of point q, inside ABC and orthogonally
  // projected along z+
  primal::Point<double, 3> q =
    primal::Point<double, 3>::make_point(1.5, 0.5, 0.5);
  dist = primal::squared_distance(q, tri);
  EXPECT_DOUBLE_EQ(0.25f, dist);

  // STEP 4: Check distance of point xA, within the voronoi region of A
  primal::Point<double, 3> xA =
    primal::Point<double, 3>::make_point(-1.0, -1.0, 0.0);
  dist = primal::squared_distance(xA, tri);
  EXPECT_DOUBLE_EQ(primal::squared_distance(xA, A), dist);

  // STEP 5: Check distance of point xB, within the voronoi region B
  primal::Point<double, 3> xB =
    primal::Point<double, 3>::make_point(1.5, 2.0, 0.0);
  dist = primal::squared_distance(xB, tri);
  EXPECT_DOUBLE_EQ(primal::squared_distance(xB, B), dist);

  // STEP 6: Check distance of point xC within the voronoi region C
  primal::Point<double, 3> xC =
    primal::Point<double, 3>::make_point(3.0, -1.0, 0.5);
  dist = primal::squared_distance(xC, tri);
  EXPECT_DOUBLE_EQ(primal::squared_distance(xC, C), dist);
}

//------------------------------------------------------------------------------
TEST(primal_squared_distance, point_to_segment)
{
  primal::Point<double, 2> A(0.0);
  primal::Point<double, 2> B(1.0, 1);
  primal::Segment<double, 2> S(A, B);

  primal::Point<double, 2> q0 = primal::Point<double, 2>::make_point(0.5, 10);
  double dist = primal::squared_distance(q0, S);
  EXPECT_DOUBLE_EQ(100.0f, dist);

  // STEP 0: check source point
  dist = primal::squared_distance(A, S);
  EXPECT_DOUBLE_EQ(0.0f, dist);

  // STEP 1: check target point
  dist = primal::squared_distance(B, S);
  EXPECT_DOUBLE_EQ(0.0f, dist);

  // STEP 2: check midpoint
  primal::Point<double, 2> Q = primal::Point<double, 2>::midpoint(A, B);
  dist = primal::squared_distance(Q, S);
  EXPECT_DOUBLE_EQ(0.0f, dist);

  // STEP 3: check projection to target point
  primal::Point<double, 2> q1(2.0);
  dist = primal::squared_distance(q1, S);
  EXPECT_DOUBLE_EQ(primal::squared_distance(q1, B), dist);

  // STEP 4: check projection source point
  primal::Point<double, 2> q2(-1.0);
  dist = primal::squared_distance(q2, S);
  EXPECT_DOUBLE_EQ(primal::squared_distance(q2, A), dist);
}

// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/Triangle.hpp"

#include "axom/primal/operators/orientation.hpp"

using namespace axom;

TEST(primal_orientation, orient3D)
{
  // STEP 0: Setup triangle ABC in 3D
  primal::Point<double, 3> A =
    primal::Point<double, 3>::make_point(0.0, 0.0, 0.0);
  primal::Point<double, 3> B =
    primal::Point<double, 3>::make_point(1.5, 1.5, 0.0);
  primal::Point<double, 3> C =
    primal::Point<double, 3>::make_point(2.5, 0.0, 0.0);
  primal::Triangle<double, 3> tri(A, B, C);

  // STEP 1: Setup test points, q0, q1, q2 => boundary, positive, negative
  primal::Point<double, 3> q0 =
    primal::Point<double, 3>::make_point(1.5, 0.5, 0.0);
  primal::Point<double, 3> q1 =
    primal::Point<double, 3>::make_point(1.5, 0.5, 0.5);
  primal::Point<double, 3> q2 =
    primal::Point<double, 3>::make_point(1.5, 0.5, -0.5);

  // STEP 2: test orientation
  int orient = primal::orientation(q0, tri);
  EXPECT_EQ(primal::ON_BOUNDARY, orient);

  orient = primal::orientation(q1, tri);
  EXPECT_EQ(primal::ON_NEGATIVE_SIDE, orient);

  orient = primal::orientation(q2, tri);
  EXPECT_EQ(primal::ON_POSITIVE_SIDE, orient);
}

//------------------------------------------------------------------------------
TEST(primal_orientation, orient2D)
{
  // STEP 0: create test segment
  primal::Point<double, 2> A(0.0);
  primal::Point<double, 2> B(1.0);
  primal::Segment<double, 2> S(A, B);

  // STEP 1: setup test points, q0, q1, q2 => boundary, positive, negative
  primal::Point<double, 2> q0(0.5);
  primal::Point<double, 2> q1 = primal::Point<double, 2>::make_point(-0.5, 0.5);
  primal::Point<double, 2> q2 = primal::Point<double, 2>::make_point(2.0, 0.5);

  // STEP 2: test orientation
  int orient = primal::orientation(q0, S);
  EXPECT_EQ(primal::ON_BOUNDARY, orient);

  orient = primal::orientation(q1, S);
  EXPECT_EQ(primal::ON_NEGATIVE_SIDE, orient);

  orient = primal::orientation(q2, S);
  EXPECT_EQ(primal::ON_POSITIVE_SIDE, orient);
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
#include "axom/slic/core/SimpleLogger.hpp"
using axom::slic::SimpleLogger;

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  SimpleLogger logger;  // create & initialize test logger,

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}

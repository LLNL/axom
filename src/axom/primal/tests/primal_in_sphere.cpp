// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"

#include "axom/slic/interface/slic.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/operators/in_sphere.hpp"

using namespace axom;

TEST(primal_in_sphere, test_in_sphere_2d)
{
  const int DIM = 2;

  // Define the three vertices of a triangle
  typedef primal::Point<double, DIM> PointType;
  PointType p0 = PointType::make_point(0, 0);
  PointType p1 = PointType::make_point(1, 0);
  PointType p2 = PointType::make_point(0, 1);

  // Test some points that are inside the circumcircle
  {
    // inside triangle
    PointType q1 = PointType::make_point(0.1, 0.1);
    EXPECT_TRUE(in_sphere(q1, p0, p1, p2));

    // outside triangle
    PointType q2 = PointType::make_point(0.78, 0.6);
    EXPECT_TRUE(in_sphere(q2, p0, p1, p2));
  }

  // Test some points that are on the circumcircle
  {
    PointType q1 = PointType::make_point(1, 1);
    EXPECT_FALSE(in_sphere(q1, p0, p1, p2));

    PointType q2 = PointType::make_point(0, 1);
    EXPECT_FALSE(in_sphere(q2, p0, p1, p2));
  }

  // Test some points that are outside the circumcircle
  {
    PointType q1 = PointType::make_point(1.1, 0);
    EXPECT_FALSE(in_sphere(q1, p0, p1, p2));

    PointType q2 = PointType::make_point(-5, -10);
    EXPECT_FALSE(in_sphere(q2, p0, p1, p2));
  }
}

TEST(primal_in_sphere, test_in_sphere_3d)
{
  const int DIM = 3;

  // Define the four vertices of a tetrahedron
  typedef primal::Point<double, DIM> PointType;
  PointType p0 = PointType::make_point(-1, -1, 1);
  PointType p1 = PointType::make_point(1, -1, -1);
  PointType p2 = PointType::make_point(-1, 1, -1);
  PointType p3 = PointType::make_point(1, 1, 1);

  // Test some points that are inside the circumsphere
  {
    PointType q1 = PointType::make_point(0.5, 0.5, 0.5);
    EXPECT_TRUE(in_sphere(q1, p0, p1, p2, p3));

    PointType q2 = PointType::make_point(0., 0., 0.);
    EXPECT_TRUE(in_sphere(q2, p0, p1, p2, p3));
  }

  // Test some points that are on the circumsphere
  {
    PointType q1 = PointType::make_point(1, 1, 1);
    EXPECT_FALSE(in_sphere(q1, p0, p1, p2, p3));

    PointType q2 = PointType::make_point(-1, 1, 1);
    EXPECT_FALSE(in_sphere(q2, p0, p1, p2, p3));
  }

  // Test some points that are outside the circumsphere
  {
    PointType q1 = PointType::make_point(1.1, 1, 1);
    EXPECT_FALSE(in_sphere(q1, p0, p1, p2, p3));

    PointType q2 = PointType::make_point(-1.1, 1, 1);
    EXPECT_FALSE(in_sphere(q2, p0, p1, p2, p3));
  }
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
#include "axom/slic/core/SimpleLogger.hpp"
using axom::slic::SimpleLogger;

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);

  SimpleLogger logger;  // create & initialize test logger,
  axom::slic::setLoggingMsgLevel(axom::slic::message::Warning);

  int result = RUN_ALL_TESTS();
  return result;
}

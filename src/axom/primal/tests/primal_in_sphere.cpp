// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/slic.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/Tetrahedron.hpp"
#include "axom/primal/operators/in_sphere.hpp"

namespace primal = axom::primal;

TEST(primal_in_sphere, test_in_sphere_2d)
{
  const int DIM = 2;
  using PointType = primal::Point<double, DIM>;
  using TriangleType = primal::Triangle<double, DIM>;

  // Define the three vertices of a triangle
  PointType p0 {0, 0};
  PointType p1 {1, 0};
  PointType p2 {0, 1};
  TriangleType tri(p0, p1, p2);

  // Test some points that are inside the circumcircle
  {
    // inside triangle
    PointType q1 {0.1, 0.1};
    EXPECT_TRUE(in_sphere(q1, p0, p1, p2));
    EXPECT_TRUE(in_sphere(q1, tri));

    // outside triangle
    PointType q2 {0.78, 0.6};
    EXPECT_TRUE(in_sphere(q2, p0, p1, p2));
    EXPECT_TRUE(in_sphere(q2, tri));
  }

  // Test some points that are on the circumcircle
  {
    PointType q1 {1, 1};
    EXPECT_FALSE(in_sphere(q1, p0, p1, p2));
    EXPECT_FALSE(in_sphere(q1, tri));

    PointType q2 {0, 1};
    EXPECT_FALSE(in_sphere(q2, p0, p1, p2));
    EXPECT_FALSE(in_sphere(q2, tri));
  }

  // Test some points that are outside the circumcircle
  {
    PointType q1 {1.1, 0};
    EXPECT_FALSE(in_sphere(q1, p0, p1, p2));
    EXPECT_FALSE(in_sphere(q1, tri));

    PointType q2 {-5, -10};
    EXPECT_FALSE(in_sphere(q2, p0, p1, p2));
    EXPECT_FALSE(in_sphere(q2, tri));
  }
}

TEST(primal_in_sphere, test_in_sphere_3d)
{
  const int DIM = 3;
  using PointType = primal::Point<double, DIM>;
  using TetrahedronType = primal::Tetrahedron<double, DIM>;

  // Define the four vertices of a tetrahedron
  PointType p0 {-1, -1, 1};
  PointType p1 {1, -1, -1};
  PointType p2 {-1, 1, -1};
  PointType p3 {1, 1, 1};
  TetrahedronType tet(p0, p1, p2, p3);

  // Test some points that are inside the circumsphere
  {
    PointType q1 {0.5, 0.5, 0.5};
    EXPECT_TRUE(in_sphere(q1, p0, p1, p2, p3));
    EXPECT_TRUE(in_sphere(q1, tet));

    PointType q2 {0., 0., 0.};
    EXPECT_TRUE(in_sphere(q2, p0, p1, p2, p3));
    EXPECT_TRUE(in_sphere(q2, tet));
  }

  // Test some points that are on the circumsphere
  {
    PointType q1 {1, 1, 1};
    EXPECT_FALSE(in_sphere(q1, p0, p1, p2, p3));
    EXPECT_FALSE(in_sphere(q1, tet));

    PointType q2 {-1, 1, 1};
    EXPECT_FALSE(in_sphere(q2, p0, p1, p2, p3));
    EXPECT_FALSE(in_sphere(q2, tet));
  }

  // Test some points that are outside the circumsphere
  {
    PointType q1 {1.1, 1, 1};
    EXPECT_FALSE(in_sphere(q1, p0, p1, p2, p3));
    EXPECT_FALSE(in_sphere(q1, tet));

    PointType q2 {-1.1, 1, 1};
    EXPECT_FALSE(in_sphere(q2, p0, p1, p2, p3));
    EXPECT_FALSE(in_sphere(q2, tet));
  }
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;
  axom::slic::setLoggingMsgLevel(axom::slic::message::Warning);

  int result = RUN_ALL_TESTS();
  return result;
}

// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/primal.hpp"
#include "axom/slic.hpp"

//------------------------------------------------------------------------------
TEST(primal_polygon, polygon_empty)
{
  using PolygonType = axom::primal::Polygon<double, 3>;
  PolygonType poly;
  EXPECT_FALSE(poly.isValid());
}

//------------------------------------------------------------------------------
TEST(primal_polygon, polygon_containment)
{
  using PolygonType = axom::primal::Polygon<double, 2>;
  using PointType = axom::primal::Point<double, 2>;

  axom::Array<PointType> vertices(
    {PointType {0, 0}, PointType {1, 1}, PointType {1, 0}, PointType {0, 1}});
  PolygonType poly(vertices);

  // Simple cases
  EXPECT_TRUE(in_polygon(PointType({0.25, 0.5}), poly));
  EXPECT_TRUE(in_polygon(PointType({0.75, 0.5}), poly));
  EXPECT_FALSE(in_polygon(PointType({0.5, 0.25}), poly));
  EXPECT_FALSE(in_polygon(PointType({0.5, 0.75}), poly));

  // Edge cases. Assume points on the boundary are "inside"
  EXPECT_TRUE(in_polygon(PointType({0, 0.5}), poly));
  EXPECT_TRUE(in_polygon(PointType({0.5, 0.5}), poly));
  EXPECT_TRUE(in_polygon(PointType({0, 0}), poly));

  // Test strict inclusion
  EXPECT_FALSE(in_polygon(PointType({0, 0.5}), poly, true));
  EXPECT_FALSE(in_polygon(PointType({0.5, 0.5}), poly, true));
  EXPECT_FALSE(in_polygon(PointType({0, 0}), poly, true));

  // Corner cases, where edge is aligned with casted ray
  vertices = axom::Array<PointType>({PointType {0, 0},
                                     PointType {0.5, 0},
                                     PointType {1, 0},
                                     PointType {1, 1},
                                     PointType {0, 1}});
  poly = PolygonType(vertices);

  EXPECT_TRUE(in_polygon(PointType({0, 0}), poly));
  EXPECT_TRUE(in_polygon(PointType({0.5, 0}), poly));
  EXPECT_TRUE(in_polygon(PointType({0, 1}), poly));
  EXPECT_TRUE(in_polygon(PointType({1, 0}), poly));
  EXPECT_TRUE(in_polygon(PointType({0, 1}), poly));

  // Test strict inclusion
  EXPECT_FALSE(in_polygon(PointType({0, 0}), poly, true));
  EXPECT_FALSE(in_polygon(PointType({0.5, 0}), poly, true));
  EXPECT_FALSE(in_polygon(PointType({0, 1}), poly, true));
  EXPECT_FALSE(in_polygon(PointType({1, 0}), poly, true));
  EXPECT_FALSE(in_polygon(PointType({0, 1}), poly, true));

  // Verify invariance to orientation
  vertices = axom::Array<PointType>(
    {PointType {0, 1}, PointType {1, 0}, PointType {1, 1}, PointType {0, 0}});
  poly = PolygonType(vertices);

  EXPECT_TRUE(in_polygon(PointType({0.25, 0.5}), poly));
  EXPECT_TRUE(in_polygon(PointType({0.75, 0.5}), poly));
  EXPECT_FALSE(in_polygon(PointType({0.5, 0.25}), poly));
  EXPECT_FALSE(in_polygon(PointType({0.5, 0.75}), poly));
  EXPECT_TRUE(in_polygon(PointType({0, 0.5}), poly));
  EXPECT_TRUE(in_polygon(PointType({0.5, 0.5}), poly));
  EXPECT_TRUE(in_polygon(PointType({0, 0}), poly));

  // Verify invariance to duplicated vertices/degenerate shapes
  vertices = axom::Array<PointType>({PointType {0, 0},
                                     PointType {0, 0},
                                     PointType {1, 1},
                                     PointType {1, 1},
                                     PointType {1, 0},
                                     PointType {1, 0},
                                     PointType {0, 1},
                                     PointType {0, 1}});
  poly = PolygonType(vertices);

  EXPECT_TRUE(in_polygon(PointType({0.25, 0.5}), poly));
  EXPECT_TRUE(in_polygon(PointType({0.75, 0.5}), poly));
  EXPECT_FALSE(in_polygon(PointType({0.5, 0.25}), poly));
  EXPECT_FALSE(in_polygon(PointType({0.5, 0.75}), poly));
  EXPECT_TRUE(in_polygon(PointType({0, 0.5}), poly));
  EXPECT_TRUE(in_polygon(PointType({0.5, 0.5}), poly));
  EXPECT_TRUE(in_polygon(PointType({0, 0}), poly));
}

//------------------------------------------------------------------------------
TEST(primal_polygon, polygon_convexity)
{
  using PolygonType = axom::primal::Polygon<double, 2>;
  using PointType = axom::primal::Point<double, 2>;

  axom::Array<PointType> vertices({PointType {0, 0}, PointType {1, 1}});
  PolygonType poly(vertices);

  // Segments and Triangles are always convex
  EXPECT_TRUE(is_convex(poly));

  poly.addVertex(PointType {0, 1});
  EXPECT_TRUE(is_convex(poly));

  // Duplicate points should not affect convexity
  vertices = axom::Array<PointType>(
    {PointType {0, 0}, PointType {0, 1}, PointType {1, 1}, PointType {1, 0}});
  poly.clear();
  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 3; j++)  // Duplicate each element 3 times
      poly.addVertex(vertices[i]);

  EXPECT_TRUE(is_convex(poly));

  // Verify checks up to rotation of vertices
  vertices = axom::Array<PointType>(
    {PointType {0, 0}, PointType {1, 1}, PointType {1, 0}, PointType {0, 1}});

  for(int i = 0; i < 4; i++)
  {
    poly.clear();
    for(int j = 0; j < 4; j++) poly.addVertex(vertices[(j + i) % 4]);
    EXPECT_FALSE(is_convex(poly));
  }

  vertices = axom::Array<PointType>(
    {PointType {0, 0}, PointType {0, 1}, PointType {1, 1}, PointType {1, 0}});

  for(int i = 0; i < 4; i++)
  {
    poly.clear();
    for(int j = 0; j < 4; j++) poly.addVertex(vertices[(j + i) % 4]);
    EXPECT_TRUE(is_convex(poly));
  }
}

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger(axom::slic::message::Info);

  result = RUN_ALL_TESTS();

  return result;
}

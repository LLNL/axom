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
TEST(primal_polygon, polygon_winding_number)
{
  using PolygonType = axom::primal::Polygon<double, 2>;
  using PointType = axom::primal::Point<double, 2>;

  axom::Array<PointType> vertices(
    {PointType {0, 0}, PointType {1, 1}, PointType {1, 0}, PointType {0, 1}});
  PolygonType poly(vertices);

  // Test the specific winding numbers
  EXPECT_EQ(winding_number(PointType({0.25, 0.5}), poly), 1);
  EXPECT_EQ(winding_number(PointType({0.75, 0.5}), poly), -1);
  EXPECT_EQ(winding_number(PointType({0.5, 0.25}), poly), 0);
  EXPECT_EQ(winding_number(PointType({0.5, 0.75}), poly), 0);

  vertices = axom::Array<PointType>({PointType {0, 1},
                                     PointType {0, -1},
                                     PointType {-2, 2},
                                     PointType {2, 2},
                                     PointType {2, -2},
                                     PointType {-2, -2}});

  poly = PolygonType(vertices);
  EXPECT_EQ(winding_number(PointType({-0.1, 0.0}), poly), -2);
  EXPECT_EQ(winding_number(PointType({0.1, 0.0}), poly), -1);
  EXPECT_EQ(winding_number(PointType({-2.0, 0.0}), poly), 0);
  EXPECT_EQ(winding_number(PointType({2.5, 0.0}), poly), 0);

  // Current policy is to return 1 on edges without strict inclusion,
  //  0 on edges with strict inclusion
  EXPECT_EQ(winding_number(PointType({0.0, 0.0}), poly, false), 1);
  EXPECT_EQ(winding_number(PointType({0.0, 0.0}), poly, true), 0);
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

  // Edge cases, where vertex is aligned with edge
  EXPECT_FALSE(in_polygon(PointType({-0.25, -0.25}), poly));
  EXPECT_FALSE(in_polygon(PointType({1.25, 1.25}), poly));
  EXPECT_FALSE(in_polygon(PointType({-0.25, 1.25}), poly));
  EXPECT_FALSE(in_polygon(PointType({1.25, -0.25}), poly));
  EXPECT_FALSE(in_polygon(PointType({0.0, 1.25}), poly));
  EXPECT_FALSE(in_polygon(PointType({0.0, -0.25}), poly));
  EXPECT_FALSE(in_polygon(PointType({1.0, 1.25}), poly));
  EXPECT_FALSE(in_polygon(PointType({1.0, -0.25}), poly));

  // Test different in/out protocols
  vertices = axom::Array<PointType>({PointType {0, 1},
                                     PointType {0, -1},
                                     PointType {-2, 2},
                                     PointType {2, 2},
                                     PointType {2, -2},
                                     PointType {-2, -2}});

  poly = PolygonType(vertices);

  // true denotes nonzero protocol. Default for SVG
  EXPECT_TRUE(in_polygon(PointType({-0.1, 0.0}), poly));
  EXPECT_TRUE(in_polygon(PointType({-0.1, 0.0}), poly, true));

  // false denotes evenodd protocol
  EXPECT_FALSE(in_polygon(PointType({-0.1, 0.0}), poly, false));
}

TEST(primal_polygon, polygon_containment_invariants)
{
  using PolygonType = axom::primal::Polygon<double, 2>;
  using PointType = axom::primal::Point<double, 2>;

  // Invariant to duplicate points
  axom::Array<PointType> vertices(
    {PointType {0, 0}, PointType {1, 1}, PointType {1, 0}, PointType {0, 1}});
  PolygonType poly;
  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 3; j++)  // Duplicate each element 3 times
      poly.addVertex(vertices[i]);

  EXPECT_TRUE(in_polygon(PointType({0.25, 0.5}), poly));
  EXPECT_TRUE(in_polygon(PointType({0.75, 0.5}), poly));
  EXPECT_FALSE(in_polygon(PointType({0.5, 0.25}), poly));
  EXPECT_FALSE(in_polygon(PointType({0.5, 0.75}), poly));

  // Verify checks up to rotation of vertices
  vertices = axom::Array<PointType>(
    {PointType {0, 0}, PointType {1, 1}, PointType {1, 0}, PointType {0, 1}});

  for(int i = 0; i < 4; i++)
  {
    poly.clear();
    for(int j = 0; j < 4; j++) poly.addVertex(vertices[(j + i) % 4]);
    EXPECT_TRUE(in_polygon(PointType({0.25, 0.5}), poly));
    EXPECT_TRUE(in_polygon(PointType({0.75, 0.5}), poly));
    EXPECT_FALSE(in_polygon(PointType({0.5, 0.25}), poly));
    EXPECT_FALSE(in_polygon(PointType({0.5, 0.75}), poly));
  }
}

TEST(primal_polygon, polygon_containment_edge)
{
  using PolygonType = axom::primal::Polygon<double, 2>;
  using PointType = axom::primal::Point<double, 2>;

  axom::Array<PointType> vertices({PointType {0, 0},
                                   PointType {0.5, 0.5},
                                   PointType {1, 1},
                                   PointType {1, 0},
                                   PointType {0, 1}});
  PolygonType poly(vertices);

  // Edge cases. Default is that points on edges are "inside".
  //  Should work in either in/out protocol
  bool strict = false;
  for(bool nonzero : {true, false})
  {
    EXPECT_TRUE(in_polygon(PointType({0, 0.5}), poly, nonzero, strict));
    EXPECT_TRUE(in_polygon(PointType({0.5, 0.5}), poly, nonzero, strict));
    EXPECT_TRUE(in_polygon(PointType({0.25, 0.25}), poly, nonzero, strict));
    EXPECT_TRUE(in_polygon(PointType({0.25, 0.75}), poly, nonzero, strict));
    EXPECT_TRUE(in_polygon(PointType({1, 0.5}), poly, nonzero, strict));
  }

  strict = true;
  for(bool nonzero : {true, false})
  {
    EXPECT_FALSE(in_polygon(PointType({0, 0.5}), poly, nonzero, strict));
    EXPECT_FALSE(in_polygon(PointType({0.5, 0.5}), poly, nonzero, strict));
    EXPECT_FALSE(in_polygon(PointType({0.25, 0.25}), poly, nonzero, strict));
    EXPECT_FALSE(in_polygon(PointType({0.25, 0.75}), poly, nonzero, strict));
    EXPECT_FALSE(in_polygon(PointType({1, 0.5}), poly, nonzero, strict));
  }

  // Corner cases, where query is on a vertex
  vertices = axom::Array<PointType>({PointType {0, 0},
                                     PointType {0.5, 0},
                                     PointType {1, 0},
                                     PointType {1, 1},
                                     PointType {0, 1}});
  poly = PolygonType(vertices);

  for(auto& vtx : vertices)
  {
    // Nonzero in/out protocol
    EXPECT_TRUE(in_polygon(vtx, poly, true, false));
    EXPECT_FALSE(in_polygon(vtx, poly, true, true));

    // Evenodd in/out protocol
    EXPECT_TRUE(in_polygon(vtx, poly, false, false));
    EXPECT_FALSE(in_polygon(vtx, poly, false, true));
  }
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

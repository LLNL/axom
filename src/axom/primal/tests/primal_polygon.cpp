// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/primal.hpp"
#include "axom/slic.hpp"

//------------------------------------------------------------------------------
TEST(primal_polygon, empty)
{
  using PolygonType = axom::primal::Polygon<double, 3>;
  PolygonType poly;
  EXPECT_FALSE(poly.isValid());
}

//------------------------------------------------------------------------------
TEST(primal_polygon, winding_number)
{
  using PolygonType = axom::primal::Polygon<double, 2>;
  using PointType = axom::primal::Point<double, 2>;

  const bool useStrictInclusion = true;

  axom::Array<PointType> vertices(
    {PointType {0, 0}, PointType {1, 1}, PointType {1, 0}, PointType {0, 1}});
  PolygonType poly(vertices);

  // Test the specific winding numbers
  EXPECT_EQ(winding_number(PointType {0.25, 0.5}, poly), 1);
  EXPECT_EQ(winding_number(PointType {0.75, 0.5}, poly), -1);
  EXPECT_EQ(winding_number(PointType {0.5, 0.25}, poly), 0);
  EXPECT_EQ(winding_number(PointType {0.5, 0.75}, poly), 0);

  vertices = axom::Array<PointType>({PointType {0, 1},
                                     PointType {0, -1},
                                     PointType {-2, 2},
                                     PointType {2, 2},
                                     PointType {2, -2},
                                     PointType {-2, -2}});

  poly = PolygonType(vertices);
  EXPECT_EQ(winding_number(PointType {-0.1, 0.0}, poly), -2);
  EXPECT_EQ(winding_number(PointType {0.1, 0.0}, poly), -1);
  EXPECT_EQ(winding_number(PointType {-2.0, 0.0}, poly), 0);
  EXPECT_EQ(winding_number(PointType {2.5, 0.0}, poly), 0);

  // Current policy is to return 1 on edges without strict inclusion,
  //  0 on edges with strict inclusion, as 0 always indicates "interior"
  EXPECT_EQ(winding_number(PointType {0.0, 0.0}, poly, !useStrictInclusion), 1);
  EXPECT_EQ(winding_number(PointType {0.0, 0.0}, poly, useStrictInclusion), 0);
}

//------------------------------------------------------------------------------
TEST(primal_polygon, containment)
{
  using PolygonType = axom::primal::Polygon<double, 2>;
  using PointType = axom::primal::Point<double, 2>;

  const bool useNonzeroRule = true;
  const bool useEvenOddRule = false;

  axom::Array<PointType> vertices(
    {PointType {0, 0}, PointType {1, 1}, PointType {1, 0}, PointType {0, 1}});
  PolygonType poly(vertices);

  // Simple cases on non-convex polygon
  EXPECT_TRUE(in_polygon(PointType {0.25, 0.5}, poly));
  EXPECT_TRUE(in_polygon(PointType {0.75, 0.5}, poly));
  EXPECT_FALSE(in_polygon(PointType {0.5, 0.25}, poly));
  EXPECT_FALSE(in_polygon(PointType {0.5, 0.75}, poly));

  // Edge cases, where vertex is aligned with edge
  EXPECT_FALSE(in_polygon(PointType {-0.25, -0.25}, poly));
  EXPECT_FALSE(in_polygon(PointType {1.25, 1.25}, poly));
  EXPECT_FALSE(in_polygon(PointType {-0.25, 1.25}, poly));
  EXPECT_FALSE(in_polygon(PointType {1.25, -0.25}, poly));
  EXPECT_FALSE(in_polygon(PointType {0.0, 1.25}, poly));
  EXPECT_FALSE(in_polygon(PointType {0.0, -0.25}, poly));
  EXPECT_FALSE(in_polygon(PointType {-1.0, 0.0}, poly));
  EXPECT_FALSE(in_polygon(PointType {1.0, 1.25}, poly));
  EXPECT_FALSE(in_polygon(PointType {1.0, -0.25}, poly));

  // Test different in/out protocols for polygon with extra loop
  vertices = axom::Array<PointType>({PointType {0, 1},
                                     PointType {0, -1},
                                     PointType {-2, 2},
                                     PointType {2, 2},
                                     PointType {2, -2},
                                     PointType {-2, -2}});

  poly = PolygonType(vertices);

  // true denotes nonzero protocol. Default for SVG
  EXPECT_TRUE(in_polygon(PointType {-0.1, 0.0}, poly));
  EXPECT_TRUE(in_polygon(PointType {-0.1, 0.0}, poly, useNonzeroRule));

  // false denotes evenodd protocol
  EXPECT_FALSE(in_polygon(PointType {-0.1, 0.0}, poly, useEvenOddRule));

  // Test in/out on degenerate example in Hormann2001, Figure 6
  axom::Array<PointType> init_vertices(
    {PointType {1, 1}, PointType {1, -2}, PointType {10, -2}, PointType {10, -1}});

  poly = PolygonType(init_vertices);

  axom::Array<PointType> degenerate_vertices({PointType {9, 1},
                                              PointType {8, 0},
                                              PointType {6, 1},
                                              PointType {7, 0},
                                              PointType {5, -1},
                                              PointType {4, 0},
                                              PointType {3, 0}});

  // This point should remain interior when each vertex is added
  EXPECT_TRUE(in_polygon(PointType({2.0, 0.0}), poly));
  for(auto& vertex : degenerate_vertices)
  {
    poly.addVertex(vertex);
    EXPECT_TRUE(in_polygon(PointType({2.0, 0.0}), poly));
  }
}

TEST(primal_polygon, containment_invariants)
{
  using PolygonType = axom::primal::Polygon<double, 2>;
  using PointType = axom::primal::Point<double, 2>;

  // Invariant to duplicate points
  axom::Array<PointType> vertices(
    {PointType {0, 0}, PointType {1, 1}, PointType {1, 0}, PointType {0, 1}});
  PolygonType poly;
  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 3; j++)  // Duplicate each element 3 times
    {
      poly.addVertex(vertices[i]);
    }

  EXPECT_TRUE(in_polygon(PointType {0.25, 0.5}, poly));
  EXPECT_TRUE(in_polygon(PointType {0.75, 0.5}, poly));
  EXPECT_FALSE(in_polygon(PointType {0.5, 0.25}, poly));
  EXPECT_FALSE(in_polygon(PointType {0.5, 0.75}, poly));

  // Verify checks up to rotation of vertices
  vertices = axom::Array<PointType>(
    {PointType {0, 0}, PointType {1, 1}, PointType {1, 0}, PointType {0, 1}});

  for(int i = 0; i < 4; i++)
  {
    poly.clear();
    for(int j = 0; j < 4; j++) poly.addVertex(vertices[(j + i) % 4]);
    EXPECT_TRUE(in_polygon(PointType {0.25, 0.5}, poly));
    EXPECT_TRUE(in_polygon(PointType {0.75, 0.5}, poly));
    EXPECT_FALSE(in_polygon(PointType {0.5, 0.25}, poly));
    EXPECT_FALSE(in_polygon(PointType {0.5, 0.75}, poly));
  }
}

TEST(primal_polygon, containment_edge)
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
  const bool useStrictInclusion = true;
  for(bool useNonzeroRule : {true, false})
  {
    EXPECT_TRUE(
      in_polygon(PointType {0, 0.5}, poly, useNonzeroRule, !useStrictInclusion));
    EXPECT_TRUE(
      in_polygon(PointType {0.5, 0.5}, poly, useNonzeroRule, !useStrictInclusion));
    EXPECT_TRUE(
      in_polygon(PointType {.25, .25}, poly, useNonzeroRule, !useStrictInclusion));
    EXPECT_TRUE(
      in_polygon(PointType {.25, .75}, poly, useNonzeroRule, !useStrictInclusion));
    EXPECT_TRUE(
      in_polygon(PointType {1, 0.5}, poly, useNonzeroRule, !useStrictInclusion));
  }

  for(bool useNonzeroRule : {true, false})
  {
    EXPECT_FALSE(
      in_polygon(PointType {0, 0.5}, poly, useNonzeroRule, useStrictInclusion));
    EXPECT_FALSE(
      in_polygon(PointType {0.5, 0.5}, poly, useNonzeroRule, useStrictInclusion));
    EXPECT_FALSE(
      in_polygon(PointType {.25, .25}, poly, useNonzeroRule, useStrictInclusion));
    EXPECT_FALSE(
      in_polygon(PointType {.25, .75}, poly, useNonzeroRule, useStrictInclusion));
    EXPECT_FALSE(
      in_polygon(PointType {1, 0.5}, poly, useNonzeroRule, useStrictInclusion));
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
    EXPECT_TRUE(in_polygon(vtx, poly, true, !useStrictInclusion));
    EXPECT_FALSE(in_polygon(vtx, poly, true, useStrictInclusion));

    // Evenodd in/out protocol
    EXPECT_TRUE(in_polygon(vtx, poly, false, !useStrictInclusion));
    EXPECT_FALSE(in_polygon(vtx, poly, false, useStrictInclusion));
  }
}

//------------------------------------------------------------------------------
TEST(primal_polygon, convexity)
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
    {
      poly.addVertex(vertices[i]);
    }

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
TEST(primal_polygon, signed_area_2d)
{
  using Polygon2D = axom::primal::Polygon<double, 2>;
  using Point2D = axom::primal::Point<double, 2>;
  using axom::utilities::abs;

  // test a simple right triangle in CW and CCW orientations
  {
    Polygon2D poly2D_ccw({Point2D {0, 0}, Point2D {1, 0}, Point2D {1, 1}});
    Polygon2D poly2D_cw({Point2D {0, 0}, Point2D {1, 1}, Point2D {1, 0}});

    // Signed area is positive for CCW and negative for CW
    EXPECT_DOUBLE_EQ(.5, poly2D_ccw.signedArea());
    EXPECT_DOUBLE_EQ(-.5, poly2D_cw.signedArea());

    // The two triangles have reverse orientations (signedArea)
    // but the same (unsigned) area
    EXPECT_DOUBLE_EQ(-poly2D_ccw.signedArea(), poly2D_cw.signedArea());
    EXPECT_DOUBLE_EQ(poly2D_ccw.area(), poly2D_cw.area());

    // compare signed and unsigned areas
    EXPECT_DOUBLE_EQ(poly2D_ccw.area(), poly2D_ccw.signedArea());
    EXPECT_DOUBLE_EQ(-poly2D_cw.area(), poly2D_cw.signedArea());
  }

  // test regular polygons with CW and CCW orienations
  for(int nSides = 3; nSides < 10; ++nSides)
  {
    Polygon2D poly2D_ccw(nSides);
    Polygon2D poly2D_cw(nSides);

    for(int i = 0; i < nSides; ++i)
    {
      const double angle = 2. * M_PI * i / nSides;
      poly2D_ccw.addVertex(Point2D {cos(angle), sin(angle)});
      poly2D_cw.addVertex(Point2D {sin(angle), cos(angle)});
    }

    const double expected_area = nSides / 2. * sin(2 * M_PI / nSides);

    // The areas are the same; signed areas are opposite
    EXPECT_DOUBLE_EQ(expected_area, poly2D_ccw.area());
    EXPECT_DOUBLE_EQ(expected_area, poly2D_cw.area());
    EXPECT_DOUBLE_EQ(poly2D_cw.area(), poly2D_ccw.area());
    EXPECT_DOUBLE_EQ(-poly2D_cw.signedArea(), poly2D_ccw.signedArea());

    EXPECT_DOUBLE_EQ(poly2D_ccw.signedArea(), poly2D_ccw.area());
    EXPECT_DOUBLE_EQ(-poly2D_cw.signedArea(), poly2D_cw.area());
  }
}

//------------------------------------------------------------------------------
TEST(primal_polygon, area_2d_3d)
{
  using Polygon2D = axom::primal::Polygon<double, 2>;
  using Point2D = axom::primal::Point<double, 2>;

  using Polygon3D = axom::primal::Polygon<double, 3>;
  using Point3D = axom::primal::Point<double, 3>;

  // test a simple right triangle
  // use same xy-data in 2D and 3D
  {
    Polygon2D poly2D({Point2D {0, 0}, Point2D {1, 0}, Point2D {1, 1}});
    EXPECT_DOUBLE_EQ(.5, poly2D.area());

    Polygon3D poly3Da({Point3D {0, 0, 0}, Point3D {1, 0, 0}, Point3D {1, 1, 0}});
    EXPECT_DOUBLE_EQ(.5, poly3Da.area());

    Polygon3D poly3Db({Point3D {0, 0, 1}, Point3D {1, 0, 1}, Point3D {1, 1, 1}});
    EXPECT_DOUBLE_EQ(.5, poly3Db.area());
  }

  // test regular polygons
  // use same xy-data in 2D and 3D
  for(int nSides = 3; nSides < 10; ++nSides)
  {
    Polygon2D poly2D(nSides);
    Polygon3D poly3D(nSides);

    // choose an arbitrary z_offset for 3D polygon
    const double z_offset = 5.;

    for(int i = 0; i < nSides; ++i)
    {
      const double angle = 2. * M_PI * i / nSides;
      poly2D.addVertex(Point2D {cos(angle), sin(angle)});
      poly3D.addVertex(Point3D {cos(angle), sin(angle), z_offset});
    }

    const double expected_area = nSides / 2. * sin(2 * M_PI / nSides);

    EXPECT_DOUBLE_EQ(expected_area, poly2D.area());
    EXPECT_DOUBLE_EQ(expected_area, poly3D.area());
    EXPECT_DOUBLE_EQ(poly2D.area(), poly3D.area());
  }
}

//------------------------------------------------------------------------------
TEST(primal_polygon, normal)
{
  using Polygon3D = axom::primal::Polygon<double, 3>;
  using Point3D = axom::primal::Point<double, 3>;

  using Triangle3D = axom::primal::Triangle<double, 3>;
  using Vector3D = axom::primal::Vector<double, 3>;

  // Test a simple right triangle
  {
    Polygon3D poly3D({Point3D {0, 0, 1}, Point3D {1, 0, 0}, Point3D {1, 1, 1}});
    Vector3D poly_normal = poly3D.normal().unitVector();

    Triangle3D tri3D({Point3D {0, 0, 1}, Point3D {1, 0, 0}, Point3D {1, 1, 1}});
    Vector3D tri_normal = tri3D.normal().unitVector();

    for(int i = 0; i < 3; ++i)
      EXPECT_NEAR(poly_normal[i], tri_normal[i], 1e-10);
  }

  // Test a planar 5-gon inscribed in a circle
  {
    Vector3D v1 = Vector3D({0.0, 1.0, 2.0}).unitVector();
    Vector3D v2 = Vector3D({2.0, -1.0, 0.5}).unitVector();
    Vector3D exp_normal = Vector3D::cross_product(v1, v2).unitVector();

    Polygon3D poly(5);

    double angles[5] = {0.0, 0.5, 1.2, 3.0, 5.0};

    for(int i = 0; i < 5; ++i)
      poly.addVertex(Point3D {cos(angles[i]) * v1[0] + sin(angles[i]) * v2[0],
                              cos(angles[i]) * v1[1] + sin(angles[i]) * v2[1],
                              cos(angles[i]) * v1[2] + sin(angles[i]) * v2[2]});

    Vector3D obs_normal = poly.normal().unitVector();

    for(int i = 0; i < 3; ++i) EXPECT_NEAR(exp_normal[i], obs_normal[i], 1e-10);
  }

  // Test a planar "5-gon" with duplicate vertices and colinear points (robustness)
  {
    Vector3D v1 = Vector3D({0.0, 1.0, 2.0}).unitVector();
    Vector3D v2 = Vector3D({2.0, -1.0, 0.5}).unitVector();
    Vector3D exp_normal = Vector3D::cross_product(v1, v2).unitVector();

    Polygon3D poly(10);

    double angles[8] = {0.0, 0.5, 0.5, 1.2, 3.0, 3.0, 3.0, 5.0};

    // Add point at angle 0
    poly.addVertex(Point3D {cos(angles[0]) * v1[0] + sin(angles[0]) * v2[0],
                            cos(angles[0]) * v1[1] + sin(angles[0]) * v2[1],
                            cos(angles[0]) * v1[2] + sin(angles[0]) * v2[2]});
    // Add a midpoint between angles 0 and 1
    poly.addVertex(Point3D::midpoint(
      poly[0],
      Point3D {cos(angles[1]) * v1[0] + sin(angles[1]) * v2[0],
               cos(angles[1]) * v1[1] + sin(angles[1]) * v2[1],
               cos(angles[1]) * v1[2] + sin(angles[1]) * v2[2]}));

    // Add the rest of the vertices
    for(int i = 1; i < 8; ++i)
      poly.addVertex(Point3D {cos(angles[i]) * v1[0] + sin(angles[i]) * v2[0],
                              cos(angles[i]) * v1[1] + sin(angles[i]) * v2[1],
                              cos(angles[i]) * v1[2] + sin(angles[i]) * v2[2]});

    // Add another midpoint
    poly.addVertex(Point3D::midpoint(poly[0], poly[8]));

    // Verify degeneracies
    EXPECT_NEAR(0.0, Triangle3D(poly[0], poly[1], poly[2]).area(), 1e-10);
    EXPECT_NEAR(0.0, Triangle3D(poly[2], poly[3], poly[4]).area(), 1e-10);

    Vector3D obs_normal = poly.normal().unitVector();

    EXPECT_EQ(10, poly.numVertices());
    for(int i = 0; i < 3; ++i) EXPECT_NEAR(exp_normal[i], obs_normal[i], 1e-10);
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

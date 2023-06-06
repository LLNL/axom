// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/primal.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/fmt.hpp"

#include "gtest/gtest.h"

// C++ headers
#include <cmath>
#include <iostream>
#include <fstream>

namespace primal = axom::primal;

TEST(primal_solid_angle, triangle)
{
  using Point3D = primal::Point<double, 3>;
  using Vector3D = primal::Vector<double, 3>;
  using Triangle = primal::Triangle<double, 3>;

  Point3D origin {0.0, 0.0, 0.0};

  // Test on octant
  Triangle octant(Point3D {1.0, 0.0, 0.0},
                  Point3D {0.0, 1.0, 0.0},
                  Point3D {0.0, 0.0, 1.0});
  EXPECT_NEAR(1.0 / 8.0, winding_number(origin, octant), 1e-10);

  // Test on various points with alternate formula (L'Huilier's)
  Triangle tri(Point3D {2, 4, 3},
               Point3D {0.0, -1.0, 0.0},
               Point3D {-1.0, 0.0, 0.0});

  Point3D queries[10] = {
    Point3D {0.0, 0.0, 0.0},
    Point3D {1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0},
    Point3D {0.2, 0.2, 0.2},
    Point3D {0.4, 0.4, 0.4},
    Point3D {0.5, 0.5, 0.5},
    Point3D {1.0, 1.0, 1.0},
    Point3D {1.0, 2.0, 0.0},
    Point3D {0.0, 5.0, 0.0},
    Point3D {0.0, 0.0, 5.0},
    Point3D {3.0, 2.0, 1.0},
  };

  for(int n = 0; n < 10; ++n)
  {
    // Use direct formula for unsigned solid angle
    double theta_a = Triangle(queries[n], tri[0], tri[1]).angle(0);
    double theta_b = Triangle(queries[n], tri[1], tri[2]).angle(0);
    double theta_c = Triangle(queries[n], tri[2], tri[0]).angle(0);

    double theta_s = 0.5 * (theta_a + theta_b + theta_c);

    double solid_angle = 4 *
      atan(sqrt(tan(0.5 * theta_s) * tan(0.5 * (theta_s - theta_a)) *
                tan(0.5 * (theta_s - theta_b)) * tan(0.5 * (theta_s - theta_c))));

    // Get sign from orientation.
    if(primal::orientation(queries[n], tri) == primal::ON_NEGATIVE_SIDE)
      // Means query point is interior
      EXPECT_NEAR(0.25 * M_1_PI * solid_angle,
                  winding_number(queries[n], tri),
                  1e-10);
    else
      // Means query point is exterior
      EXPECT_NEAR(-0.25 * M_1_PI * solid_angle,
                  winding_number(queries[n], tri),
                  1e-10);
  }

  // Test with simple degenerate triangle
  Triangle deg(Point3D {1.0, 2.0, 3.0},
               Point3D {2.0, 3.0, 4.0},
               Point3D {1.0, 2.0, 3.0});
  for(int n = 0; n < 5; ++n)
    EXPECT_DOUBLE_EQ(0.0, winding_number(queries[n], deg));

  // Test with nondegenerate triangle, zero winding number
  Point3D coplanar {-1.0, 1.0, 1.0};
  EXPECT_NEAR(0.0, winding_number(coplanar, octant), 1e-10);

  // Test with coincident point
  Point3D coincident {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
  EXPECT_NEAR(0.0, winding_number(coincident, octant), 1e-10);

  // Test with point at vertex
  Point3D vertex {0.0, 1.0, 0.0};
  EXPECT_DOUBLE_EQ(0.0, winding_number(vertex, octant), 1e-10);
}

//------------------------------------------------------------------------------
TEST(primal_solid_angle, simple_polygon)
{
  using Point3D = primal::Point<double, 3>;
  using Triangle = primal::Triangle<double, 3>;
  using Polygon = primal::Polygon<double, 3>;

  Point3D origin {0.0, 0.0, 0.0};

  // Test on a triangle
  axom::Array<Point3D> octant_vertices(
    {Point3D {1.0, 0.0, 0.0}, Point3D {0.0, 1.0, 0.0}, Point3D {0.0, 0.0, 1.0}});
  Polygon octant(octant_vertices);

  EXPECT_NEAR(1.0 / 8.0, winding_number(origin, octant), 1e-10);

  // Reverse the orientation
  axom::Array<Point3D> reversed_octant_vertices(
    {Point3D {0.0, 0.0, 1.0}, Point3D {0.0, 1.0, 0.0}, Point3D {1.0, 0.0, 0.0}});
  Polygon reversed_octant(reversed_octant_vertices);

  EXPECT_NEAR(-1.0 / 8.0, winding_number(origin, reversed_octant), 1e-10);

  // Test on a simple square shape, compare with triangularization
  axom::Array<Point3D> square_vertices({Point3D {1.0, 0.0, 0.0},
                                        Point3D {1.0, 1.0, 0.0},
                                        Point3D {0.0, 1.0, 0.0},
                                        Point3D {0.0, 0.0, 0.0}});
  Polygon square(square_vertices);
  Triangle square_tri1(square[0], square[1], square[2]);
  Triangle square_tri2(square[2], square[3], square[0]);

  Point3D queries[5] = {
    Point3D {0.0, 4.0, 1.0},
    Point3D {-1.0, 2.0, 2.0},
    Point3D {0.0, -5.0, 3.0},
    Point3D {0.0, 0.0, 4.0},
    Point3D {3.0, 2.0, 5.0},
  };

  for(int n = 0; n < 5; ++n)
    EXPECT_NEAR(winding_number(queries[n], square),
                winding_number(queries[n], square_tri1) +
                  winding_number(queries[n], square_tri2),
                1e-10);

  // Test on coplanar points
  Point3D troubled_queries[5] = {
    Point3D {0.5, 0.5, 0.0},
    Point3D {1.0, 0.5, 0.0},
    Point3D {0.5, 1.0, 0.0},
    Point3D {1.5, 1.5, 0.0},
    Point3D {0.0, 0.0, 0.0},
  };

  for(int n = 0; n < 5; ++n)
    EXPECT_NEAR(winding_number(troubled_queries[n], square), 0, 1e-10);
}

//------------------------------------------------------------------------------
TEST(primal_solid_angle, nonconvex_polygon)
{
  using Point3D = primal::Point<double, 3>;
  using Triangle = primal::Triangle<double, 3>;
  using Polygon = primal::Polygon<double, 3>;

  Point3D origin {0.0, 0.0, 0.0};

  // Test on a nonconvex pentagon shape, compare with triangularization
  axom::Array<Point3D> pentagon_vertices({Point3D {1.0, 0.0, 0.0},
                                          Point3D {0.5, 0.5, 0.0},
                                          Point3D {1.0, 1.0, 0.0},
                                          Point3D {0.0, 1.0, 0.0},
                                          Point3D {0.0, 0.0, 0.0}});
  Polygon pentagon(pentagon_vertices);
  Triangle pentagon_tri1(pentagon[0], pentagon[1], pentagon[4]);
  Triangle pentagon_tri2(pentagon[2], pentagon[3], pentagon[4]);

  Point3D queries[5] = {
    Point3D {0.0, 4.0, 1.0},
    Point3D {-1.0, 2.0, 2.0},
    Point3D {0.0, -5.0, 3.0},
    Point3D {0.0, 0.0, 4.0},
    Point3D {3.0, 2.0, 5.0},
  };

  for(int n = 0; n < 5; ++n)
    EXPECT_NEAR(winding_number(queries[n], pentagon),
                winding_number(queries[n], pentagon_tri1) +
                  winding_number(queries[n], pentagon_tri2),
                1e-10);

  // Test on coplanar points
  Point3D troubled_queries[5] = {
    Point3D {0.5, 0.5, 0.0},
    Point3D {1.0, 0.5, 0.0},
    Point3D {0.5, 1.0, 0.0},
    Point3D {1.5, 1.5, 0.0},
    Point3D {0.0, 0.0, 0.0},
  };

  for(int n = 0; n < 5; ++n)
    EXPECT_NEAR(winding_number(troubled_queries[n], pentagon), 0, 1e-10);
}

//------------------------------------------------------------------------------
TEST(primal_solid_angle, degenerate_polygon)
{
  using Point3D = primal::Point<double, 3>;
  using Vector3D = primal::Vector<double, 3>;
  using Triangle = primal::Triangle<double, 3>;
  using Polygon = primal::Polygon<double, 3>;

  Vector3D v1 = Vector3D({0.0, 1.0, 2.0}).unitVector();
  Vector3D v2 = Vector3D({2.0, -1.0, 0.5}).unitVector();

  Polygon good_pentagon(5), bad_pentagon(9);
  double good_angles[5] = {0.0, 0.5, 1.2, 3.0, 5.0};

  // Test with final vertex coincident with initial vertex
  double bad_angles[9] = {0.0, 0.5, 0.5, 1.2, 3.0, 3.0, 3.0, 5.0, 0.0};

  // Add vertices to good polygon
  for(int i = 0; i < 5; ++i)
    good_pentagon.addVertex(
      Point3D {cos(good_angles[i]) * v1[0] + sin(good_angles[i]) * v2[0],
               cos(good_angles[i]) * v1[1] + sin(good_angles[i]) * v2[1],
               cos(good_angles[i]) * v1[2] + sin(good_angles[i]) * v2[2]});

  // Create bad polygon with degeneracies
  // Add point at angle 0
  bad_pentagon.addVertex(
    Point3D {cos(bad_angles[0]) * v1[0] + sin(bad_angles[0]) * v2[0],
             cos(bad_angles[0]) * v1[1] + sin(bad_angles[0]) * v2[1],
             cos(bad_angles[0]) * v1[2] + sin(bad_angles[0]) * v2[2]});
  // Add a midpoint between angles 0 and 1
  bad_pentagon.addVertex(Point3D::midpoint(
    bad_pentagon[0],
    Point3D {cos(bad_angles[1]) * v1[0] + sin(bad_angles[1]) * v2[0],
             cos(bad_angles[1]) * v1[1] + sin(bad_angles[1]) * v2[1],
             cos(bad_angles[1]) * v1[2] + sin(bad_angles[1]) * v2[2]}));
  // Add the rest of the vertices
  for(int i = 1; i < 9; ++i)
    bad_pentagon.addVertex(
      Point3D {cos(bad_angles[i]) * v1[0] + sin(bad_angles[i]) * v2[0],
               cos(bad_angles[i]) * v1[1] + sin(bad_angles[i]) * v2[1],
               cos(bad_angles[i]) * v1[2] + sin(bad_angles[i]) * v2[2]});

  Point3D queries[5] = {
    Point3D {0.0, 4.0, 1.0},
    Point3D {-1.0, 2.0, 2.0},
    Point3D {0.0, -5.0, 3.0},
    Point3D {0.0, 0.0, 4.0},
    Point3D {3.0, 2.0, 5.0},
  };

  EXPECT_NEAR(good_pentagon.area(), bad_pentagon.area(), 1e-10);
  for(int n = 0; n < 5; ++n)
    EXPECT_NEAR(winding_number(queries[n], good_pentagon),
                winding_number(queries[n], bad_pentagon),
                1e-10);

  // Test with duplicate point at endpoint
  axom::Array<Point3D> good_square_vertices({Point3D {1.0, 0.0, 0.0},
                                             Point3D {1.0, 1.0, 0.0},
                                             Point3D {0.0, 1.0, 0.0},
                                             Point3D {0.0, 0.0, 0.0}});

  axom::Array<Point3D> bad_square_vertices({Point3D {1.0, 0.0, 0.0},
                                            Point3D {1.0, 1.0, 0.0},
                                            Point3D {0.0, 1.0, 0.0},
                                            Point3D {0.0, 0.0, 0.0},
                                            Point3D {0.0, 0.0, 0.0}});

  for(int n = 0; n < 5; ++n)
    EXPECT_NEAR(winding_number(queries[n], Polygon(good_square_vertices)),
                winding_number(queries[n], Polygon(bad_square_vertices)),
                1e-10);
}

//------------------------------------------------------------------------------
TEST(primal_solid_angle, selfintersecting_star)
{
  using Point3D = primal::Point<double, 3>;
  using Vector3D = primal::Vector<double, 3>;
  using Triangle = primal::Triangle<double, 3>;
  using Polygon = primal::Polygon<double, 3>;

  Point3D query {0, 0, 10};

  Vector3D v1 = Vector3D({0.0, 1.0, 2.0}).unitVector();
  Vector3D v2 = Vector3D({2.0, -1.0, 0.5}).unitVector();

  Polygon pentagram(5);
  double outer_angles[5] = {1 * M_PI / 10,
                            9 * M_PI / 10,
                            17 * M_PI / 10,
                            5 * M_PI / 10,
                            13 * M_PI / 10};

  Polygon pentagon(5);
  double inner_angles[5] {3 * M_PI / 10,
                          7 * M_PI / 10,
                          11 * M_PI / 10,
                          15 * M_PI / 10,
                          19 * M_PI / 10};

  double r0 = sin(M_PI / 10) / sin(7 * M_PI / 10);

  // Add vertices to pentagram
  for(int i = 0; i < 5; ++i)
    pentagram.addVertex(
      Point3D {cos(outer_angles[i]) * v1[0] + sin(outer_angles[i]) * v2[0],
               cos(outer_angles[i]) * v1[1] + sin(outer_angles[i]) * v2[1],
               cos(outer_angles[i]) * v1[2] + sin(outer_angles[i]) * v2[2]});

  // Construct the inner pentagon
  for(int i = 0; i < 5; ++i)
    pentagon.addVertex(Point3D {
      r0 * cos(inner_angles[i]) * v1[0] + r0 * sin(inner_angles[i]) * v2[0],
      r0 * cos(inner_angles[i]) * v1[1] + r0 * sin(inner_angles[i]) * v2[1],
      r0 * cos(inner_angles[i]) * v1[2] + r0 * sin(inner_angles[i]) * v2[2]});

  // Construct the stars of the pentagram
  Polygon star_tips[5];
  star_tips[0] =
    Polygon(axom::Array<Point3D>({pentagram[0], pentagon[0], pentagon[4]}));
  star_tips[1] =
    Polygon(axom::Array<Point3D>({pentagram[1], pentagon[2], pentagon[1]}));
  star_tips[2] =
    Polygon(axom::Array<Point3D>({pentagram[2], pentagon[4], pentagon[3]}));
  star_tips[3] =
    Polygon(axom::Array<Point3D>({pentagram[3], pentagon[1], pentagon[0]}));
  star_tips[4] =
    Polygon(axom::Array<Point3D>({pentagram[4], pentagon[3], pentagon[2]}));

  // Get a non-convex, but non-self intersecting version of the star
  // clang-format off
  axom::Array<Point3D> full_star_vertices({pentagram[0], pentagon[0],
                                           pentagram[3], pentagon[1],
                                           pentagram[1], pentagon[2],
                                           pentagram[4], pentagon[3],
                                           pentagram[2], pentagon[4],
  });
  // clang-format on
  Polygon full_star(full_star_vertices);

  // Add up components of the star
  double star_components = 0;
  star_components += winding_number(query, pentagon);
  for(int i = 0; i < 5; ++i)
    star_components += winding_number(query, star_tips[i]);

  // Triangulate the pentagram directly
  Polygon pentagram_tris[5 - 2];
  for(int i = 0; i < 5 - 2; ++i)
    pentagram_tris[i] = Polygon(
      axom::Array<Point3D>({pentagram[0], pentagram[i + 1], pentagram[i + 2]}));

  // Add up components of the pentagram triangulation
  double pentagram_triangulation = 0;
  for(int i = 0; i < 5 - 2; ++i)
    pentagram_triangulation += winding_number(query, pentagram_tris[i]);

  // Test equality of various iterations of star and pentagram
  Point3D queries[5] = {
    Point3D {0.0, 4.0, 1.0},
    Point3D {-1.0, 2.0, 2.0},
    Point3D {0.0, -5.0, 3.0},
    Point3D {0.0, 0.0, 4.0},
    Point3D {3.0, 2.0, 5.0},
  };

  for(int n = 0; n < 5; ++n)
  {
    // Nonoverlapping star should equal sum of tips + pentagon
    EXPECT_NEAR(star_components, winding_number(query, full_star), 1e-10);

    // Pentagram should be full star + extra pentagon
    EXPECT_NEAR(
      winding_number(query, full_star) + winding_number(query, pentagon),
      pentagram_triangulation,
      1e-10);

    // Pentagram should equal it's triangulation
    EXPECT_NEAR(winding_number(query, pentagram), pentagram_triangulation, 1e-10);
  }
}

//------------------------------------------------------------------------------
TEST(primal_solid_angle, selfintersecting_quadrilateral)
{
  using Point3D = primal::Point<double, 3>;
  using Vector3D = primal::Vector<double, 3>;
  using Triangle = primal::Triangle<double, 3>;
  using Polygon = primal::Polygon<double, 3>;

  Point3D queries[5] = {
    Point3D {0.5, 0.5, 1.0},
    Point3D {-1.0, 2.0, 2.0},
    Point3D {0.0, -5.0, 3.0},
    Point3D {0.0, 0.0, 4.0},
    Point3D {3.0, 2.0, 5.0},
  };

  axom::Array<Point3D> square_vertices({
    Point3D {1.0, 0.0, 0.0},
    Point3D {1.0, 1.0, 0.0},
    Point3D {0.0, 1.0, 0.0},
    Point3D {0.0, 0.0, 0.0},
  });
  Polygon square(square_vertices);

  // Test on a self-intersecting "square" shape
  axom::Array<Point3D> squareish_vertices({
    Point3D {1.0, 0.0, 0.0},
    Point3D {1.0, 1.0, 0.0},
    Point3D {0.0, 1.0, 0.0},
    Point3D {0.0, 0.0, 0.0},
    Point3D {1.0, 0.0, 0.0},
    Point3D {1.0, 1.0, 0.0},
    Point3D {0.0, 1.0, 0.0},
    Point3D {0.0, 0.0, 0.0},
  });
  Polygon squareish(squareish_vertices);

  // Area of overlapping square should be twice that of nonoverlapping
  for(int n = 0; n < 5; ++n)
    EXPECT_NEAR(2 * winding_number(queries[n], square),
                winding_number(queries[n], squareish),
                1e-10);

  // Test on non-uniformly oriented hourglass shape
  axom::Array<Point3D> hourglass_vertices({
    Point3D {1.0, 0.0, 0.0},
    Point3D {0.0, 1.0, 0.0},
    Point3D {1.0, 1.0, 0.0},
    Point3D {0.0, 0.0, 0.0},
  });
  Polygon hourglass(hourglass_vertices);

  Triangle top_bulb(Point3D {1.0, 0.0, 0.0},
                    Point3D {0.5, 0.5, 0},
                    Point3D {0.0, 0.0, 0.0});
  Triangle bot_bulb(Point3D {0.0, 1.0, 0.0},
                    Point3D {1.0, 1.0, 0},
                    Point3D {0.5, 0.5, 0.0});

  // Area of overlapping square should be twice that of nonoverlapping
  for(int n = 0; n < 5; ++n)
    EXPECT_NEAR(winding_number(queries[n], top_bulb) +
                  winding_number(queries[n], bot_bulb),
                winding_number(queries[n], hourglass),
                1e-10);
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  int result = RUN_ALL_TESTS();

  return result;
}

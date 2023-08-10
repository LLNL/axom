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
    {
      // Means query point is interior
      EXPECT_NEAR(0.25 * M_1_PI * solid_angle,
                  winding_number(queries[n], tri),
                  1e-10);
    }
    else
    {
      // Means query point is exterior
      EXPECT_NEAR(-0.25 * M_1_PI * solid_angle,
                  winding_number(queries[n], tri),
                  1e-10);
    }
  }

  // Test with simple degenerate triangle
  Triangle deg(Point3D {1.0, 2.0, 3.0},
               Point3D {2.0, 3.0, 4.0},
               Point3D {1.0, 2.0, 3.0});
  for(int n = 0; n < 5; ++n)
  {
    EXPECT_DOUBLE_EQ(0.0, winding_number(queries[n], deg));
  }

  // Test with nondegenerate triangle, zero winding number
  Point3D coplanar {-1.0, 1.0, 1.0};
  EXPECT_NEAR(0.0, winding_number(coplanar, octant), 1e-10);

  // Test with coincident point
  Point3D coincident {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
  EXPECT_NEAR(0.0, winding_number(coincident, octant), 1e-10);

  // Test with point at vertex
  Point3D vertex {0.0, 1.0, 0.0};
  EXPECT_NEAR(0.0, winding_number(vertex, octant), 1e-10);
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
  {
    EXPECT_NEAR(winding_number(queries[n], square),
                winding_number(queries[n], square_tri1) +
                  winding_number(queries[n], square_tri2),
                1e-10);
  }

  // Test on coplanar points
  Point3D troubled_queries[5] = {
    Point3D {0.5, 0.5, 0.0},
    Point3D {1.0, 0.5, 0.0},
    Point3D {0.5, 1.0, 0.0},
    Point3D {1.5, 1.5, 0.0},
    Point3D {0.0, 0.0, 0.0},
  };

  for(int n = 0; n < 5; ++n)
  {
    EXPECT_NEAR(winding_number(troubled_queries[n], square), 0, 1e-10);
  }
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
  {
    EXPECT_NEAR(winding_number(queries[n], pentagon),
                winding_number(queries[n], pentagon_tri1) +
                  winding_number(queries[n], pentagon_tri2),
                1e-10);
  }

  // Test on coplanar points
  Point3D troubled_queries[5] = {
    Point3D {0.5, 0.5, 0.0},
    Point3D {1.0, 0.5, 0.0},
    Point3D {0.5, 1.0, 0.0},
    Point3D {1.5, 1.5, 0.0},
    Point3D {0.0, 0.0, 0.0},
  };

  for(int n = 0; n < 5; ++n)
  {
    EXPECT_NEAR(winding_number(troubled_queries[n], pentagon), 0, 1e-10);
  }
}

//------------------------------------------------------------------------------
TEST(primal_solid_angle, degenerate_polygon)
{
  using Point3D = primal::Point<double, 3>;
  using Vector3D = primal::Vector<double, 3>;
  using Polygon = primal::Polygon<double, 3>;

  Vector3D v1 = Vector3D({0.0, 1.0, 2.0}).unitVector();
  Vector3D v2 = Vector3D({2.0, -1.0, 0.5}).unitVector();

  Polygon good_pentagon(5), bad_pentagon(9);
  double good_angles[5] = {0.0, 0.5, 1.2, 3.0, 5.0};

  // Test with final vertex coincident with initial vertex
  double bad_angles[9] = {0.0, 0.5, 0.5, 1.2, 3.0, 3.0, 3.0, 5.0, 0.0};

  // Add vertices to good polygon
  for(int i = 0; i < 5; ++i)
  {
    good_pentagon.addVertex(
      Point3D {cos(good_angles[i]) * v1[0] + sin(good_angles[i]) * v2[0],
               cos(good_angles[i]) * v1[1] + sin(good_angles[i]) * v2[1],
               cos(good_angles[i]) * v1[2] + sin(good_angles[i]) * v2[2]});
  }

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
  {
    bad_pentagon.addVertex(
      Point3D {cos(bad_angles[i]) * v1[0] + sin(bad_angles[i]) * v2[0],
               cos(bad_angles[i]) * v1[1] + sin(bad_angles[i]) * v2[1],
               cos(bad_angles[i]) * v1[2] + sin(bad_angles[i]) * v2[2]});
  }

  Point3D queries[5] = {
    Point3D {0.0, 4.0, 1.0},
    Point3D {-1.0, 2.0, 2.0},
    Point3D {0.0, -5.0, 3.0},
    Point3D {0.0, 0.0, 4.0},
    Point3D {3.0, 2.0, 5.0},
  };

  EXPECT_NEAR(good_pentagon.area(), bad_pentagon.area(), 1e-10);
  for(int n = 0; n < 5; ++n)
  {
    EXPECT_NEAR(winding_number(queries[n], good_pentagon),
                winding_number(queries[n], bad_pentagon),
                1e-10);
  }

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
  {
    EXPECT_NEAR(winding_number(queries[n], Polygon(good_square_vertices)),
                winding_number(queries[n], Polygon(bad_square_vertices)),
                1e-10);
  }
}

//------------------------------------------------------------------------------
TEST(primal_solid_angle, selfintersecting_star)
{
  using Point3D = primal::Point<double, 3>;
  using Vector3D = primal::Vector<double, 3>;
  using Polygon = primal::Polygon<double, 3>;

  Point3D single_query {1, 2, 3};

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
  {
    pentagram.addVertex(
      Point3D {cos(outer_angles[i]) * v1[0] + sin(outer_angles[i]) * v2[0],
               cos(outer_angles[i]) * v1[1] + sin(outer_angles[i]) * v2[1],
               cos(outer_angles[i]) * v1[2] + sin(outer_angles[i]) * v2[2]});
  }

  // Construct the inner pentagon
  for(int i = 0; i < 5; ++i)
  {
    pentagon.addVertex(Point3D {
      r0 * cos(inner_angles[i]) * v1[0] + r0 * sin(inner_angles[i]) * v2[0],
      r0 * cos(inner_angles[i]) * v1[1] + r0 * sin(inner_angles[i]) * v2[1],
      r0 * cos(inner_angles[i]) * v1[2] + r0 * sin(inner_angles[i]) * v2[2]});
  }

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
  star_components += winding_number(single_query, pentagon);
  for(int i = 0; i < 5; ++i)
  {
    star_components += winding_number(single_query, star_tips[i]);
  }

  // Triangulate the pentagram directly
  Polygon pentagram_tris[5 - 2];
  for(int i = 0; i < 5 - 2; ++i)
  {
    pentagram_tris[i] = Polygon(
      axom::Array<Point3D>({pentagram[0], pentagram[i + 1], pentagram[i + 2]}));
  }

  // Add up components of the pentagram triangulation
  double pentagram_triangulation = 0;
  for(int i = 0; i < 5 - 2; ++i)
  {
    pentagram_triangulation += winding_number(single_query, pentagram_tris[i]);
  }

  for(int n = 0; n < 5; ++n)
  {
    // Nonoverlapping star should equal sum of tips + pentagon
    EXPECT_NEAR(star_components, winding_number(single_query, full_star), 1e-10);

    // Pentagram should be full star + extra pentagon
    EXPECT_NEAR(winding_number(single_query, full_star) +
                  winding_number(single_query, pentagon),
                pentagram_triangulation,
                1e-10);

    // Pentagram should equal it's triangulation
    EXPECT_NEAR(winding_number(single_query, pentagram),
                pentagram_triangulation,
                1e-10);
  }
}

//------------------------------------------------------------------------------
TEST(primal_solid_angle, selfintersecting_quadrilateral)
{
  using Point3D = primal::Point<double, 3>;
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
  {
    EXPECT_NEAR(2 * winding_number(queries[n], square),
                winding_number(queries[n], squareish),
                1e-10);
  }

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
  {
    EXPECT_NEAR(winding_number(queries[n], top_bulb) +
                  winding_number(queries[n], bot_bulb),
                winding_number(queries[n], hourglass),
                1e-10);
  }
}

//------------------------------------------------------------------------------
TEST(primal_solid_angle, planar_bezierpatch)
{
  return;

  using Point3D = primal::Point<double, 3>;
  using Vector3D = primal::Vector<double, 3>;
  using Polygon = primal::Polygon<double, 3>;
  using BezierPatch = primal::BezierPatch<double, 3>;

  // Define normal vector for the quadrilateral
  Vector3D v1 = Vector3D({0.0, 1.0, 2.0}).unitVector();
  Vector3D v2 = Vector3D({2.0, -1.0, 0.5}).unitVector();
  Vector3D v3 = Vector3D::cross_product(v1, v2);

  Polygon quad(4);

  double angles[5] = {1.0, 1.5, 2.5, 3.0};
  // Add vertices to quadrilateral
  for(int i = 0; i < 4; ++i)
  {
    quad.addVertex(Point3D {cos(angles[i]) * v1[0] + sin(angles[i]) * v2[0],
                            cos(angles[i]) * v1[1] + sin(angles[i]) * v2[1],
                            cos(angles[i]) * v1[2] + sin(angles[i]) * v2[2]});
  }

  // Construct a first order Bezier patch out of the same vertices
  Point3D controlPoints[4] = {quad[1], quad[0], quad[2], quad[3]};
  BezierPatch quad_patch(controlPoints, 1, 1);
  const bool isValidPatch = false;

  Point3D queries[8] = {Point3D {0.0, 4.0, 1.0},
                        Point3D {-1.0, 2.0, 2.0},
                        Point3D {0.0, -5.0, 3.0},
                        Point3D {0.0, 0.0, 4.0},
                        Point3D {3.0, 2.0, 5.0},
                        Point3D((v3 + 0.1 * v1 + 0.1 * v2).array()),
                        Point3D((-v3 + 0.1 * v1 + 0.1 * v2).array()),
                        Point3D((10.0 * v1 + 10.0 * v2).array())};

  // Should be equal with both kinds of primitive
  for(int n = 0; n < 8; ++n)
  {
    EXPECT_NEAR(winding_number(queries[n], quad),
                winding_number_recursive(queries[n], quad_patch, isValidPatch),
                1e-10);
  }

  // The winding numbers computed should be the same even if
  //  we use quadrature instead of resorting to the direct formula.
  // Ensure this by shrinking our tolerances.
  const double quad_tol = 1e-10;
  const double edge_tol = 1e-10;
  const double EPS = 0;
  for(int n = 0; n < 8; ++n)
  {
    EXPECT_NEAR(winding_number_recursive(queries[n], quad_patch, isValidPatch),
                winding_number_recursive(queries[n],
                                         quad_patch,
                                         isValidPatch,
                                         quad_tol,
                                         edge_tol,
                                         EPS),
                1e-10);
  }
}

//------------------------------------------------------------------------------
TEST(primal_solid_angle, bezierpatch_sphere)
{
  return;

  using Point3D = primal::Point<double, 3>;
  using Vector3D = primal::Vector<double, 3>;
  using BPatch = primal::BezierPatch<double, 3>;

  double rt2 = sqrt(2), rt3 = sqrt(3), rt6 = sqrt(6);

  // Define the nodes and weights for one of six rational, biquartic Bezier patches
  //  that compose the unit sphere. These will be rotated to form the other 5.
  // Nodes and weights obtained from the technical report
  // "Tiling the Sphere with Rational Bezier Patches",
  //  James E. Cobb, University of Utah, 1988

  // clang-format off
  axom::Array<Point3D> node_data = {
    Point3D {4*(1-rt3),     4*(1-rt3),     4*(1-rt3)}, Point3D {rt2*(rt3-4),            -rt2, rt2*(rt3-4)}, Point3D {4*(1-2*rt3)/3,   0, 4*(1-2*rt3)/3}, Point3D {rt2*(rt3-4),           rt2,   rt2*(rt3-4)}, Point3D {4*(1-rt3),     4*(rt3-1),     4*(1-rt3)},
    Point3D {     -rt2, rt2*(rt3 - 4), rt2*(rt3 - 4)}, Point3D {(2-3*rt3)/2,     (2-3*rt3)/2,  -(rt3+6)/2}, Point3D {rt2*(2*rt3-7)/3, 0,      -5*rt6/3}, Point3D {(2-3*rt3)/2,   (3*rt3-2)/2,    -(rt3+6)/2}, Point3D {     -rt2,   rt2*(4-rt3),   rt2*(rt3-4)},
    Point3D {        0, 4*(1-2*rt3)/3, 4*(1-2*rt3)/3}, Point3D {          0, rt2*(2*rt3-7)/3,    -5*rt6/3}, Point3D {0,               0,   4*(rt3-5)/3}, Point3D {          0, rt2*(7-2*rt3)/3,    -5*rt6/3}, Point3D {        0, 4*(2*rt3-1)/3, 4*(1-2*rt3)/3},
    Point3D {      rt2, rt2*(rt3 - 4), rt2*(rt3 - 4)}, Point3D {(3*rt3-2)/2,     (2-3*rt3)/2,  -(rt3+6)/2}, Point3D {rt2*(7-2*rt3)/3, 0,      -5*rt6/3}, Point3D {(3*rt3-2)/2,   (3*rt3-2)/2,    -(rt3+6)/2}, Point3D {      rt2,   rt2*(4-rt3),   rt2*(rt3-4)},
    Point3D {4*(rt3-1),     4*(1-rt3),     4*(1-rt3)}, Point3D {rt2*(4-rt3),            -rt2, rt2*(rt3-4)}, Point3D {4*(2*rt3-1)/3,   0, 4*(1-2*rt3)/3}, Point3D {rt2*(4-rt3),           rt2,   rt2*(rt3-4)}, Point3D {4*(rt3-1),     4*(rt3-1),     4*(1-rt3)}};

  axom::Array<double> weight_data = {
         4*(3-rt3), rt2*(3*rt3-2),   4*(5-rt3)/3, rt2*(3*rt3-2),     4*(3-rt3),
     rt2*(3*rt3-2),     (rt3+6)/2, rt2*(rt3+6)/3,     (rt3+6)/2, rt2*(3*rt3-2),
       4*(5-rt3)/3, rt2*(rt3+6)/3, 4*(5*rt3-1)/9, rt2*(rt3+6)/3,   4*(5-rt3)/3,
     rt2*(3*rt3-2),     (rt3+6)/2, rt2*(rt3+6)/3,     (rt3+6)/2, rt2*(3*rt3-2),
         4*(3-rt3), rt2*(3*rt3-2),   4*(5-rt3)/3, rt2*(3*rt3-2),     4*(3-rt3)};
  // clang-format on

  BPatch sphere_faces[6];
  for(int n = 0; n < 6; ++n)
  {
    sphere_faces[n].setOrder(4, 4);
    sphere_faces[n].makeRational();
  }

  sphere_faces[0].setOrder(4, 4);
  for(int i = 0; i < 5; ++i)
  {
    for(int j = 0; j < 5; ++j)
    {
      const int idx = 5 * i + j;
      for(int n = 0; n < 6; ++n)
        sphere_faces[n].setWeight(i, j, weight_data[idx]);

      // Set up each face by rotating one of the patch faces
      sphere_faces[0](i, j)[0] = node_data[idx][1];
      sphere_faces[0](i, j)[1] = node_data[idx][0];
      sphere_faces[0](i, j)[2] = node_data[idx][2];
      sphere_faces[0](i, j).array() /= weight_data[idx];

      sphere_faces[1](i, j)[0] = -node_data[idx][0];
      sphere_faces[1](i, j)[1] = -node_data[idx][1];
      sphere_faces[1](i, j)[2] = -node_data[idx][2];
      sphere_faces[1](i, j).array() /= weight_data[idx];

      sphere_faces[2](i, j)[0] = node_data[idx][2];
      sphere_faces[2](i, j)[1] = node_data[idx][1];
      sphere_faces[2](i, j)[2] = node_data[idx][0];
      sphere_faces[2](i, j).array() /= weight_data[idx];

      sphere_faces[3](i, j)[0] = -node_data[idx][1];
      sphere_faces[3](i, j)[1] = -node_data[idx][2];
      sphere_faces[3](i, j)[2] = -node_data[idx][0];
      sphere_faces[3](i, j).array() /= weight_data[idx];

      sphere_faces[4](i, j)[0] = node_data[idx][0];
      sphere_faces[4](i, j)[1] = node_data[idx][2];
      sphere_faces[4](i, j)[2] = node_data[idx][1];
      sphere_faces[4](i, j).array() /= weight_data[idx];

      sphere_faces[5](i, j)[0] = -node_data[idx][2];
      sphere_faces[5](i, j)[1] = -node_data[idx][0];
      sphere_faces[5](i, j)[2] = -node_data[idx][1];
      sphere_faces[5](i, j).array() /= weight_data[idx];
    }
  }

  // Split each surface into valid subsurfaces
  axom::Array<BPatch> valid_subsurfaces;
  for(int n = 0; n < 6; ++n)
  {
    split_to_convex_shallow(sphere_faces[n], valid_subsurfaces);
  }
  const int num_subsurfaces = valid_subsurfaces.size();
  const bool isValidSubsurface = true;

  // Iterate over points of interest, i.e. close to a boundary.
  //  Specifically need to exercise the ray casting steps
  const int idx = 3;  // Select a difficult face
  Vector3D query_directions[11] = {
    Vector3D(valid_subsurfaces[idx].evaluate(0.50, 0.50)),        // 0
    Vector3D(valid_subsurfaces[idx].evaluate(0.99, 0.95)),        // 1
    Vector3D(valid_subsurfaces[idx].evaluate(0.999, 0.995)),      // 2
    Vector3D(valid_subsurfaces[idx].evaluate(0.9999, 0.9995)),    // 3
    Vector3D(valid_subsurfaces[idx].evaluate(0.99999, 0.99995)),  // 4
    Vector3D(valid_subsurfaces[idx].evaluate(1.0, 0.9999)),       // 5
    Vector3D(valid_subsurfaces[idx].evaluate(0.9999, 1.0)),       // 6
    Vector3D(valid_subsurfaces[idx].evaluate(1.0, 1.0)),          // 7
    Vector3D(valid_subsurfaces[idx].evaluate(0.01, 0.95)),        // 8
    Vector3D(valid_subsurfaces[idx].evaluate(0.99, 0.05)),        // 9
    Vector3D(valid_subsurfaces[idx].evaluate(0.01, 0.05))};       // 10

  const double edge_tol = 1e-6;
  const double edge_offset = 1e-5;

  const double quad_tol = 1e-5;
  const double EPS = 1e-8;

  // Test some easy cases
  auto origin = Point3D({0.0, 0.0, 0.0});
  auto near_origin = Point3D({0.1, -0.2, 0.15});

  double origin_wn = 0.0, near_origin_wn = 0.0;
  for(int k = 0; k < num_subsurfaces; ++k)
  {
    origin_wn += winding_number_recursive(origin,
                                          valid_subsurfaces[k],
                                          isValidSubsurface,
                                          edge_tol,
                                          quad_tol,
                                          EPS);
    near_origin_wn += winding_number_recursive(near_origin,
                                               valid_subsurfaces[k],
                                               isValidSubsurface,
                                               edge_tol,
                                               quad_tol,
                                               EPS);
  }
  EXPECT_NEAR(origin_wn, 1.0, num_subsurfaces * quad_tol);
  EXPECT_NEAR(near_origin_wn, 1.0, num_subsurfaces * quad_tol);

  for(int i = 0; i < 11; ++i)
  {
    // Pick point close to the surface
    auto far_query = Point3D(10 * query_directions[i].array());

    double far_wn = 0.0;
    for(int k = 0; k < num_subsurfaces; ++k)
    {
      far_wn += winding_number_recursive(far_query,
                                         valid_subsurfaces[k],
                                         isValidSubsurface,
                                         edge_tol,
                                         quad_tol,
                                         EPS);
    }
    EXPECT_NEAR(far_wn, 0.0, num_subsurfaces * quad_tol);
  }

  // Iterate over difficult query directions for very close points
  for(int i = 0; i < 11; ++i)
  {
    std::cout << std::endl << i << std::endl;
    // Pick point close to the surface
    auto inner_query = Point3D((1.0 - edge_offset) * query_directions[i].array());
    auto outer_query = Point3D((1.0 + edge_offset) * query_directions[i].array());

    // Iterate over the patches that compose the sphere
    double inner_wn = 0;
    double inner_wn_old = 0.0;
    for(int k = 0; k < num_subsurfaces; ++k)
    {
      inner_wn += winding_number_recursive(inner_query,
                                           valid_subsurfaces[k],
                                           isValidSubsurface,
                                           edge_tol,
                                           quad_tol,
                                           EPS);

      //inner_wn_old += winding_number_old(inner_query,
      //                                   valid_subsurfaces[k],
      //                                   edge_tol,
      //                                   quad_tol,
      //                                   EPS);

      //if(!axom::utilities::isNearlyEqual(inner_wn, inner_wn_old, 0.1))
      //{
      //  printf("BPOAOFIJAEOFHFAAOEOF\n");
      //  //printf("%f\n", inner_wn - inner_wn_old);
      //}
    }
    EXPECT_NEAR(inner_wn, 1.0, num_subsurfaces * quad_tol);

    // Iterate over the patches that compose the sphere
    double outer_wn = 0;
    for(int k = 0; k < num_subsurfaces; ++k)
    {
      outer_wn += winding_number_recursive(outer_query,
                                           valid_subsurfaces[k],
                                           isValidSubsurface,
                                           edge_tol,
                                           quad_tol,
                                           EPS);
    }
    EXPECT_NEAR(outer_wn, 0.0, num_subsurfaces * quad_tol);

    // Pick a point on the surface too.
    //  Regardless of what tolerances are picked, the winding number
    //  should lie between the values on either side when rounded
    auto coincident_query = Point3D(query_directions[i].array());
    double coincident_wn = 0.0;
    for(int k = 0; k < valid_subsurfaces.size(); ++k)
    {
      coincident_wn += winding_number_recursive(coincident_query,
                                                valid_subsurfaces[k],
                                                isValidSubsurface,
                                                edge_tol,
                                                quad_tol,
                                                EPS);
    }
    EXPECT_NEAR(coincident_wn, 0.5, num_subsurfaces * quad_tol);
  }
}

TEST(primal_solid_angle, bezierpatch_degenerate_sphere)
{
  return;

  using Point2D = primal::Point<double, 2>;
  using Point3D = primal::Point<double, 3>;
  using Vector3D = primal::Vector<double, 3>;
  using BCurve = primal::BezierCurve<double, 2>;
  using BPatch = primal::BezierPatch<double, 3>;

  auto rotate_curve_to_patches = [](BCurve curve, BPatch patches[]) -> void {
    const int ord = curve.getOrder();
    curve.makeRational();

    for(int n = 0; n < 4; ++n)
    {
      patches[n].setOrder(ord, 2);
      patches[n].makeRational();
    }

    for(int p = 0; p <= ord; ++p)
    {
      auto node = curve[p];
      patches[0](p, 0) = Point3D {node[0], 0.0, node[1]};
      patches[0](p, 1) = Point3D {node[0], node[0], node[1]};
      patches[0](p, 2) = Point3D {0.0, node[0], node[1]};

      patches[1](p, 0) = Point3D {0.0, node[0], node[1]};
      patches[1](p, 1) = Point3D {-node[0], node[0], node[1]};
      patches[1](p, 2) = Point3D {-node[0], 0.0, node[1]};

      patches[2](p, 0) = Point3D {-node[0], 0.0, node[1]};
      patches[2](p, 1) = Point3D {-node[0], -node[0], node[1]};
      patches[2](p, 2) = Point3D {0.0, -node[0], node[1]};

      patches[3](p, 0) = Point3D {0.0, -node[0], node[1]};
      patches[3](p, 1) = Point3D {node[0], -node[0], node[1]};
      patches[3](p, 2) = Point3D {node[0], 0.0, node[1]};

      for(int n = 0; n < 4; ++n)
      {
        patches[n].setWeight(p, 0, curve.getWeight(p));
        patches[n].setWeight(p, 1, curve.getWeight(p) / std::sqrt(2));
        patches[n].setWeight(p, 2, curve.getWeight(p));
      }
    }
  };

  // Make a sphere out of 8 degenerate Bezier Patches

  // Define rational Beziers curve along the x-z axis and rotate them
  Point2D curve1_nodes[] = {Point2D {0, 1}, Point2D {1, 1}, Point2D {1, 0}};
  double curve1_weights[] = {1, 1 / std::sqrt(2), 1};
  BCurve curve1(curve1_nodes, curve1_weights, 2);
  BPatch patches1[4];
  rotate_curve_to_patches(curve1, patches1);

  Point2D curve2_nodes[] = {Point2D {1, 0}, Point2D {1, -1}, Point2D {0, -1}};
  double curve2_weights[] = {1, 1 / std::sqrt(2), 1};
  BCurve curve2(curve2_nodes, curve2_weights, 2);
  BPatch patches2[4];
  rotate_curve_to_patches(curve2, patches2);

  // Split each surface into valid subsurfaces
  axom::Array<BPatch> valid_subsurfaces;
  for(int n = 0; n < 4; ++n)
  {
    split_to_convex_shallow(patches1[n], valid_subsurfaces);
    split_to_convex_shallow(patches2[n], valid_subsurfaces);
  }

  const int num_subsurfaces = valid_subsurfaces.size();
  const bool isValidSubsurface = true;

  // Iterate over points of interest, i.e. close to a boundary.
  //  Specifically need to exercise the ray casting steps
  int idx = 6;  // Select a difficult face
  Vector3D query_directions[11] = {
    Vector3D(valid_subsurfaces[idx].evaluate(0.50, 0.50)),        // 0
    Vector3D(valid_subsurfaces[idx].evaluate(0.99, 0.95)),        // 1
    Vector3D(valid_subsurfaces[idx].evaluate(0.999, 0.995)),      // 2
    Vector3D(valid_subsurfaces[idx].evaluate(0.9999, 0.9995)),    // 3
    Vector3D(valid_subsurfaces[idx].evaluate(0.99999, 0.99995)),  // 4
    Vector3D(valid_subsurfaces[idx].evaluate(1.0, 0.9999)),       // 5
    Vector3D(valid_subsurfaces[idx].evaluate(0.9999, 1.0)),       // 6
    Vector3D(valid_subsurfaces[idx].evaluate(1.0, 1.0)),          // 7
    Vector3D(valid_subsurfaces[idx].evaluate(0.01, 0.95)),        // 8
    Vector3D(valid_subsurfaces[idx].evaluate(0.99, 0.05)),        // 9
    Vector3D(valid_subsurfaces[idx].evaluate(0.01, 0.05))};       // 10

  const double edge_tol = 1e-6;
  const double edge_offset = 1e-5;

  const double quad_tol = 1e-5;
  const double EPS = 1e-8;

  // Test some easy cases
  auto origin = Point3D({0.0, 0.0, 0.0});
  auto near_origin = Point3D({0.1, -0.2, 0.15});

  double origin_wn = 0.0, near_origin_wn = 0.0;
  for(int k = 0; k < num_subsurfaces; ++k)
  {
    origin_wn += winding_number_recursive(origin,
                                          valid_subsurfaces[k],
                                          isValidSubsurface,
                                          edge_tol,
                                          quad_tol,
                                          EPS);
    near_origin_wn += winding_number_recursive(near_origin,
                                               valid_subsurfaces[k],
                                               isValidSubsurface,
                                               edge_tol,
                                               quad_tol,
                                               EPS);
  }
  EXPECT_NEAR(origin_wn, 1.0, num_subsurfaces * quad_tol);
  EXPECT_NEAR(near_origin_wn, 1.0, num_subsurfaces * quad_tol);

  for(int i = 0; i < 11; ++i)
  {
    auto far_query = Point3D(10 * query_directions[i].array());

    double far_wn = 0.0;
    for(int k = 0; k < num_subsurfaces; ++k)
    {
      far_wn += winding_number_recursive(far_query,
                                         valid_subsurfaces[k],
                                         isValidSubsurface,
                                         edge_tol,
                                         quad_tol,
                                         EPS);
    }
    EXPECT_NEAR(far_wn, 0.0, num_subsurfaces * quad_tol);
  }

  // Iterate over difficult query directions for very close points
  for(int i = 0; i < 11; ++i)
  {
    std::cout << i << std::endl;

    // Pick point close to the surface
    auto inner_query = Point3D((1.0 - edge_offset) * query_directions[i].array());
    auto outer_query = Point3D((1.0 + edge_offset) * query_directions[i].array());

    // Iterate over the patches that compose the sphere
    double inner_wn = 0;
    for(int k = 0; k < num_subsurfaces; ++k)
    {
      inner_wn += winding_number_recursive(inner_query,
                                           valid_subsurfaces[k],
                                           isValidSubsurface,
                                           edge_tol,
                                           quad_tol,
                                           EPS);
    }
    EXPECT_NEAR(inner_wn, 1.0, num_subsurfaces * quad_tol);

    // Iterate over the patches that compose the sphere
    double outer_wn = 0;
    for(int k = 0; k < num_subsurfaces; ++k)
    {
      outer_wn += winding_number_recursive(outer_query,
                                           valid_subsurfaces[k],
                                           isValidSubsurface,
                                           edge_tol,
                                           quad_tol,
                                           EPS);
    }
    EXPECT_NEAR(outer_wn, 0.0, num_subsurfaces * quad_tol);

    // Pick a point on the surface too.
    auto coincident_query = Point3D(query_directions[i].array());
    double coincident_wn = 0.0;
    for(int k = 0; k < valid_subsurfaces.size(); ++k)
    {
      coincident_wn += winding_number_recursive(coincident_query,
                                                valid_subsurfaces[k],
                                                isValidSubsurface,
                                                edge_tol,
                                                quad_tol,
                                                EPS);
    }
    EXPECT_NEAR(coincident_wn, 0.5, num_subsurfaces * quad_tol);
  }
}

TEST(primal_solid_angle, bezierpatch_empty_interior)
{
  return;

  using Point2D = primal::Point<double, 2>;
  using Point3D = primal::Point<double, 3>;
  using Vector3D = primal::Vector<double, 3>;
  using BPatch = primal::BezierPatch<double, 3>;

  // clang-format off
  axom::Array<Point3D> point_data = {
    Point3D {2.5,0.2500,0.000}, Point3D {2.5,0.2500,0.000}, Point3D {2.5,0.2500,0.000}, Point3D{2.5,0.2500,0.000},
    Point3D {2.7,1.1250,0.750}, Point3D {2.7,1.1250,0.750}, Point3D {2.7,1.1250,0.750}, Point3D{2.7,1.1250,0.750},
    Point3D {2.6,1.8125,1.125}, Point3D {2.6,1.8125,1.125}, Point3D {2.6,1.8125,1.125}, Point3D{2.6,1.8125,1.125},
    Point3D {2.5,2.5000,1.125}, Point3D {2.5,2.5000,1.125}, Point3D {2.5,2.5000,1.125}, Point3D{2.5,2.5000,1.125}
  };
  // clang-format on

  // Degnerate patch with zero surface area
  BPatch bad_patch(point_data, 3, 3);

  Point3D surface_point = bad_patch.evaluate(0.2, 0.8);
  Point3D other_point = Point3D {surface_point[0] + 0.01,
                                 surface_point[1] - 0.01,
                                 surface_point[2] + 0.01};

  EXPECT_NEAR(winding_number_recursive(surface_point, bad_patch), 0.0, 1e-5);
  EXPECT_NEAR(winding_number_recursive(other_point, bad_patch), 0.0, 1e-5);
}

TEST(primal_solid_angle, bezierpatch_nonsense_geometry)
{
  using Point2D = primal::Point<double, 2>;
  using Point3D = primal::Point<double, 3>;
  using Vector3D = primal::Vector<double, 3>;
  using BPatch = primal::BezierPatch<double, 3>;

  // Get a hand-coded control polygon
  BPatch test_patch(1, 5);
  test_patch(0, 0) = Point3D {-5.0, 3.0, 0.0};
  test_patch(0, 1) = Point3D {-4.0, 4.0, 0.0};
  test_patch(0, 2) = Point3D {-2.0, 5.0, 0.0};
  test_patch(0, 3) = Point3D {0.0, 5.5, 0.0};
  test_patch(0, 4) = Point3D {2.0, 5.5, 0.0};
  test_patch(0, 5) = Point3D {4.5, 4.5, 0.0};

  test_patch(1, 0) = Point3D {-7.0, 3.0, 2.0};
  test_patch(1, 1) = Point3D {-4.0 - 2.0, 4.0, 2.0};
  test_patch(1, 2) = Point3D {-2.0 - 2.0, 5.0, 2.0};
  test_patch(1, 3) = Point3D {0.0 - 2.0, 5.5, 2.0};
  test_patch(1, 4) = Point3D {2.0 - 2.0, 5.5, 2.0};
  test_patch(1, 5) = Point3D {4.5 - 2.0, 4.5, 2.0};

  // clang-format off
  // Create literally random control points
  axom::Array<Point3D> point_data = {
    Point3D {0.4, 0.1, 0.8},
    Point3D {0.5, 0. , 0.3},
    Point3D {0.8, 0.8, 0.7},
    Point3D {0.2, 0.7, 0.6},
    Point3D {0.1, 0. , 0.6},
    Point3D {0.6, 0.1, 0.9},
    Point3D {0.3, 0.3, 0. },
    Point3D {0.1, 0.5, 0.3},
    Point3D {0.2, 0.4, 0.6},
    Point3D {0.5, 0.8, 0.1},
    Point3D {0.1, 0.9, 0.5},
    Point3D {0.1, 0.6, 0.7},
    Point3D {0.4, 1. , 0.4},
    Point3D {0.6, 0.4, 0.7},
    Point3D {0.6, 0.8, 0.2},
    Point3D {0.8, 0.3, 0.3}
  };
  // clang-format on

  BPatch nonsense_patch(point_data, 3, 3);

  // Split each surface into valid subsurfaces
  axom::Array<BPatch> valid_subsurfaces;
  split_to_convex_shallow(nonsense_patch, valid_subsurfaces);

  //python_print(nonsense_patch);
  //for(int n = 0; n < valid_subsurfaces.size(); ++n)
  //  python_print(valid_subsurfaces[n]);

  FILE* cmap_file = fopen(
    "../../../../../axom_aux/cmap_data/nonsense_shape/"
    "single_surface.csv",
    "w");

  int npts_u = 150;
  int npts_v = 150;
  double umin = -0.1, umax = 0.1;
  double vmin = -0.1, vmax = 0.1;

  // Define a linear patch that represents this plane
  Vector3D e1 = Vector3D({1, 0, 0}).unitVector();
  Vector3D e2 = Vector3D({0, 1, 0}).unitVector();
  Vector3D center = Vector3D {0.3733333333333331, 0.4466666666666663, 0.5};

  Point3D controlPoints[4] = {Point3D((center + umin * e1 + vmin * e2).array()),
                              Point3D((center + umax * e1 + vmin * e2).array()),
                              Point3D((center + umin * e1 + vmax * e2).array()),
                              Point3D((center + umax * e1 + vmax * e2).array())};
  BPatch quad_patch(controlPoints, 1, 1);

  Point3D bad_query {0.37333333333333307, 0.44666666666666632, 0.50000000000000000};
  //winding_number(bad_query, valid_subsurfaces[13], true);

  for(double v = -0.10000000000000001; v <= vmax; v += (vmax - vmin) / npts_v)
  {
    printf("(u, v) = (u, %g)\n", v);
    for(double u = 0.034666666666666783; u <= umax; u += (umax - umin) / npts_u)
    {
      Point3D query((center + u * e1 + v * e2).array());

      double wn_casting = 0.0;
      double wn_recursive = 0.0;
      for(int k = 5; k < valid_subsurfaces.size(); ++k)
      {
        wn_casting += winding_number_casting(query, valid_subsurfaces[k], true);
        wn_recursive += winding_number_recursive(query, valid_subsurfaces[k]);

        if(!axom::utilities::isNearlyEqual(wn_casting, wn_recursive, 0.1))
        {
          printf("Casting:   %f\n", wn_casting);
          printf("Recursive: %f\n", wn_recursive);
          break;
        }
      }

      fprintf(cmap_file,
              "%.16f, %.16f, %.16f, %.16f, %.16f, %.16f, %.17f\n",
              u,
              v,
              query[0],
              query[1],
              query[2],
              wn_casting,
              wn_recursive);
    }
  }

  fclose(cmap_file);
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  int result = RUN_ALL_TESTS();

  return result;
}

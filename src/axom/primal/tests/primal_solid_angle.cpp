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
  Triangle tri(Point3D {1.0, 2.0, 3.0},
               Point3D {2.0, 3.0, 4.0},
               Point3D {-1.0, 2.0, 3.0});

  Point3D queries[5] = {
    Point3D {0.0, 0.0, 0.0},
    Point3D {1.0, 2.0, 0.0},
    Point3D {0.0, 5.0, 0.0},
    Point3D {0.0, 0.0, 5.0},
    Point3D {3.0, 2.0, 1.0},
  };

  for(int n = 0; n < 5; ++n)
  {
    double theta_a = Triangle(queries[n], tri[0], tri[1]).angle(0);
    double theta_b = Triangle(queries[n], tri[1], tri[2]).angle(0);
    double theta_c = Triangle(queries[n], tri[2], tri[0]).angle(0);

    double theta_s = 0.5 * (theta_a + theta_b + theta_c);

    // Solid angle from this formula is unsigned
    double solid_angle = M_1_PI *
      atan(sqrt(tan(0.5 * theta_s) * tan(0.5 * (theta_s - theta_a)) *
                tan(0.5 * (theta_s - theta_b)) * tan(0.5 * (theta_s - theta_c))));

    // Get sign from orientation.
    if(primal::orientation(queries[n], tri) == primal::ON_NEGATIVE_SIDE)
      // Means query point is interior
      EXPECT_NEAR(solid_angle, winding_number(queries[n], tri), 1e-10);
    else
      // Means query point is exterior
      EXPECT_NEAR(-solid_angle, winding_number(queries[n], tri), 1e-10);
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
TEST(primal_solid_angle, polygon)
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
TEST(primal_solid_angle, degenerate_polygon)
{
  using Point3D = primal::Point<double, 3>;
  using Vector3D = primal::Vector<double, 3>;
  using Triangle = primal::Triangle<double, 3>;
  using Polygon = primal::Polygon<double, 3>;

  Vector3D v1 = Vector3D({0.0, 1.0, 2.0}).unitVector();
  Vector3D v2 = Vector3D({2.0, -1.0, 0.5}).unitVector();

  Polygon good_poly(5), bad_poly(9);
  double good_angles[5] = {0.0, 0.5, 1.2, 3.0, 5.0};
  double bad_angles[9] = {0.0, 0.5, 0.5, 1.2, 3.0, 3.0, 3.0, 5.0, 0.0};

  // Add vertices to good polygon
  for(int i = 0; i < 5; ++i)
    good_poly.addVertex(
      Point3D {cos(good_angles[i]) * v1[0] + sin(good_angles[i]) * v2[0],
               cos(good_angles[i]) * v1[1] + sin(good_angles[i]) * v2[1],
               cos(good_angles[i]) * v1[2] + sin(good_angles[i]) * v2[2]});

  // Create bad polygon with degeneracies
  // Add point at angle 0
  bad_poly.addVertex(
    Point3D {cos(bad_angles[0]) * v1[0] + sin(bad_angles[0]) * v2[0],
             cos(bad_angles[0]) * v1[1] + sin(bad_angles[0]) * v2[1],
             cos(bad_angles[0]) * v1[2] + sin(bad_angles[0]) * v2[2]});
  // Add a midpoint between angles 0 and 1
  bad_poly.addVertex(Point3D::midpoint(
    bad_poly[0],
    Point3D {cos(bad_angles[1]) * v1[0] + sin(bad_angles[1]) * v2[0],
             cos(bad_angles[1]) * v1[1] + sin(bad_angles[1]) * v2[1],
             cos(bad_angles[1]) * v1[2] + sin(bad_angles[1]) * v2[2]}));
  // Add the rest of the vertices
  for(int i = 1; i < 9; ++i)
    bad_poly.addVertex(
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

  EXPECT_NEAR(good_poly.area(), bad_poly.area(), 1e-10);

  for(int n = 0; n < 5; ++n)
    EXPECT_NEAR(winding_number(queries[n], good_poly),
                winding_number(queries[n], bad_poly),
                1e-10);
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  int result = RUN_ALL_TESTS();

  return result;
}

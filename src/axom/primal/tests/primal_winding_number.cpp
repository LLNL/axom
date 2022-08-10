// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
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

TEST(primal_winding_number, containment_protocol)
{
  using Point2D = primal::Point<double, 2>;
  using Bezier = primal::BezierCurve<double, 2>;
  using CPolygon = primal::CurvedPolygon<double, 2>;

  // 8th order, closed curve with internal loop
  Point2D loop_nodes[] = {Point2D {0.0, 0.0},
                          Point2D {1.0, 0.0},
                          Point2D {1.0, 1.0},
                          Point2D {0.0, 1.0},
                          Point2D {0.0, 0.0},
                          Point2D {1.0, 0.0},
                          Point2D {1.0, 1.0},
                          Point2D {0.0, 1.0},
                          Point2D {0.0, 0.0}};
  Bezier loop_curve(loop_nodes, 8);
  CPolygon loop_poly;
  loop_poly.addEdge(loop_curve);

  // Inner loop is considered "interior" with nonzero protocol. Default behavior.
  bool nonzero = true;
  EXPECT_TRUE(in_curved_polygon(Point2D({0.5, 0.5}), loop_poly));
  EXPECT_TRUE(in_curved_polygon(Point2D({0.5, 0.5}), loop_poly, nonzero));

  // Inner loop is considered "exterior" with even/odd protocol
  nonzero = false;
  EXPECT_FALSE(in_curved_polygon(Point2D({0.5, 0.5}), loop_poly, nonzero));
}

TEST(primal_winding_number, simple_cases)
{
  using Point2D = primal::Point<double, 2>;
  using Bezier = primal::BezierCurve<double, 2>;
  using CPolygon = primal::CurvedPolygon<double, 2>;

  double abs_tol = 1e-8;
  double edge_tol = 1e-8;
  double EPS = 1e-50;

  // Simple cubic shape
  Point2D top_nodes[] = {Point2D {0.0, 0.0},
                         Point2D {0.0, 1.0},
                         Point2D {-1.0, 1.0},
                         Point2D {-1.0, 0.0}};
  Bezier top_curve(top_nodes, 3);

  Point2D bot_nodes[] = {Point2D {-1.0, 0.0},
                         Point2D {-1.0, -1.0},
                         Point2D {0.0, -1.0},
                         Point2D {0.0, 0.0}};
  Bezier bot_curve(bot_nodes, 3);
  Bezier simple_shape_edges[] = {top_curve, bot_curve};
  CPolygon simple_shape(simple_shape_edges, 2);

  // Check interior points
  for(int i = 1; i < 7; i++)
  {
    double offset = std::pow(10, -i);
    Point2D q_vert({-0.352, 0.72 - offset});
    Point2D q_horz({-0.352 - offset, 0.72});

    EXPECT_NEAR(winding_number(q_vert, simple_shape, edge_tol, EPS), 1.0, abs_tol);
    EXPECT_NEAR(winding_number(q_horz, simple_shape, edge_tol, EPS), 1.0, abs_tol);
  }

  // Check exterior points
  for(int i = 1; i < 7; i++)
  {
    double offset = std::pow(10, -i);
    Point2D q_vert({-0.352, 0.72 + offset});
    Point2D q_horz({-0.352 + offset, 0.72});

    EXPECT_NEAR(winding_number(q_vert, simple_shape, edge_tol, EPS), 0.0, abs_tol);
    EXPECT_NEAR(winding_number(q_horz, simple_shape, edge_tol, EPS), 0.0, abs_tol);
  }

  // Test that points on either side of cubic are offset by 1
  EXPECT_NEAR(
    winding_number(Point2D({-0.352, 0.72 - edge_tol * 2}), top_curve, edge_tol, EPS) -
      winding_number(Point2D({-0.352, 0.72 + edge_tol * 2}), top_curve, edge_tol, EPS),
    1,
    abs_tol);

  top_curve.reverseOrientation();
  EXPECT_NEAR(
    winding_number(Point2D({-0.352, 0.72 + edge_tol * 2}), top_curve, edge_tol, EPS) -
      winding_number(Point2D({-0.352, 0.72 - edge_tol * 2}), top_curve, edge_tol, EPS),
    1,
    abs_tol);
}

TEST(primal_winding_number, edge_cases)
{
  // Check the edge cases for the winding number closure formulation
  using Point2D = primal::Point<double, 2>;
  using Bezier = primal::BezierCurve<double, 2>;

  double abs_tol = 1e-10;
  double edge_tol = 1e-8;
  double EPS = 1e-50;

  // Line segment
  Point2D nodes1[] = {Point2D {0.0, 0.0}, Point2D {1.0, 1.0}};
  Bezier linear(nodes1, 1);

  // Cubic curve
  Point2D nodes2[] = {Point2D {0.0, 0.0},
                      Point2D {0.0, 1.0},
                      Point2D {-1.0, 1.0},
                      Point2D {-1.0, 0.0}};
  Bezier cubic(nodes2, 3);

  // Flipped Cubic curve
  Point2D fnodes2[] = {Point2D {0.0, 0.0},
                       Point2D {0.0, -1.0},
                       Point2D {-1.0, -1.0},
                       Point2D {-1.0, 0.0}};
  Bezier fcubic(fnodes2, 3);

  // Test on linear cases
  EXPECT_NEAR(winding_number(Point2D({-0.45, -0.45}), linear, edge_tol, EPS),
              0.0,
              abs_tol);
  EXPECT_NEAR(winding_number(Point2D({1.45, 1.45}), linear, edge_tol, EPS),
              0.0,
              abs_tol);

  // Test on cubic. If query is outside the range of the two endpoints, return 0
  EXPECT_NEAR(winding_number(Point2D({-2.5, 0.0}), cubic, edge_tol, EPS),
              0.0,
              abs_tol);
  EXPECT_NEAR(winding_number(Point2D({-2.5, 0.0}), fcubic, edge_tol, EPS),
              0.0,
              abs_tol);

  // If query point between the endpoints, return +0.5/-0.5
  EXPECT_NEAR(winding_number(Point2D({-0.5, 0.0}), cubic, edge_tol, EPS),
              0.5,
              abs_tol);
  EXPECT_NEAR(winding_number(Point2D({-0.5, 0.0}), fcubic, edge_tol, EPS),
              -0.5,
              abs_tol);

  cubic.reverseOrientation();
  fcubic.reverseOrientation();
  EXPECT_NEAR(winding_number(Point2D({-0.5, 0.0}), cubic, edge_tol, EPS),
              -0.5,
              abs_tol);
  EXPECT_NEAR(winding_number(Point2D({-0.5, 0.0}), fcubic, edge_tol, EPS),
              0.5,
              abs_tol);

  // Extra tests if initial and terminal tangent lines are colinear
  Point2D nodes3[] = {Point2D {0.1, 0.0},
                      Point2D {1.0, 0.0},
                      Point2D {0.0, 1.0},
                      Point2D {-1.0, 0.0},
                      Point2D {-0.1, 0.0}};
  Bezier quartic(nodes3, 4);

  // Tangent lines in opposite directions
  EXPECT_NEAR(winding_number(Point2D({0, 0}), quartic, edge_tol, EPS),
              0.5,
              abs_tol);

  // Flip the curve vertically
  quartic[2] = Point2D({0.0, -1.0});
  EXPECT_NEAR(winding_number(Point2D({0, 0}), quartic, edge_tol, EPS),
              -0.5,
              abs_tol);

  // Flip one of the tangent lines
  quartic[1] = Point2D({0.0, 0.0});
  EXPECT_NEAR(winding_number(Point2D({0, 0}), quartic, edge_tol, EPS),
              -0.5,
              abs_tol);

  // Flip vertically again
  quartic[2] = Point2D({0.0, 1.0});
  EXPECT_NEAR(winding_number(Point2D({0, 0}), quartic, edge_tol, EPS),
              0.5,
              abs_tol);
}

TEST(primal_winding_number, corner_cases)
{
  using Point2D = primal::Point<double, 2>;
  using Bezier = primal::BezierCurve<double, 2>;

  double abs_tol = 1e-4;
  double edge_tol = 1e-8;
  double EPS = 1e-50;

  // Line segment
  Point2D nodes1[] = {Point2D {0.0, 0.0}, Point2D {1.0, 1.0}};
  Bezier linear(nodes1, 1);

  // Cubic curve
  Point2D nodes2[] = {Point2D {0.0, 0.0},
                      Point2D {0.0, 1.0},
                      Point2D {-1.0, 1.0},
                      Point2D {-1.0, 0.0}};
  Bezier cubic(nodes2, 3);

  // At any point on a line, returns 0
  EXPECT_NEAR(  // Query on linear curve
    winding_number(Point2D({0.45, 0.45}), linear, edge_tol, EPS),
    0.0,
    abs_tol);
  EXPECT_NEAR(  // Query on endpoint of linear
    winding_number(Point2D({0.0, 0.0}), linear, edge_tol, EPS),
    0.0,
    abs_tol);
  EXPECT_NEAR(  // Query on endpoint of linear
    winding_number(Point2D({1.0, 1.0}), linear, edge_tol, EPS),
    0.0,
    abs_tol);

  EXPECT_NEAR(  // Query on initial endpoint of cubic
    winding_number(Point2D({-1.0, 0.0}), cubic, edge_tol, EPS),
    0.25,
    abs_tol);
  EXPECT_NEAR(  // Query on terminal endpoint of cubic
    winding_number(Point2D({-1.0, 0.0}), cubic, edge_tol, EPS),
    0.25,
    abs_tol);

  // The query is on the endpoint after one bisection
  EXPECT_NEAR(winding_number(Point2D({-0.5, 0.75}), cubic, edge_tol, EPS),
              0.312832962673,
              abs_tol);

  // The query is not on an endpoint after any number of bisections
  EXPECT_NEAR(winding_number(cubic.evaluate(0.4), cubic, edge_tol, EPS),
              0.310998027957,
              abs_tol);

  // Extra tests if initial and terminal tangent lines are colinear
  Point2D nodes3[] = {Point2D {0.1, 0.0},
                      Point2D {1.0, 0.0},
                      Point2D {0.0, 1.0},
                      Point2D {-1.0, 0.0},
                      Point2D {-0.1, 0.0}};
  Bezier quartic(nodes3, 4);

  // Bugs out because endpoint tangent lines point in opposite directions
  EXPECT_NEAR(winding_number(Point2D({0.1, 0}), quartic, edge_tol, EPS),
              0.5,
              abs_tol);
  EXPECT_NEAR(winding_number(Point2D({-0.1, 0}), quartic, edge_tol, EPS),
              0.5,
              abs_tol);

  // Flip the curve vertically
  quartic[2] = Point2D({0.0, -1.0});
  EXPECT_NEAR(winding_number(Point2D({0.1, 0}), quartic, edge_tol, EPS),
              -0.5,
              abs_tol);
  EXPECT_NEAR(winding_number(Point2D({-0.1, 0}), quartic, edge_tol, EPS),
              -0.5,
              abs_tol);

  // Flip one of the tangent lines
  quartic[1] = Point2D({0.0, 0.0});
  EXPECT_NEAR(winding_number(Point2D({0.1, 0}), quartic, edge_tol, EPS),
              0,
              abs_tol);
  EXPECT_NEAR(winding_number(Point2D({-0.1, 0}), quartic, edge_tol, EPS),
              -0.5,
              abs_tol);

  // Flip vertically again
  quartic[2] = Point2D({0.0, 1.0});
  EXPECT_NEAR(winding_number(Point2D({0.1, 0}), quartic, edge_tol, EPS),
              0,
              abs_tol);
  EXPECT_NEAR(winding_number(Point2D({-0.1, 0}), quartic, edge_tol, EPS),
              0.5,
              abs_tol);
}

TEST(primal_winding_number, self_intersecting_cases)
{
  using Point2D = primal::Point<double, 2>;
  using Bezier = primal::BezierCurve<double, 2>;

  double abs_tol = 1e-4;
  double edge_tol = 1e-8;
  double EPS = 1e-50;

  // Cubic curve with cusp
  Point2D nodes1[] = {Point2D {0.0, 0.0},
                      Point2D {1.0, 1.0},
                      Point2D {0.0, 1.0},
                      Point2D {1.0, 0.0}};
  Bezier cubic_cusp(nodes1, 3);

  EXPECT_NEAR(winding_number(Point2D({0.5, 0.5}), cubic_cusp, edge_tol, EPS),
              -0.75,
              abs_tol);

  // Bugs out because it has a duplicate endpoint after one split.
  EXPECT_NEAR(winding_number(Point2D({0.5, 0.75}), cubic_cusp, edge_tol, EPS),
              0.18716704181099889,
              abs_tol);
  EXPECT_NEAR(winding_number(Point2D({0.5, 1.5}), cubic_cusp, edge_tol, EPS),
              0.10241638235,
              abs_tol);

  // Cubic curve with loop
  Point2D nodes2[] = {Point2D {0.0, 0.0},
                      Point2D {2.0, 1.0},
                      Point2D {-1.0, 1.0},
                      Point2D {1.0, 0.0}};
  Bezier cubic_loop(nodes2, 3);
  EXPECT_NEAR(winding_number(Point2D({0.5, 0.3}), cubic_loop, edge_tol, EPS),
              0.327979130377,
              abs_tol);
  EXPECT_NEAR(winding_number(Point2D({0.5, 0.75}), cubic_loop, edge_tol, EPS),
              0.687167046798,
              abs_tol);

  // Closed cubic curve
  Point2D nodes3[] = {Point2D {0.0, 0.0},
                      Point2D {2.0, 1.0},
                      Point2D {-1.0, 1.0},
                      Point2D {0.0, 0.0}};
  Bezier cubic_closed(nodes3, 3);
  EXPECT_NEAR(winding_number(Point2D({0.0, 0.0}), cubic_closed, edge_tol, EPS),
              0.301208191175,
              abs_tol);
}

TEST(primal_winding_number, degenerate_cases)
{
  using Point2D = primal::Point<double, 2>;
  using Bezier = primal::BezierCurve<double, 2>;

  double abs_tol = 1e-4;
  double edge_tol = 1e-8;
  double EPS = 1e-50;

  // Backtracking quadratic curve
  Point2D bt_nodes[] = {Point2D {0.0, 0.0},
                        Point2D {0.0, 1.0},
                        Point2D {0.0, 0.5}};
  Bezier bt_curve(bt_nodes, 2);
  EXPECT_NEAR(winding_number(Point2D({0.0, -0.25}), bt_curve, edge_tol, EPS),
              0,
              abs_tol);
  EXPECT_NEAR(winding_number(Point2D({0.0, 0.25}), bt_curve, edge_tol, EPS),
              0,
              abs_tol);
  EXPECT_NEAR(winding_number(Point2D({0.0, 0.75}), bt_curve, edge_tol, EPS),
              0,
              abs_tol);
  EXPECT_NEAR(winding_number(Point2D({0.0, 1.25}), bt_curve, edge_tol, EPS),
              0,
              abs_tol);
  EXPECT_NEAR(winding_number(Point2D({0.25, 0.25}), bt_curve, edge_tol, EPS),
              -0.25,
              abs_tol);

  // empty curve, high order. Gets caught by isLinear
  Point2D empty_nodes[] = {Point2D {0.0, 0.0},
                           Point2D {0.0, 0.0},
                           Point2D {0.0, 0.0},
                           Point2D {0.0, 0.0}};
  Bezier empty_curve(empty_nodes, 3);
  axom::Array<Point2D> test_points = {Point2D {0.0, 1.0},
                                      Point2D {1.0, 0.0},
                                      Point2D {1.0, 1.0},
                                      Point2D {0.0, 1.0}};

  EXPECT_NEAR(winding_number(Point2D({0, 0}), empty_curve, edge_tol, EPS),
              0,
              abs_tol);

  for(auto pt : test_points)
  {
    EXPECT_NEAR(winding_number(pt, empty_curve, edge_tol, EPS), 0, abs_tol);
    pt[0] *= -1;
    EXPECT_NEAR(winding_number(pt, empty_curve, edge_tol, EPS), 0, abs_tol);
    pt[1] *= -1;
    EXPECT_NEAR(winding_number(pt, empty_curve, edge_tol, EPS), 0, abs_tol);
    pt[0] *= -1;
    EXPECT_NEAR(winding_number(pt, empty_curve, edge_tol, EPS), 0, abs_tol);
  }

  // Cubic curve
  Point2D cubic_nodes[] = {Point2D {0.0, 0.0},
                           Point2D {0.0, 0.0},
                           Point2D {0.0, 0.0},
                           Point2D {0.0, 1.0},
                           Point2D {-1.0, 1.0},
                           Point2D {-1.0, 0.0},
                           Point2D {-1.0, 0.0},
                           Point2D {-1.0, 0.0}};
  Bezier cubic(cubic_nodes, 7);

  EXPECT_NEAR(  // Query on initial endpoint of cubic
    winding_number(Point2D({-1.0, 0.0}), cubic, edge_tol, EPS),
    0.25,
    abs_tol);
  EXPECT_NEAR(  // Query on terminal endpoint of cubic
    winding_number(Point2D({-1.0, 0.0}), cubic, edge_tol, EPS),
    0.25,
    abs_tol);
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  int result = RUN_ALL_TESTS();

  return result;
}

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

TEST(primal_winding_number, containment_protocol)
{
  // Test that containment procedures are consistent
  using Point2D = primal::Point<double, 2>;
  using Bezier = primal::BezierCurve<double, 2>;
  using CPolygon = primal::CurvedPolygon<double, 2>;

  double edge_tol = 1e-8;
  double EPS = primal::PRIMAL_TINY;
  int npts = 15;

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

  for(auto alg_type : {primal::RECURSIVE_BISECTION, primal::CONVEX_QUADRATURE})
  {
    // Inner loop is considered "interior" with nonzero protocol. Default behavior.
    bool useNonzeroRule = true;
    EXPECT_TRUE(in_curved_polygon(Point2D({0.5, 0.5}),
                                  loop_poly,
                                  useNonzeroRule,
                                  edge_tol,
                                  EPS,
                                  alg_type,
                                  npts));

    // Inner loop is considered "exterior" with even/odd protocol
    useNonzeroRule = false;
    EXPECT_FALSE(in_curved_polygon(Point2D({0.5, 0.5}),
                                   loop_poly,
                                   useNonzeroRule,
                                   edge_tol,
                                   EPS,
                                   alg_type,
                                   npts));
  }
}

TEST(primal_winding_number, simple_cases)
{
  // Test points that are straightforwardly "inside" or "outside"
  //  the closed shape
  using Point2D = primal::Point<double, 2>;
  using Vector2D = primal::Vector<double, 2>;
  using Bezier = primal::BezierCurve<double, 2>;
  using CPolygon = primal::CurvedPolygon<double, 2>;

  double abs_tol = 1e-8;
  double EPS = primal::PRIMAL_TINY;
  int npts = 15;

  // Simple closed shape with cubic edges
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

  for(auto alg_type : {primal::RECURSIVE_BISECTION})
  {
    for(int i = 3; i < 4; i++)
    {
      double edge_tol = std::pow(10, -i);

      for(double t = 2 * edge_tol; t < (1.0 - 2 * edge_tol); t += 0.1)
      {
        Point2D on_curve(top_curve.evaluate(t));
        Vector2D tangent(top_curve.dt(t));

        // Check interior points
        Point2D q_int({
          -tangent[1] * 2.0 * edge_tol / tangent.norm() + on_curve[0],
          tangent[0] * 2.0 * edge_tol / tangent.norm() + on_curve[1],
        });

        EXPECT_NEAR(
          winding_number(q_int, simple_shape, edge_tol, EPS, alg_type, npts),
          1.0,
          abs_tol);

        // Check exterior points
        Point2D q_ext({
          -tangent[1] * 2.0 * -edge_tol / tangent.norm() + on_curve[0],
          tangent[0] * 2.0 * -edge_tol / tangent.norm() + on_curve[1],
        });

        EXPECT_NEAR(
          winding_number(q_ext, simple_shape, edge_tol, EPS, alg_type, npts),
          0.0,
          abs_tol);
      }
    }
  }

  return;
}

TEST(primal_winding_number, closure_edge_cases)
{
  // Tests for when query is on the linear closure
  using Point2D = primal::Point<double, 2>;
  using Bezier = primal::BezierCurve<double, 2>;

  double abs_tol = 1e-8;
  double edge_tol = 1e-8;
  double EPS = primal::PRIMAL_TINY;
  int npts = 15;

  // Test on linear cases
  Point2D linear_nodes[] = {Point2D {0.0, 0.0}, Point2D {1.0, 1.0}};
  Bezier linear(linear_nodes, 1);

  for(auto alg_type : {primal::RECURSIVE_BISECTION})
  {
    EXPECT_NEAR(
      winding_number(Point2D({-0.45, -0.45}), linear, edge_tol, EPS, alg_type, npts),
      0.0,
      abs_tol);
    EXPECT_NEAR(
      winding_number(Point2D({1.45, 1.45}), linear, edge_tol, EPS, alg_type, npts),
      0.0,
      abs_tol);
  }

  // Extra tests if initial and terminal tangent lines are collinear
  Point2D quartic_nodes[] = {Point2D {0.1, 0.0},
                             Point2D {1.0, 0.0},
                             Point2D {0.0, 1.0},
                             Point2D {-1.0, 0.0},
                             Point2D {-0.1, 0.0}};
  Bezier quartic(quartic_nodes, 4);

  for(auto alg_type : {primal::RECURSIVE_BISECTION})
  {
    // Tangent lines in opposite directions
    EXPECT_NEAR(
      winding_number(Point2D({0, 0}), quartic, edge_tol, EPS, alg_type, npts),
      0.5,
      abs_tol);
    EXPECT_NEAR(
      winding_number(Point2D({2.5, 0}), quartic, edge_tol, EPS, alg_type, npts),
      0.0,
      abs_tol);
    EXPECT_NEAR(
      winding_number(Point2D({-2.5, 0}), quartic, edge_tol, EPS, alg_type, npts),
      0.0,
      abs_tol);
  }

  // Flip the curve vertically
  quartic[2] = Point2D({0.0, -1.0});
  for(auto alg_type : {primal::RECURSIVE_BISECTION})
  {
    EXPECT_NEAR(
      winding_number(Point2D({0, 0}), quartic, edge_tol, EPS, alg_type, npts),
      -0.5,
      abs_tol);
    EXPECT_NEAR(
      winding_number(Point2D({2.5, 0}), quartic, edge_tol, EPS, alg_type, npts),
      0.0,
      abs_tol);
    EXPECT_NEAR(
      winding_number(Point2D({-2.5, 0}), quartic, edge_tol, EPS, alg_type, npts),
      0.0,
      abs_tol);
  }

  // Flip one of the tangent lines
  quartic[1] = Point2D({0.0, 0.0});
  for(auto alg_type : {primal::RECURSIVE_BISECTION})
  {
    EXPECT_NEAR(
      winding_number(Point2D({0, 0}), quartic, edge_tol, EPS, alg_type, npts),
      -0.5,
      abs_tol);
    EXPECT_NEAR(
      winding_number(Point2D({2.5, 0}), quartic, edge_tol, EPS, alg_type, npts),
      0.0,
      abs_tol);
    EXPECT_NEAR(
      winding_number(Point2D({-2.5, 0}), quartic, edge_tol, EPS, alg_type, npts),
      0.0,
      abs_tol);
  }

  // Flip vertically again
  quartic[2] = Point2D({0.0, 1.0});
  for(auto alg_type : {primal::RECURSIVE_BISECTION})
  {
    EXPECT_NEAR(
      winding_number(Point2D({0, 0}), quartic, edge_tol, EPS, alg_type, npts),
      0.5,
      abs_tol);
    EXPECT_NEAR(
      winding_number(Point2D({2.5, 0}), quartic, edge_tol, EPS, alg_type, npts),
      0.0,
      abs_tol);
    EXPECT_NEAR(
      winding_number(Point2D({-2.5, 0}), quartic, edge_tol, EPS, alg_type, npts),
      0.0,
      abs_tol);
  }
}

TEST(primal_winding_number, corner_cases)
{
  // Tests for when query is identically on either endpoint of the Bezier curve.
  //  Conventionally undefined mathematically, we return the limiting value,
  //  which depends on tangent lines at the query point
  using Point2D = primal::Point<double, 2>;
  using Bezier = primal::BezierCurve<double, 2>;

  double abs_tol = 1e-4;
  double edge_tol = 1e-8;
  double EPS = primal::PRIMAL_TINY;

  // Line segment
  Point2D linear_nodes[] = {Point2D {0.0, 0.0}, Point2D {1.0, 1.0}};
  Bezier linear(linear_nodes, 1);

  // Cubic curve
  Point2D cubic_nodes[] = {Point2D {0.0, 0.0},
                           Point2D {0.0, 1.0},
                           Point2D {-1.0, 1.0},
                           Point2D {-1.0, 0.0}};
  Bezier cubic(cubic_nodes, 3);

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

  // Query point on both endpoints
  Point2D closed_cubic_nodes[] = {Point2D {0.0, 0.0},
                                  Point2D {2.0, 1.0},
                                  Point2D {-1.0, 1.0},
                                  Point2D {0.0, 0.0}};
  Bezier cubic_closed(closed_cubic_nodes, 3);
  EXPECT_NEAR(winding_number(Point2D({0.0, 0.0}), cubic_closed, edge_tol, EPS),
              0.301208191175,
              abs_tol);

  // Extra tests if initial and terminal tangent lines are collinear
  Point2D quartic_nodes[] = {Point2D {0.1, 0.0},
                             Point2D {1.0, 0.0},
                             Point2D {0.0, 1.0},
                             Point2D {-1.0, 0.0},
                             Point2D {-0.1, 0.0}};
  Bezier quartic(quartic_nodes, 4);

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

TEST(primal_winding_number, edge_cases)
{
  // Tests for when query is identically on interior of Bezier curve.
  //   Conventionally undefined mathematically. Uses endpoint formulas
  //   to determine value after some number of bisections
  using Point2D = primal::Point<double, 2>;
  using Bezier = primal::BezierCurve<double, 2>;

  double abs_tol = 1e-4;
  double edge_tol = 1e-8;
  double EPS = primal::PRIMAL_TINY;

  // Line segment
  Point2D linear_nodes[] = {Point2D {0.0, 0.0}, Point2D {1.0, 1.0}};
  Bezier linear(linear_nodes, 1);

  // At any point on a line, returns 0
  for(double t = 0.1; t < 1; t += 0.1)
    EXPECT_NEAR(winding_number(Point2D({t, t}), linear, edge_tol, EPS),
                0.0,
                abs_tol);

  // Cubic curve, where query is not on an endpoint after any number of bisections
  Point2D cubic_nodes[] = {Point2D {0.0, 0.0},
                           Point2D {0.0, 1.0},
                           Point2D {-1.0, 1.0},
                           Point2D {-1.0, 0.0}};
  Bezier cubic(cubic_nodes, 3);

  EXPECT_NEAR(winding_number(cubic.evaluate(0.1), cubic, edge_tol, EPS),
              0.276676361896,
              abs_tol);
  EXPECT_NEAR(winding_number(cubic.evaluate(0.4), cubic, edge_tol, EPS),
              0.310998033871,
              abs_tol);
  EXPECT_NEAR(winding_number(cubic.evaluate(0.7), cubic, edge_tol, EPS),
              0.305165888012,
              abs_tol);

  // Cubic curve with internal loop
  Point2D cubic_loop_nodes[] = {Point2D {0.0, 0.0},
                                Point2D {2.0, 1.0},
                                Point2D {-1.0, 1.0},
                                Point2D {1.0, 0.0}};
  Bezier cubic_loop(cubic_loop_nodes, 3);
  EXPECT_NEAR(winding_number(Point2D({0.5, 0.3}), cubic_loop, edge_tol, EPS),
              0.327979130377,
              abs_tol);
  EXPECT_NEAR(winding_number(Point2D({0.5, 0.75}), cubic_loop, edge_tol, EPS),
              0.687167046798,
              abs_tol);
}

TEST(primal_winding_number, degenerate_cases)
{
  // Tests for when Bezier curves are defined with duplicate nodes
  using Point2D = primal::Point<double, 2>;
  using Bezier = primal::BezierCurve<double, 2>;

  double abs_tol = 1e-8;
  double edge_tol = 1e-8;
  double EPS = primal::PRIMAL_TINY;

  // Flat curve with anti-parallel tangent lines
  Point2D double_bt_nodes[] = {Point2D {0.0, 1.0},
                               Point2D {0.0, 2.0},
                               Point2D {0.0, 0.0},
                               Point2D {0.0, -2.0},
                               Point2D {0.0, -1.0}};
  Bezier double_bt(double_bt_nodes, 4);

  Point2D linear_nodes[] = {Point2D {0.0, 1.0}, Point2D {0.0, -1.0}};
  Bezier linear(linear_nodes, 1);

  for(double t = -3.0; t <= 3.0; t += 0.1)
  {
    EXPECT_NEAR(winding_number(Point2D({0.0, t}), double_bt, edge_tol, EPS),
                0.0,
                abs_tol);
    EXPECT_NEAR(winding_number(Point2D({1.0, t}), double_bt, edge_tol, EPS),
                winding_number(Point2D({1.0, t}), linear, edge_tol, EPS),
                abs_tol);
  }

  // Check endpoints specifically
  EXPECT_NEAR(winding_number(Point2D({0.0, 1.0}), double_bt, edge_tol, EPS),
              0.0,
              abs_tol);
  EXPECT_NEAR(winding_number(Point2D({0.0, -1.0}), double_bt, edge_tol, EPS),
              0.0,
              abs_tol);

  // empty curve, high order.
  Point2D empty_nodes[] = {Point2D {0.0, 0.0},
                           Point2D {0.0, 0.0},
                           Point2D {0.0, 0.0},
                           Point2D {0.0, 0.0}};
  Bezier empty_curve(empty_nodes, 3);
  axom::Array<Point2D> test_points = {Point2D {0.0, 0.0},
                                      Point2D {0.0, 1.0},
                                      Point2D {1.0, 0.0},
                                      Point2D {1.0, 1.0},
                                      Point2D {0.0, 1.0}};

  for(auto pt : test_points)
    EXPECT_NEAR(winding_number(pt, empty_curve, edge_tol, EPS), 0, abs_tol);

  // Check default empty Bezier curves
  Bezier very_empty_curve(-1);
  for(auto pt : test_points)
    EXPECT_NEAR(winding_number(pt, very_empty_curve, edge_tol, EPS), 0, abs_tol);

  very_empty_curve.setOrder(0);
  for(auto pt : test_points)
    EXPECT_NEAR(winding_number(pt, very_empty_curve, edge_tol, EPS), 0, abs_tol);

  // Cubic curve with many duplicated endpoints
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

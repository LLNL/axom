// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/primal.hpp"
#include "axom/slic.hpp"
#include "axom/fmt.hpp"

#include "gtest/gtest.h"

namespace primal = axom::primal;

TEST(primal_integral, evaluate_area_integral)
{
  using Point2D = primal::Point<double, 2>;
  using Bezier = primal::BezierCurve<double, 2>;
  double abs_tol = 1e-10;

  // Quadrature nodes. Should be sufficiently high to pass tests
  int npts = 20;

  // Define anonymous functions for testing
  auto const_integrand = [](Point2D /*x*/) -> double { return 1.0; };
  auto poly_integrand = [](Point2D x) -> double { return x[0] * x[1] * x[1]; };
  auto transc_integrand = [](Point2D x) -> double {
    return std::sin(x[0] * x[1]);
  };

  // Test on triangular domain
  Point2D trinodes1[] = {Point2D {0.0, 0.0}, Point2D {1.0, 0.0}};
  Bezier tri1(trinodes1, 1);

  Point2D trinodes2[] = {Point2D {1.0, 0.0}, Point2D {0.0, 1.0}};
  Bezier tri2(trinodes2, 1);

  Point2D trinodes3[] = {Point2D {0.0, 1.0}, Point2D {0.0, 0.0}};
  Bezier tri3(trinodes3, 1);

  axom::Array<Bezier> triangle({tri1, tri2, tri3});

  // Compare against hand computed/high-precision calculated values
  EXPECT_NEAR(evaluate_area_integral(triangle, const_integrand, npts),
              0.5,
              abs_tol);
  EXPECT_NEAR(evaluate_area_integral(triangle, poly_integrand, npts),
              1.0 / 60.0,
              abs_tol);
  EXPECT_NEAR(evaluate_area_integral(triangle, transc_integrand, npts),
              0.0415181074232,
              abs_tol);

  // Test on parabolic domain (between f(x) = 1-x^2 and g(x) = x^2-1, shifted to the right 1 unit)
  Point2D paranodes1[] = {Point2D {2.0, 0.0},
                          Point2D {1.0, 2.0},
                          Point2D {0.0, 0.0}};
  Bezier para1(paranodes1, 2);

  Point2D paranodes2[] = {Point2D {0.0, 0.0},
                          Point2D {1.0, -2.0},
                          Point2D {2.0, 0.0}};
  Bezier para2(paranodes2, 2);

  axom::Array<Bezier> parabola_shape({para1, para2});

  // Compare against hand computed/high-precision calculated values
  EXPECT_NEAR(evaluate_area_integral(parabola_shape, const_integrand, npts),
              8.0 / 3.0,
              abs_tol);
  EXPECT_NEAR(evaluate_area_integral(parabola_shape, poly_integrand, npts),
              64.0 / 105.0,
              abs_tol);
  EXPECT_NEAR(evaluate_area_integral(parabola_shape, transc_integrand, npts),
              0.0,
              abs_tol);

  // Ensure compatibility with curved polygons
  Bezier pedges[2] = {para1, para2};
  primal::CurvedPolygon<double, 2> parabola_polygon(pedges, 2);
  EXPECT_NEAR(evaluate_area_integral(parabola_polygon, const_integrand, npts),
              8.0 / 3.0,
              abs_tol);
  EXPECT_NEAR(evaluate_area_integral(parabola_polygon, poly_integrand, npts),
              64.0 / 105.0,
              abs_tol);
  EXPECT_NEAR(evaluate_area_integral(parabola_polygon, transc_integrand, npts),
              0.0,
              abs_tol);
}

TEST(primal_integral, evaluate_line_integral_scalar)
{
  using Point2D = primal::Point<double, 2>;
  using Bezier = primal::BezierCurve<double, 2>;
  double abs_tol = 1e-10;

  // Quadrature nodes. Should be sufficiently high to pass tests
  int npts = 30;

  // Define anonymous functions for testing
  auto const_integrand = [](Point2D /*x*/) -> double { return 1.0; };
  auto poly_integrand = [](Point2D x) -> double { return x[0] * x[1] * x[1]; };
  auto transc_integrand = [](Point2D x) -> double {
    return std::sin(x[0] * x[1]);
  };

  // Test on single parabolic segment
  Point2D paranodes[] = {Point2D {-1.0, 1.0},
                         Point2D {0.5, -2.0},
                         Point2D {2.0, 4.0}};
  Bezier parabola_segment(paranodes, 2);

  // Compare against hand computed/high-precision calculated values.

  // Constant integrand line integral is equivalent to arc-length calculation
  EXPECT_NEAR(evaluate_line_integral(parabola_segment, const_integrand, npts),
              6.12572661998,
              abs_tol);

  EXPECT_NEAR(evaluate_line_integral(parabola_segment, poly_integrand, npts),
              37.8010703669,
              abs_tol);
  EXPECT_NEAR(evaluate_line_integral(parabola_segment, transc_integrand, npts),
              0.495907795678,
              abs_tol);

  // Test on a collection of Bezier curves
  Point2D segnodes1[] = {Point2D {-1.0, -1.0},
                         Point2D {-1.0 / 3.0, 1.0},
                         Point2D {1.0 / 3.0, -1.0},
                         Point2D {1.0, 1.0}};
  Bezier cubic_segment(segnodes1, 3);

  Point2D segnodes2[] = {Point2D {1.0, 1.0}, Point2D {-1.0, 0.0}};
  Bezier linear_segment(segnodes2, 1);

  Point2D segnodes3[] = {Point2D {-1.0, 0.0},
                         Point2D {-3.0, 1.0},
                         Point2D {-1.0, 2.0}};
  Bezier quadratic_segment(segnodes3, 2);

  axom::Array<Bezier> connected_curve(
    {cubic_segment, linear_segment, quadratic_segment});
  EXPECT_NEAR(evaluate_line_integral(connected_curve, const_integrand, npts),
              8.28968500196,
              abs_tol);
  EXPECT_NEAR(evaluate_line_integral(connected_curve, poly_integrand, npts),
              -5.97565740064,
              abs_tol);
  EXPECT_NEAR(evaluate_line_integral(connected_curve, transc_integrand, npts),
              -0.574992518405,
              abs_tol);

  // Test algorithm on disconnected curves
  axom::Array<Bezier> disconnected_curve({cubic_segment, quadratic_segment});
  EXPECT_NEAR(evaluate_line_integral(disconnected_curve, const_integrand, npts),
              6.05361702446,
              abs_tol);
  EXPECT_NEAR(evaluate_line_integral(disconnected_curve, poly_integrand, npts),
              -6.34833539689,
              abs_tol);
  EXPECT_NEAR(evaluate_line_integral(disconnected_curve, transc_integrand, npts),
              -0.914161242161,
              abs_tol);
}

TEST(primal_integral, evaluate_line_integral_vector)
{
  using Point2D = primal::Point<double, 2>;
  using Vector2D = primal::Vector<double, 2>;
  using Bezier = primal::BezierCurve<double, 2>;
  double abs_tol = 1e-10;

  // Quadrature nodes. Should be sufficiently high to pass tests
  int npts = 30;

  // Test on a single line segment
  auto vec_field = [](Point2D x) -> Vector2D {
    return Vector2D({x[1] * x[1], 3 * x[0] - 6 * x[1]});
  };

  Point2D segnodes[] = {Point2D {3.0, 7.0}, Point2D {0.0, 12.0}};
  Bezier linear_segment(segnodes, 1);

  // Compare against hand computed values
  EXPECT_NEAR(evaluate_line_integral(linear_segment, vec_field, npts),
              -1079.0 / 2.0,
              abs_tol);

  // Test on a closed curve
  auto area_field = [](Point2D x) -> Vector2D {
    return Vector2D({-0.5 * x[1], 0.5 * x[0]});
  };
  auto conservative_field = [](Point2D x) -> Vector2D {
    return Vector2D({2 * x[0] * x[1] * x[1], 2 * x[0] * x[0] * x[1]});
  };
  auto winding_field = [](Point2D x) -> Vector2D {
    double denom = 2 * M_PI * (x[0] * x[0] + x[1] * x[1]);
    return Vector2D({-x[1] / denom, x[0] / denom});
  };

  Point2D paranodes1[] = {Point2D {1.0, 0.0},
                          Point2D {0.0, 2.0},
                          Point2D {-1.0, 0.0}};
  Bezier para1(paranodes1, 2);

  Point2D paranodes2[] = {Point2D {-1.0, 0.0},
                          Point2D {0.0, -2.0},
                          Point2D {1.0, 0.0}};
  Bezier para2(paranodes2, 2);

  axom::Array<Bezier> parabola_shape({para1, para2});

  // This vector field calculates the area of the region
  EXPECT_NEAR(evaluate_line_integral(parabola_shape, area_field, npts),
              8.0 / 3.0,
              abs_tol);

  // This vector field is conservative, so it should evaluate to zero
  EXPECT_NEAR(evaluate_line_integral(parabola_shape, conservative_field, npts),
              0.0,
              abs_tol);

  // This vector field is generated by a in/out query, should return 1 (inside)
  EXPECT_NEAR(evaluate_line_integral(parabola_shape, winding_field, npts),
              1.0,
              abs_tol);

  // Test algorithm on disconnected curves
  Point2D paranodes2_shifted[] = {Point2D {-1.0, -1.0},
                                  Point2D {0.0, -3.0},
                                  Point2D {1.0, -1.0}};
  Bezier para2_shift(paranodes2_shifted, 2);
  axom::Array<Bezier> disconnected_parabola_shape({para1, para2_shift});

  EXPECT_NEAR(
    evaluate_line_integral(disconnected_parabola_shape, area_field, npts),
    11.0 / 3.0,
    abs_tol);
  EXPECT_NEAR(
    evaluate_line_integral(disconnected_parabola_shape, conservative_field, npts),
    0.0,
    abs_tol);
  EXPECT_NEAR(
    evaluate_line_integral(disconnected_parabola_shape, winding_field, npts),
    0.75,
    abs_tol);
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  int result = RUN_ALL_TESTS();

  return result;
}

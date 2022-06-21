// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/primal.hpp"
#include "axom/slic.hpp"
#include "axom/fmt.hpp"

#include "gtest/gtest.h"

#include <iostream>

namespace primal = axom::primal;
namespace slic = axom::slic;

/*!
 * \brief Utility function to initialize the logger
 */
void initializeLogger()
{
  // Initialize Logger
  slic::initialize();
  slic::setLoggingMsgLevel(axom::slic::message::Info);

  slic::LogStream* logStream;

#ifdef AXOM_USE_MPI
  std::string fmt = "[<RANK>][<LEVEL>]: <MESSAGE>\n";
  #ifdef AXOM_USE_LUMBERJACK
  const int RLIMIT = 8;
  logStream = new slic::LumberjackStream(&std::cout, MPI_COMM_WORLD, RLIMIT, fmt);
  #else
  logStream = new slic::SynchronizedStream(&std::cout, MPI_COMM_WORLD, fmt);
  #endif
#else
  std::string fmt = "[<LEVEL>]: <MESSAGE>\n";
  logStream = new slic::GenericOutputStream(&std::cout, fmt);
#endif  // AXOM_USE_MPI

  slic::addStreamToAllMsgLevels(logStream);
}

/*!
 * \brief Utility function to finalize the logger
 */
void finalizeLogger()
{
  if(slic::isInitialized())
  {
    slic::flushStreams();
    slic::finalize();
  }
}

TEST(primal_integral, evaluate_area_integral)
{
  using Point2D = primal::Point<double, 2>;
  using Bezier = primal::BezierCurve<double, 2>;
  double abs_tol = 1e-10;

  // Quadrature nodes. Should be sufficiently high to pass tests
  int npts = 15;

  // Define anonymous functions for testing
  auto const_integrand = [](Point2D x) -> double { return 1.0; };
  auto poly_integrand = [](Point2D x) -> double { return x[0] * x[1] * x[1]; };
  auto transc_integrand = [](Point2D x) -> double {
    return std::sin(x[0] * x[1]);
  };

  // Test on triangular domain
  double trinodes1[] = {0.0, 1.0, 0.0, 0.0};
  Bezier tri1(trinodes1, 1);

  double trinodes2[] = {1.0, 0.0, 0.0, 1.0};
  Bezier tri2(trinodes2, 1);

  double trinodes3[] = {0.0, 0.0, 1.0, 0.0};
  Bezier tri3(trinodes3, 1);

  axom::Array<Bezier> triangle({tri1, tri2, tri3});
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
  double paranodes1[] = {2.0, 1.0, 0.0, 0.0, 2.0, 0.0};
  Bezier para1(paranodes1, 2);

  double paranodes2[] = {0.0, 1.0, 2.0, 0.0, -2.0, 0.0};
  Bezier para2(paranodes2, 2);

  axom::Array<Bezier> parabola_shape({para1, para2});
  EXPECT_NEAR(evaluate_area_integral(parabola_shape, const_integrand, npts),
              8.0 / 3.0,
              abs_tol);
  EXPECT_NEAR(evaluate_area_integral(parabola_shape, poly_integrand, npts),
              64.0 / 105.0,
              abs_tol);
  EXPECT_NEAR(evaluate_area_integral(parabola_shape, transc_integrand, npts),
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
  auto const_integrand = [](Point2D x) -> double { return 1.0; };
  auto poly_integrand = [](Point2D x) -> double { return x[0] * x[1] * x[1]; };
  auto transc_integrand = [](Point2D x) -> double {
    return std::sin(x[0] * x[1]);
  };

  // Test on single parabolic segment
  double paranodes[] = {-1.0, 0.5, 2.0, 1.0, -2.0, 4.0};
  Bezier parabola_segment(paranodes, 2);

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
  double segnodes1[] = {-1.0, -1.0 / 3.0, 1.0 / 3.0, 1.0, -1.0, 1.0, -1.0, 1.0};
  Bezier cubic_segment(segnodes1, 3);

  double segnodes2[] = {1.0, -1.0, 1.0, 0.0};
  Bezier linear_segment(segnodes2, 1);

  double segnodes3[] = {-1.0, -3.0, -1.0, 0.0, 1.0, 2.0};
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

  double segnodes[] = {3.0, 0.0, 7.0, 12.0};
  Bezier linear_segment(segnodes, 1);

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

  double paranodes1[] = {1.0, 0.0, -1.0, 0.0, 2.0, 0.0};
  Bezier para1(paranodes1, 2);

  double paranodes2[] = {-1.0, 0.0, 1.0, 0.0, -2.0, 0.0};
  Bezier para2(paranodes2, 2);

  axom::Array<Bezier> parabola_shape({para1, para2});
  EXPECT_NEAR(evaluate_line_integral(parabola_shape, area_field, npts),
              8.0 / 3.0,
              abs_tol);
  EXPECT_NEAR(evaluate_line_integral(parabola_shape, conservative_field, npts),
              0.0,
              abs_tol);
  EXPECT_NEAR(evaluate_line_integral(parabola_shape, winding_field, npts),
              1.0,
              abs_tol);
}

int main(int argc, char* argv[])
{
  // -- Initialize logger
  initializeLogger();

  ::testing::InitGoogleTest(&argc, argv);

  int result = RUN_ALL_TESTS();

  // -- Finalize logger
  finalizeLogger();

  return 0;
}
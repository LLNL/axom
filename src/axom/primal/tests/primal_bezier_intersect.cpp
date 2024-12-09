// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file primal_bezier_intersect.cpp
 * \brief This file tests the Bezier curve intersection routines
 */

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/slic.hpp"

#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/operators/intersect.hpp"

#include "axom/core/numerics/matvecops.hpp"

#include <cmath>
#include <iomanip>
#include <algorithm>
#include <sstream>

namespace primal = axom::primal;

/**
 * Helper function to compute the intersections of two curves and check that
 * their intersection points match our expectations, stored in \a exp_s
 * and \a exp_t. Intersections are computed within tolerance \a eps
 * and our checks use \a test_eps.
 *
 * Param \a shouldPrintIntersections is used for debugging and for generating
 * the initial array of expected intersections.
 */
template <typename CoordType>
void checkIntersections(const primal::BezierCurve<CoordType, 2>& curve1,
                        const primal::BezierCurve<CoordType, 2>& curve2,
                        const axom::Array<CoordType>& exp_s,
                        const axom::Array<CoordType>& exp_t,
                        double eps,
                        double test_eps,
                        bool shouldPrintIntersections = false)
{
  constexpr int DIM = 2;
  using Array = axom::Array<CoordType>;

  // Check validity of input data exp_s and exp_t.
  // They should have the same size and be sorted
  EXPECT_EQ(exp_s.size(), exp_t.size());
  EXPECT_TRUE(std::is_sorted(exp_s.begin(), exp_s.end()));
  EXPECT_TRUE(std::is_sorted(exp_t.begin(), exp_t.end()));

  const int num_exp_intersections = static_cast<int>(exp_s.size());
  const bool exp_intersect = (num_exp_intersections > 0);

  // Intersect the two curves, intersection parameters will be
  // in arrays s and t, for curve1 and curve2, respectively
  Array s, t;
  bool curves_intersect = intersect(curve1, curve2, s, t, eps);
  EXPECT_EQ(exp_intersect, curves_intersect);
  EXPECT_EQ(s.size(), t.size());

  // check that we found the expected number of intersection points
  const int num_actual_intersections = static_cast<int>(s.size());
  EXPECT_EQ(num_exp_intersections, num_actual_intersections);

  // check that the evaluated intersection points are identical
  for(int i = 0; i < num_actual_intersections; ++i)
  {
    auto p1 = curve1.evaluate(s[i]);
    auto p2 = curve2.evaluate(t[i]);

    EXPECT_NEAR(0., primal::squared_distance(p1, p2), test_eps);

    for(int d = 0; d < DIM; ++d)
    {
      EXPECT_NEAR(p1[d], p2[d], test_eps);
    }
  }

  // check that the intersections match our precomputed values
  std::sort(s.begin(), s.end());
  std::sort(t.begin(), t.end());

  if(shouldPrintIntersections)
  {
    std::stringstream sstr;

    sstr << "Intersections for curves: "
         << "\n\t" << curve1 << "\n\t" << curve2;

    sstr << "\ns (" << s.size() << "): ";
    for(auto i = 0; i < s.size(); ++i)
    {
      sstr << std::setprecision(16) << s[i] << ",";
    }

    sstr << "\nt (" << t.size() << "): ";
    for(auto i = 0; i < t.size(); ++i)
    {
      sstr << std::setprecision(16) << t[i] << ",";
    }

    SLIC_INFO(sstr.str());
  }

  for(int i = 0; i < num_actual_intersections; ++i)
  {
    EXPECT_NEAR(exp_s[i], s[i], test_eps);
    EXPECT_NEAR(exp_t[i], t[i], test_eps);

    if(shouldPrintIntersections)
    {
      SLIC_INFO("\t" << i << ": {s:" << s[i] << ", t:" << t[i]
                     << std::setprecision(16) << ", s_actual:" << exp_s[i]
                     << ", t_actual:" << exp_t[i] << "}");
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_bezier_inter, linear_bezier)
{
  static const int DIM = 2;

  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  const int order = 1;

  // case 1 -- two lines that intersect at their midpoint
  {
    SCOPED_TRACE("linear bezier simple");

    PointType data1[order + 1] = {PointType {0.0, 0.0}, PointType {1.0, 1.0}};

    BezierCurveType curve1(data1, order);

    PointType data2[order + 1] = {PointType {0.0, 1.0}, PointType {1.0, 0.0}};
    BezierCurveType curve2(data2, order);

    axom::Array<CoordType> exp_intersections1 = {0.5};
    axom::Array<CoordType> exp_intersections2 = {0.5};

    const double eps = 1E-3;
    checkIntersections(curve1,
                       curve2,
                       exp_intersections1,
                       exp_intersections2,
                       eps,
                       eps);
  }

  // case 2 -- two lines that intersect at their endpoints
  {
    SCOPED_TRACE("linear bezier endpoints");

    PointType data1[order + 1] = {PointType {1.0, 1.0}, PointType {0.0, 0.0}};

    BezierCurveType curve1(data1, order);

    PointType data2[order + 1] = {PointType {1.0, 1.0}, PointType {2.0, 0.0}};
    BezierCurveType curve2(data2, order);

    axom::Array<CoordType> exp_intersections1 = {0.0};
    axom::Array<CoordType> exp_intersections2 = {0.0};

    const double eps = 1E-3;
    checkIntersections(curve1,
                       curve2,
                       exp_intersections1,
                       exp_intersections2,
                       eps,
                       eps);
  }

  // case 3 -- two lines that intersect at an interior point
  {
    SCOPED_TRACE("linear bezier non-midpoint");

    PointType data1[order + 1] = {PointType {0.0, 0.0}, PointType {4.0, 2.0}};

    BezierCurveType curve1(data1, order);

    PointType data2[order + 1] = {PointType {-2.0, 2.0}, PointType {2.0, 0.0}};
    BezierCurveType curve2(data2, order);

    axom::Array<CoordType> exp_intersections1 = {.25};
    axom::Array<CoordType> exp_intersections2 = {.75};

    const double eps = 1E-3;
    checkIntersections(curve1,
                       curve2,
                       exp_intersections1,
                       exp_intersections2,
                       eps,
                       eps);
  }
}

TEST(primal_bezier_inter, linear_bezier_interp_params)
{
  constexpr int DIM = 2;

  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  const int order = 1;
  const double eps = 1E-3;

  const int num_i_samples = 37;
  const int num_j_samples = 13;
  // NOTE: Skipping endpoints for now.
  for(int i = 0; i < num_i_samples; ++i)
  {
    double t = i / static_cast<double>(num_i_samples + 1);

    for(int j = 0; j < num_j_samples; ++j)
    {
      double s = j / static_cast<double>(num_j_samples + 1);

      std::stringstream sstr;
      sstr << "linear bezier perpendicular (s,t) = (" << s << "," << t << ")";
      SCOPED_TRACE(sstr.str());

      PointType data1[order + 1] = {PointType {0.0, s}, PointType {1.0, s}};
      BezierCurveType curve1(data1, order);

      PointType data2[order + 1] = {PointType {t, 0.0}, PointType {t, 1.0}};
      BezierCurveType curve2(data2, order);

      axom::Array<CoordType> exp_intersections1 = {t};
      axom::Array<CoordType> exp_intersections2 = {s};

      // test for intersections
      checkIntersections(curve1,
                         curve2,
                         exp_intersections1,
                         exp_intersections2,
                         eps,
                         eps);

      // test for intersections after swapping order of curves
      checkIntersections(curve2,
                         curve1,
                         exp_intersections2,
                         exp_intersections1,
                         eps,
                         eps);
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_bezier_inter, no_intersections_bezier)
{
  static const int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  SLIC_INFO("primal: testing bezier intersection");
  SCOPED_TRACE("no intersections");

  const int order = 3;

  // cubic line
  PointType data1[order + 1] = {PointType {0.0, 0.0},
                                PointType {1.0, 0.0},
                                PointType {2.0, 0.0},
                                PointType {3.0, 0.0}};

  BezierCurveType curve1(data1, order);

  // Cubic curve
  PointType data2[order + 1] = {PointType {0.0, 0.5},
                                PointType {1.0, 1.0},
                                PointType {2.0, 3.0},
                                PointType {3.0, 1.5}};
  BezierCurveType curve2(data2, order);

  axom::Array<CoordType> exp_intersections;

  const double eps = 1E-16;
  const double eps_test = 1E-10;

  checkIntersections(curve1,
                     curve2,
                     exp_intersections,
                     exp_intersections,
                     eps,
                     eps_test);
}

//------------------------------------------------------------------------------
TEST(primal_bezier_inter, cubic_quadratic_bezier)
{
  static const int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  SLIC_INFO("primal: testing bezier intersection");
  SLIC_INFO("different orders");

  const int order2 = 3;

  // cubic line
  PointType data1 = PointType::make_point(0.0, 0.0);

  BezierCurveType curve1(0);
  curve1[0] = data1;

  // Cubic curve
  PointType data2[order2 + 1] = {PointType {0.0, 0.5},
                                 PointType {1.0, -1.0},
                                 PointType {2.0, 1.0},
                                 PointType {3.0, -0.5}};
  BezierCurveType curve2(data2, order2);

  // Note: same intersection params for curve and line
  axom::Array<CoordType> exp_intersections = {0.17267316464601146,
                                              0.5,
                                              0.827326835353989};

  const double eps = 1E-16;
  const double eps_test = 1E-10;

  for(int otherorder = 1; otherorder <= 20; ++otherorder)
  {
    curve1.setOrder(otherorder);
    for(int i = 0; i < otherorder; ++i)
    {
      curve1[i][0] = curve1[i][0] * (otherorder - 1) / (1.0 * otherorder);
    }
    curve1[otherorder] = PointType {3.0, 0};
    SLIC_INFO("Testing w/ order 3 and " << otherorder);

    std::stringstream sstr;
    sstr << "different orders study " << otherorder;
    SCOPED_TRACE(sstr.str());

    checkIntersections(curve1,
                       curve2,
                       exp_intersections,
                       exp_intersections,
                       eps,
                       eps_test);
  }
}

//------------------------------------------------------------------------------
TEST(primal_bezier_inter, cubic_bezier_varying_eps)
{
  static const int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  SLIC_INFO("primal: testing bezier intersection");

  const int order = 3;

  // cubic line
  PointType data1[order + 1] = {PointType {0.0, 0.0},
                                PointType {1.0, 0.0},
                                PointType {2.0, 0.0},
                                PointType {3.0, 0.0}};

  BezierCurveType curve1(data1, order);

  // Cubic curve
  PointType data2[order + 1] = {PointType {0.0, 0.5},
                                PointType {1.0, -1.0},
                                PointType {2.0, 1.0},
                                PointType {3.0, -0.5}};
  BezierCurveType curve2(data2, order);

  // Note: same intersection params for curve and line
  axom::Array<CoordType> exp_intersections = {0.17267316464601146,
                                              0.5,
                                              0.827326835353989};

  for(int exp = 1; exp <= 16; ++exp)
  {
    const double eps = std::pow(10, -exp);
    const double eps_test = std::pow(10, -exp + 1);
    SLIC_INFO("Testing w/ eps = " << eps);
    std::stringstream sstr;
    sstr << "cubic eps study " << eps;
    SCOPED_TRACE(sstr.str());

    checkIntersections(curve1,
                       curve2,
                       exp_intersections,
                       exp_intersections,
                       eps,
                       eps_test);
  }
}

TEST(primal_bezier_inter, cubic_bezier_nine_intersections)
{
  static const int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  SLIC_INFO("primal: testing bezier intersection");
  SCOPED_TRACE("cubic bezier -- 9 intersections");

  // A configuration of two cubic bezier curves that intersect at nine points
  const int order = 3;
  PointType data1[order + 1] = {PointType {100, 90},
                                PointType {125, 260},
                                PointType {125, 0},
                                PointType {140, 145}};

  BezierCurveType curve1(data1, order);

  PointType data2[order + 1] = {PointType {75, 110},
                                PointType {265, 120},
                                PointType {0, 130},
                                PointType {145, 135}};
  BezierCurveType curve2(data2, order);

  const double eps = 1E-16;
  const double eps_test = 1E-10;

  axom::Array<CoordType> exp_s = {0.04832125363145223,
                                  0.09691296096235966,
                                  0.1149546049907844,
                                  0.4525395443158645,
                                  0.5130753010911715,
                                  0.5874279469012619,
                                  0.9096761427522992,
                                  0.937347647827529,
                                  0.9747492689591496};

  axom::Array<CoordType> exp_t = {0.05756422605799361,
                                  0.1237935342287382,
                                  0.1676610724548464,
                                  0.4229476589108138,
                                  0.518580927287544,
                                  0.6475905573955213,
                                  0.872202573839462,
                                  0.9370948395920837,
                                  0.9853581324900483};

  checkIntersections(curve1, curve2, exp_s, exp_t, eps, eps_test);
}

/**
 * Helper function to compute the intersections of a curve and a ray and check that
 * their intersection points match our expectations, stored in \a exp_s
 * and \a exp_t. Intersections are computed within tolerance \a eps
 * and our checks use \a test_eps.
 *
 * Param \a shouldPrintIntersections is used for debugging and for generating
 * the initial array of expected intersections.
 */
template <typename CoordType, typename CurveType>
void checkIntersectionsRay(const primal::Ray<CoordType, 2>& ray,
                           const CurveType& curve,
                           const axom::Array<CoordType>& exp_r,
                           const axom::Array<CoordType>& exp_c,
                           double eps,
                           double test_eps,
                           bool shouldPrintIntersections = false)
{
  constexpr int DIM = 2;
  using Array = axom::Array<CoordType>;

  // Check validity of input data exp_c and exp_r.
  // They should have the same size
  EXPECT_EQ(exp_c.size(), exp_r.size());

  const int num_exp_intersections = static_cast<int>(exp_c.size());
  const bool exp_intersect = (num_exp_intersections > 0);

  // Intersect the curve and ray, intersection parameters will be
  // in arrays s and t, for curve and ray, respectively
  Array r, c;
  bool curves_intersect = intersect(ray, curve, r, c, eps);
  EXPECT_EQ(exp_intersect, curves_intersect);
  EXPECT_EQ(r.size(), c.size());

  // check that we found the expected number of intersection points
  const int num_actual_intersections = static_cast<int>(c.size());
  EXPECT_EQ(num_exp_intersections, num_actual_intersections);

  // check that the evaluated intersection points are identical
  for(int i = 0; i < num_actual_intersections; ++i)
  {
    auto p1 = curve.evaluate(c[i]);
    auto p2 = ray.at(r[i]);

    EXPECT_NEAR(0., primal::squared_distance(p1, p2), test_eps);

    for(int d = 0; d < DIM; ++d)
    {
      EXPECT_NEAR(p1[d], p2[d], test_eps);
    }
  }

  if(shouldPrintIntersections)
  {
    std::stringstream sstr;

    sstr << "Intersections for curve and ray: "
         << "\n\t" << curve << "\n\t" << ray;

    sstr << "\ns (" << c.size() << "): ";
    for(auto i = 0u; i < c.size(); ++i)
    {
      sstr << std::setprecision(16) << c[i] << ",";
    }

    sstr << "\nt (" << r.size() << "): ";
    for(auto i = 0u; i < r.size(); ++i)
    {
      sstr << std::setprecision(16) << r[i] << ",";
    }

    SLIC_INFO(sstr.str());
  }

  for(int i = 0; i < num_actual_intersections; ++i)
  {
    EXPECT_NEAR(exp_c[i], c[i], test_eps);
    EXPECT_NEAR(exp_r[i], r[i], test_eps);

    if(shouldPrintIntersections)
    {
      SLIC_INFO("\t" << i << ": {r:" << r[i] << ", c:" << c[i]
                     << std::setprecision(16) << ", s_actual:" << exp_c[i]
                     << ", t_actual:" << exp_r[i] << "}");
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_bezier_inter, ray_linear_bezier)
{
  static const int DIM = 2;

  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using RayType = primal::Ray<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  const int order = 1;

  // case 1 -- Intersect the curve at the midpoint
  {
    SCOPED_TRACE("linear bezier simple");

    PointType ray_origin = PointType::zero();
    VectorType ray_direction({1.0, 1.0});
    RayType ray(ray_origin, ray_direction);

    PointType data[order + 1] = {PointType {1.0, 0.0}, PointType {0.0, 1.0}};
    BezierCurveType curve(data, order);

    axom::Array<CoordType> exp_intersections1 = {std::sqrt(0.5)};
    axom::Array<CoordType> exp_intersections2 = {0.5};

    const double eps = 1E-3;
    checkIntersectionsRay(ray,
                          curve,
                          exp_intersections1,
                          exp_intersections2,
                          eps,
                          eps);
  }

  // case 2 -- Intersect the curve at an endpoint
  {
    SCOPED_TRACE("linear bezier endpoints");

    PointType data[order + 1] = {PointType {1.0, 0.0}, PointType {0.0, 1.0}};
    BezierCurveType curve(data, order);

    // Only count intersections at the t = 0 parameter of the curve
    PointType ray_origin1 = PointType::zero();
    VectorType ray_direction1({1.0, 0.0});
    RayType ray1(ray_origin1, ray_direction1);

    const double eps = 1E-3;
    checkIntersectionsRay(ray1,
                          curve,
                          axom::Array<CoordType>({1.0}),
                          axom::Array<CoordType>({0.0}),
                          eps,
                          eps);

    // Don't count intersections at the t = 1 parameter of the curve
    PointType ray_origin2 = PointType::zero();
    VectorType ray_direction2({0.0, 1.0});
    RayType ray2(ray_origin2, ray_direction2);

    checkIntersectionsRay(ray2,
                          curve,
                          axom::Array<CoordType>(),
                          axom::Array<CoordType>(),
                          eps,
                          eps);

    // Count intersections at the t = 0 parameter of the ray
    PointType ray_origin3({0.5, 0.5});
    VectorType ray_direction3({1.0, 1.0});
    RayType ray3(ray_origin3, ray_direction3);

    checkIntersectionsRay(ray3,
                          curve,
                          axom::Array<CoordType>({0.0}),
                          axom::Array<CoordType>({0.5}),
                          eps,
                          eps);
  }

  // case 3 -- A ray that intersects a curve at an interior point
  {
    SCOPED_TRACE("linear bezier non-midpoint");

    PointType data[order + 1] = {PointType {1.0, 0.0}, PointType {0.0, 1.0}};
    BezierCurveType curve(data, order);

    PointType ray_origin({0.0, 0.0});
    VectorType ray_direction({1.0, 2.0});
    RayType ray(ray_origin, ray_direction);

    axom::Array<CoordType> exp_intersections1 = {std::sqrt(5) / 3.0};
    axom::Array<CoordType> exp_intersections2 = {2.0 / 3.0};

    const double eps = 1E-3;
    checkIntersectionsRay(ray,
                          curve,
                          exp_intersections1,
                          exp_intersections2,
                          eps,
                          eps);
  }
}

//------------------------------------------------------------------------------
TEST(primal_bezier_inter, ray_no_intersections_bezier)
{
  static const int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;
  using RayType = primal::Ray<CoordType, DIM>;

  SLIC_INFO("primal: testing bezier intersection");
  SCOPED_TRACE("no intersections");

  const int order = 3;

  // Ray
  PointType ray_origin({0.0, 0.0});
  VectorType ray_direction({1.0, 0.0});
  RayType ray(ray_origin, ray_direction);

  // Cubic curve
  PointType data[order + 1] = {PointType {0.0, 0.5},
                               PointType {1.0, 1.0},
                               PointType {2.0, 3.0},
                               PointType {3.0, 1.5}};
  BezierCurveType curve(data, order);

  axom::Array<CoordType> exp_intersections;

  const double eps = 1E-16;
  const double eps_test = 1E-10;

  checkIntersectionsRay(ray,
                        curve,
                        exp_intersections,
                        exp_intersections,
                        eps,
                        eps_test);
}

TEST(primal_bezier_inter, ray_linear_bezier_interp_params)
{
  constexpr int DIM = 2;

  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using RayType = primal::Ray<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  const int order = 1;
  const double eps = 1E-3;

  const int num_i_samples = 37;
  const int num_j_samples = 13;

  // NOTE: Skipping endpoints for now.
  for(int i = 0; i < num_i_samples; ++i)
  {
    double t = i / static_cast<double>(num_i_samples + 1);

    for(int j = 0; j < num_j_samples; ++j)
    {
      double s = j / static_cast<double>(num_j_samples + 1);

      std::stringstream sstr;
      sstr << "linear bezier perpendicular (s,t) = (" << s << "," << t << ")";
      SCOPED_TRACE(sstr.str());

      PointType data1[order + 1] = {PointType {0.0, s}, PointType {1.0, s}};
      BezierCurveType curve1(data1, order);

      PointType ray_origin1({t, 0.0});
      VectorType ray_direction1({0.0, 1.0});
      RayType ray1(ray_origin1, ray_direction1);

      axom::Array<CoordType> exp_intersections_t = {t};
      axom::Array<CoordType> exp_intersections_s = {s};

      // test for intersections
      checkIntersectionsRay(ray1,
                            curve1,
                            exp_intersections_s,
                            exp_intersections_t,
                            eps,
                            eps);

      // test for intersections after swapping the curve and ray directions
      PointType data2[order + 1] = {PointType {t, 0.0}, PointType {t, 1.0}};
      BezierCurveType curve2(data2, order);

      PointType ray_origin2({0.0, s});
      VectorType ray_direction2({1.0, 0.0});
      RayType ray2(ray_origin2, ray_direction2);

      checkIntersectionsRay(ray2,
                            curve2,
                            exp_intersections_t,
                            exp_intersections_s,
                            eps,
                            eps);
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_bezier_inter, ray_cubic_quadratic_bezier)
{
  static const int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;
  using RayType = primal::Ray<CoordType, DIM>;

  SLIC_INFO("primal: testing bezier intersection");

  const int order = 3;

  // Ray direction
  VectorType ray_direction({1.0, 0.0});

  // Cubic curve
  PointType data[order + 1] = {PointType {0.0, 0.5},
                               PointType {1.0, -1.0},
                               PointType {2.0, 1.0},
                               PointType {3.0, -0.5}};
  BezierCurveType curve(data, order);

  axom::Array<CoordType> all_intersections = {0.17267316464601146,
                                              0.5,
                                              0.827326835353989};

  const double eps = 1E-16;
  const double eps_test = 1E-10;

  for(CoordType origin = 0.0; origin <= 1.0; origin += 0.05)
  {
    PointType ray_origin({origin, 0.0});
    SLIC_INFO("Testing w/ origin at " << ray_origin);

    RayType ray(ray_origin, ray_direction);

    auto curve_pt_0 = curve.evaluate(all_intersections[0]);
    auto curve_pt_1 = curve.evaluate(all_intersections[1]);
    auto curve_pt_2 = curve.evaluate(all_intersections[2]);

    axom::Array<CoordType> exp_intersections;
    axom::Array<CoordType> ray_intersections;
    if(origin < curve_pt_0[0])
    {
      exp_intersections.push_back(all_intersections[0]);
      ray_intersections.push_back(curve_pt_0[0] - origin);
    }

    if(origin < curve_pt_1[0])
    {
      exp_intersections.push_back(all_intersections[1]);
      ray_intersections.push_back(curve_pt_1[0] - origin);
    }

    if(origin < curve_pt_2[0])
    {
      exp_intersections.push_back(all_intersections[2]);
      ray_intersections.push_back(curve_pt_2[0] - origin);
    }

    checkIntersectionsRay(ray,
                          curve,
                          ray_intersections,
                          exp_intersections,
                          eps,
                          eps_test);
  }
}

//------------------------------------------------------------------------------
TEST(primal_bezier_inter, ray_cubic_bezier_varying_eps)
{
  static const int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;
  using RayType = primal::Ray<CoordType, DIM>;

  SLIC_INFO("primal: testing bezier intersection");

  const int order = 3;

  // Ray
  PointType ray_origin({0.0, 0.0});
  VectorType ray_direction({1.0, 0.0});
  RayType ray(ray_origin, ray_direction);

  // Cubic curve
  PointType data[order + 1] = {PointType {0.0 / 3.0, 0.5},
                               PointType {1.0 / 3.0, -1.0},
                               PointType {2.0 / 3.0, 1.0},
                               PointType {3.0 / 3.0, -0.5}};
  BezierCurveType curve(data, order);

  // Note: same intersection params for curve and line
  axom::Array<CoordType> exp_intersections = {0.17267316464601146,
                                              0.5,
                                              0.827326835353989};

  for(int exp = 1; exp <= 16; ++exp)
  {
    const double eps = std::pow(10, -exp);
    const double eps_test = std::pow(10, -exp + 1);
    SLIC_INFO("Testing w/ eps = " << eps);
    std::stringstream sstr;
    sstr << "cubic eps study " << eps;
    SCOPED_TRACE(sstr.str());

    checkIntersectionsRay(ray,
                          curve,
                          exp_intersections,
                          exp_intersections,
                          eps,
                          eps_test);
  }
}

//------------------------------------------------------------------------------
TEST(primal_bezier_inter, ray_cubic_bezier_four_intersections)
{
  static const int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;
  using RayType = primal::Ray<CoordType, DIM>;

  SLIC_INFO("primal: testing bezier intersection");

  // An intersection of a ray and a high-order Bezier curve,
  //  with one intersection repeated in physical space
  const int order = 7;
  PointType data[order + 1] = {PointType {100, 90},
                               PointType {125, 260},
                               PointType {125, 0},
                               PointType {140, 145},
                               PointType {75, 110},
                               PointType {265, 120},
                               PointType {0, 130},
                               PointType {145, 135}};
  BezierCurveType curve(data, order);

  PointType ray_origin({90.0, 100.0});
  VectorType ray_direction({1.0, 1.0961665896209309});
  RayType ray(ray_origin, ray_direction);

  const double eps = 1E-16;
  const double eps_test = 1E-10;

  axom::Array<CoordType> exp_s = {21.19004780170474,
                                  45.76845689117871,
                                  35.9416068276128,
                                  45.76845689117871};

  axom::Array<CoordType> exp_t = {0.0264232742968,
                                  0.2047732691922508,
                                  0.813490954734,
                                  0.96880275626114684};

  checkIntersectionsRay(ray, curve, exp_s, exp_t, eps, eps_test);
}

//------------------------------------------------------------------------------
TEST(primal_bezier_inter, ray_nurbs_intersections)
{
  static const int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using NURBSCurveType = primal::NURBSCurve<CoordType, DIM>;
  using RayType = primal::Ray<CoordType, DIM>;

  SLIC_INFO("primal: testing NURBS intersection");

  const double eps = 1E-10;
  const double eps_test = 1E-10;

  // NURBS Curve which defines an entire circle
  PointType data[7] = {PointType {1.0, 0.0},
                       PointType {1.0, 2.0},
                       PointType {-1.0, 2.0},
                       PointType {-1.0, 0.0},
                       PointType {-1.0, -2.0},
                       PointType {1.0, -2.0},
                       PointType {1.0, 0.0}};
  double weights[7] = {1.0, 1. / 3., 1. / 3., 1.0, 1. / 3., 1. / 3., 1.0};

  double knots[11] = {-1.0, -1.0, -1.0, -1.0, 0.5, 0.5, 0.5, 2.0, 2.0, 2.0, 2.0};
  NURBSCurveType circle(data, weights, 7, knots, 11);

  // Insert some extra knots to make the Bezier extraction more interesting
  circle.insertKnot(0.0, 1);
  circle.insertKnot(1.0, 2);

  // These parameters include the extra knots at 0.0 and 1.0
  double params[10];
  axom::numerics::linspace(-1.0, 2.0, params, 10);

  PointType ray_origin({0.0, 0.0});
  for(int i = 6; i < 7; ++i)  // Skip the last parameter, which is equal to i=0
  {
    std::cout << "Testing with parameter " << i << " " << params[i] << std::endl;

    VectorType ray_direction(ray_origin, circle.evaluate(params[i]));
    RayType ray(ray_origin, ray_direction);

    checkIntersectionsRay(ray, circle, {1.0}, {params[i]}, eps, eps_test, true);
  }
}

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;  // create & initialize test logger,

  result = RUN_ALL_TESTS();

  return result;
}

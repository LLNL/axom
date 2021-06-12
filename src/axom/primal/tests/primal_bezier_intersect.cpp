// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
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
template <typename CoordType, int DIM>
void checkIntersections(const primal::BezierCurve<CoordType, DIM>& curve1,
                        const primal::BezierCurve<CoordType, DIM>& curve2,
                        const std::vector<CoordType>& exp_s,
                        const std::vector<CoordType>& exp_t,
                        double eps,
                        double test_eps,
                        bool shouldPrintIntersections = false)
{
  using Array = std::vector<CoordType>;

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
    for(auto i = 0u; i < s.size(); ++i)
      sstr << std::setprecision(16) << s[i] << ",";

    sstr << "\nt (" << t.size() << "): ";
    for(auto i = 0u; i < t.size(); ++i)
      sstr << std::setprecision(16) << t[i] << ",";

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

    std::vector<CoordType> exp_intersections1 = {0.5};
    std::vector<CoordType> exp_intersections2 = {0.5};

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

    std::vector<CoordType> exp_intersections1 = {0.0};
    std::vector<CoordType> exp_intersections2 = {0.0};

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

    std::vector<CoordType> exp_intersections1 = {.25};
    std::vector<CoordType> exp_intersections2 = {.75};

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
  static const int DIM = 2;

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

      std::vector<CoordType> exp_intersections1 = {t};
      std::vector<CoordType> exp_intersections2 = {s};

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

  std::vector<CoordType> exp_intersections;

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
  std::vector<CoordType> exp_intersections = {0.17267316464601146,
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
  std::vector<CoordType> exp_intersections = {0.17267316464601146,
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

  std::vector<CoordType> exp_s = {0.04832125363145223,
                                  0.09691296096235966,
                                  0.1149546049907844,
                                  0.4525395443158645,
                                  0.5130753010911715,
                                  0.5874279469012619,
                                  0.9096761427522992,
                                  0.937347647827529,
                                  0.9747492689591496};

  std::vector<CoordType> exp_t = {0.05756422605799361,
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

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;  // create & initialize test logger,

  result = RUN_ALL_TESTS();

  return result;
}

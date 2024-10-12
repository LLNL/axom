// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*! 
 * \file primal_bezier_curve.cpp
 * \brief This file tests primal's Bezier curve functionality
 */

#include "gtest/gtest.h"

#include "axom/slic.hpp"

#include "axom/primal/geometry/NURBSCurve.hpp"
#include "axom/primal/operators/squared_distance.hpp"

namespace primal = axom::primal;

//------------------------------------------------------------------------------
TEST(primal_nurbscurve, evaluate)
{
  SLIC_INFO("Testing NURBS evaluation");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSCurveType = primal::NURBSCurve<CoordType, DIM>;

  const int max_degree = 3;
  PointType data[max_degree + 1] = {PointType {0.6, 1.2, 1.0},
                                    PointType {1.3, 1.6, 1.8},
                                    PointType {2.9, 2.4, 2.3},
                                    PointType {3.2, 3.5, 3.0}};

  // clang-format off
  PointType exp_start_vals[4][4] =  // degree 0                  degree 1                   degree 2                   degree 3
    /* 1 pt */ {{PointType {0.6, 1.2, 1.0}, PointType {-999, -999, -999}, PointType {-999, -999, -999}, PointType {-999, -999, -999}},
    /* 2 pt */  {PointType {0.6, 1.2, 1.0}, PointType {0.6, 1.2, 1.0}, PointType {-999, -999, -999}, PointType {-999, -999, -999}},
    /* 3 pt */  {PointType {0.6, 1.2, 1.0}, PointType {0.6, 1.2, 1.0}, PointType {0.6, 1.2, 1.0}, PointType {-999, -999, -999}},
    /* 4 pt */  {PointType {0.6, 1.2, 1.0}, PointType {0.6, 1.2, 1.0}, PointType {0.6, 1.2, 1.0}, PointType {0.6, 1.2, 1.0}}};

  PointType exp_mid_vals[4][4] =  // degree 0                    degree 1                    degree 2                         degree 3
    /* 1 pt */ {{PointType {0.6, 1.2, 1.0},
                 PointType {-999, -999, -999},
                 PointType {-999, -999, -999},
                 PointType {-999, -999, -999}},
                /* 2 pt */
                {PointType {1.3, 1.6, 1.8},
                 PointType {0.95, 1.4, 1.4},
                 PointType {-999, -999, -999},
                 PointType {-999, -999, -999}},
                /* 3 pt */
                {PointType {1.3, 1.6, 1.8},
                 PointType {1.3, 1.6, 1.8},
                 PointType {1.525, 1.7, 1.725},
                 PointType {-999, -999, -999}},
                /* 4 pt */
                {PointType {2.9, 2.4, 2.3},
                 PointType {2.1, 2.0, 2.05},
                 PointType {2.1, 2.0, 2.05},
                 PointType {2.05, 2.0875, 2.0375}}};

  PointType exp_end_vals[4][4] =  // degree 0                  degree 1                   degree 2                   degree 3
    /* 1 pt */ {{PointType {0.6, 1.2, 1.0},
                 PointType {-999, -999, -999},
                 PointType {-999, -999, -999},
                 PointType {-999, -999, -999}},
                /* 2 pt */
                {PointType {1.3, 1.6, 1.8},
                 PointType {1.3, 1.6, 1.8},
                 PointType {-999, -999, -999},
                 PointType {-999, -999, -999}},
                /* 3 pt */
                {PointType {2.9, 2.4, 2.3},
                 PointType {2.9, 2.4, 2.3},
                 PointType {2.9, 2.4, 2.3},
                 PointType {-999, -999, -999}},
                /* 4 pt */
                {PointType {3.2, 3.5, 3.0},
                 PointType {3.2, 3.5, 3.0},
                 PointType {3.2, 3.5, 3.0},
                 PointType {3.2, 3.5, 3.0}}};
  // clang-format on

  // Test evaluation at various degrees and various number of control
  //  points in `data`.
  for(int npts = 1; npts <= 4; ++npts)
  {
    for(int deg = 0; deg <= npts - 1; ++deg)
    {
      NURBSCurveType curve(data, npts, deg);

      PointType calc_start = curve.evaluate(0.0);
      PointType calc_mid = curve.evaluate(0.5);
      PointType calc_end = curve.evaluate(1.0);

      for(int i = 0; i < DIM; ++i)
      {
        EXPECT_NEAR(calc_start[i], exp_start_vals[npts - 1][deg][i], 1e-15);
        EXPECT_NEAR(calc_mid[i], exp_mid_vals[npts - 1][deg][i], 1e-15);
        EXPECT_NEAR(calc_end[i], exp_end_vals[npts - 1][deg][i], 1e-15);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbscurve, first_derivatives)
{
  SLIC_INFO("Testing NURBS derivative calculation");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using NURBSCurveType = primal::NURBSCurve<CoordType, DIM>;

  const int max_degree = 3;
  PointType data[max_degree + 1] = {PointType {0.6, 1.2, 1.0},
                                    PointType {1.3, 1.6, 1.8},
                                    PointType {2.9, 2.4, 2.3},
                                    PointType {3.2, 3.5, 3.0}};

  // clang-format off
  VectorType exp_start_vals[4][4] =  // degree 0                  degree 1                   degree 2                   degree 3
    /* 1 pt */ {{VectorType {0.0, 0.0, 0.0},
                 VectorType {-999, -999, -999},
                 VectorType {-999, -999, -999},
                 VectorType {-999, -999, -999}},
                /* 2 pt */
                {VectorType {0.0, 0.0, 0.0},
                 VectorType {0.7, 0.4, 0.8},
                 VectorType {-999, -999, -999},
                 VectorType {-999, -999, -999}},
                /* 3 pt */
                {VectorType {0.0, 0.0, 0.0},
                 VectorType {1.4, 0.8, 1.6},
                 VectorType {1.4, 0.8, 1.6},
                 VectorType {-999, -999, -999}},
                /* 4 pt */
                {VectorType {0.0, 0.0, 0.0},
                 VectorType {2.1, 1.2, 2.4},
                 VectorType {2.8, 1.6, 3.2},
                 VectorType {2.1, 1.2, 2.4}}};

  VectorType exp_mid_vals[4][4] =  // degree 0                    degree 1                    degree 2                         degree 3
    /* 1 pt */ {{VectorType {0.0, 0.0, 0.0},
                 VectorType {-999, -999, -999},
                 VectorType {-999, -999, -999},
                 VectorType {-999, -999, -999}},
                /* 2 pt */
                {VectorType {0.0, 0.0, 0.0},
                 VectorType {0.7, 0.4, 0.8},
                 VectorType {-999, -999, -999},
                 VectorType {-999, -999, -999}},
                /* 3 pt */
                {VectorType {0.0, 0.0, 0.0},
                 VectorType {3.2, 1.6, 1.0},
                 VectorType {2.3, 1.2, 1.3},
                 VectorType {-999, -999, -999}},
                /* 4 pt */
                {VectorType {0.0, 0.0, 0.0},
                 VectorType {4.8, 2.4, 1.5},
                 VectorType {3.2, 1.6, 1.0},
                 VectorType {3.15, 2.325, 1.875}}};

  VectorType exp_end_vals[4][4] =  // degree 0                  degree 1                   degree 2                   degree 3
    /* 1 pt */ {{VectorType {0.0, 0.0, 0.0},
                 VectorType {-999, -999, -999},
                 VectorType {-999, -999, -999},
                 VectorType {-999, -999, -999}},
                /* 2 pt */
                {VectorType {0.0, 0.0, 0.0},
                 VectorType {0.7, 0.4, 0.8},
                 VectorType {-999, -999, -999},
                 VectorType {-999, -999, -999}},
                /* 3 pt */
                {VectorType {0.0, 0.0, 0.0},
                 VectorType {3.2, 1.6, 1.0},
                 VectorType {3.2, 1.6, 1.0},
                 VectorType {-999, -999, -999}},
                /* 4 pt */
                {VectorType {0.0, 0.0, 0.0},
                 VectorType {0.9, 3.3, 2.1},
                 VectorType {1.2, 4.4, 2.8},
                 VectorType {0.9, 3.3, 2.1}}};
  // clang-format on

  // Test evaluation at various degrees and various number of control
  //  points in `data`.
  for(int npts = 1; npts <= 4; ++npts)
  {
    for(int deg = 0; deg <= npts - 1; ++deg)
    {
      NURBSCurveType curve(data, npts, deg);

      VectorType calc_start = curve.dt(0.0);
      VectorType calc_mid = curve.dt(0.5);
      VectorType calc_end = curve.dt(1.0);

      for(int i = 0; i < DIM; ++i)
      {
        EXPECT_NEAR(calc_start[i], exp_start_vals[npts - 1][deg][i], 1e-13);
        EXPECT_NEAR(calc_mid[i], exp_mid_vals[npts - 1][deg][i], 1e-13);
        EXPECT_NEAR(calc_end[i], exp_end_vals[npts - 1][deg][i], 1e-13);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbscurve, second_derivatives)
{
  SLIC_INFO("Testing NURBS second derivative calculation");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using NURBSCurveType = primal::NURBSCurve<CoordType, DIM>;

  const int max_degree = 3;
  PointType data[max_degree + 1] = {PointType {0.6, 1.2, 1.0},
                                    PointType {1.3, 1.6, 1.8},
                                    PointType {2.9, 2.4, 2.3},
                                    PointType {3.2, 3.5, 3.0}};

  // clang-format off
  VectorType exp_start_vals[4][4] =  // degree 0                  degree 1                   degree 2                   degree 3
    /* 1 pt */ {{VectorType {0.0, 0.0, 0.0},
                 VectorType {-999, -999, -999},
                 VectorType {-999, -999, -999},
                 VectorType {-999, -999, -999}},
                /* 2 pt */
                {VectorType {0.0, 0.0, 0.0},
                 VectorType {0, 0, 0.0, 0.0},
                 VectorType {-999, -999, -999},
                 VectorType {-999, -999, -999}},
                /* 3 pt */
                {VectorType {0.0, 0.0, 0.0},
                 VectorType {0, 0, 0.0, 0.0},
                 VectorType {1.8, 0.8, -0.6},
                 VectorType {-999, -999, -999}},
                /* 4 pt */
                {VectorType {0.0, 0.0, 0.0},
                 VectorType {0, 0, 0.0, 0.0},
                 VectorType {0.8, 0.0, -4.4},
                 VectorType {5.4, 2.4, -1.8}}};

  VectorType exp_mid_vals[4][4] =  // degree 0                    degree 1                    degree 2                         degree 3
    /* 1 pt */ {{VectorType {0.0, 0.0, 0.0},
                 VectorType {-999, -999, -999},
                 VectorType {-999, -999, -999},
                 VectorType {-999, -999, -999}},
                /* 2 pt */
                {VectorType {0.0, 0.0, 0.0},
                 VectorType {0.0, 0.0, 0.0},
                 VectorType {-999, -999, -999},
                 VectorType {-999, -999, -999}},
                /* 3 pt */
                {VectorType {0.0, 0.0, 0.0},
                 VectorType {0.0, 0.0, 0.0},
                 VectorType {1.8, 0.8, -0.6},
                 VectorType {-999, -999, -999}},
                /* 4 pt */
                {VectorType {0.0, 0.0, 0.0},
                 VectorType {0.0, 0.0, 0.0},
                 VectorType {-4.0, 5.6, 3.6},
                 VectorType {-1.2, 2.1, -0.3}}};

  VectorType exp_end_vals[4][4] =  // degree 0                  degree 1                   degree 2                   degree 3
    /* 1 pt */ {{VectorType {0.0, 0.0, 0.0},
                 VectorType {-999, -999, -999},
                 VectorType {-999, -999, -999},
                 VectorType {-999, -999, -999}},
                /* 2 pt */
                {VectorType {0.0, 0.0, 0.0},
                 VectorType {0.0, 0.0, 0.0},
                 VectorType {-999, -999, -999},
                 VectorType {-999, -999, -999}},
                /* 3 pt */
                {VectorType {0.0, 0.0, 0.0},
                 VectorType {0.0, 0.0, 0.0},
                 VectorType {1.8, 0.8, -0.6},
                 VectorType {-999, -999, -999}},
                /* 4 pt */
                {VectorType {0.0, 0.0, 0.0},
                 VectorType {0.0, 0.0, 0.0},
                 VectorType {-4.0, 5.6, 3.6},
                 VectorType {-7.8, 1.8, 1.2}}};
  // clang-format on

  // Test evaluation at various degrees and various number of control
  //  points in `data`.
  for(int npts = 1; npts <= 4; ++npts)
  {
    for(int deg = 0; deg <= npts - 1; ++deg)
    {
      NURBSCurveType curve(data, npts, deg);

      VectorType calc_start = curve.dtdt(0.0);
      VectorType calc_mid = curve.dtdt(0.5);
      VectorType calc_end = curve.dtdt(1.0);

      for(int i = 0; i < DIM; ++i)
      {
        EXPECT_NEAR(calc_start[i], exp_start_vals[npts - 1][deg][i], 1e-13);
        EXPECT_NEAR(calc_mid[i], exp_mid_vals[npts - 1][deg][i], 1e-13);
        EXPECT_NEAR(calc_end[i], exp_end_vals[npts - 1][deg][i], 1e-13);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbscurve, knot_insertion)
{
  SLIC_INFO("Testing NURBS knot insertion");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSCurveType = primal::NURBSCurve<CoordType, DIM>;

  PointType data[4] = {PointType {0.6, 1.2, 1.0},
                       PointType {1.3, 1.6, 1.8},
                       PointType {2.9, 2.4, 2.3},
                       PointType {3.2, 3.5, 3.0}};

  double weights[4] = {1.0, 2.0, 3.0, 4.0};

  NURBSCurveType curve(data, weights, 4, 3);
  NURBSCurveType curve_knots(data, weights, 4, 3);

  // Knot insertion shouldn't change the paramterization or position
  //  of the curve.

  // Insert a knot at 0.5
  curve_knots.insertKnot(0.5, 3);
  for(double t = 0.0; t <= 1.0; t += 0.1)
  {
    PointType p = curve.evaluate(t);
    PointType p_knots = curve_knots.evaluate(t);
    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_NEAR(p[i], p_knots[i], 1e-13);
    }
  }

  curve_knots.insertKnot(0.7, 1);
  for(double t = 0.0; t <= 1.0; t += 0.1)
  {
    PointType p = curve.evaluate(t);
    PointType p_knots = curve_knots.evaluate(t);
    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_NEAR(p[i], p_knots[i], 1e-13);
    }
  }

  curve_knots.insertKnot(0.4, 2);
  for(double t = 0.0; t <= 1.0; t += 0.1)
  {
    PointType p = curve.evaluate(t);
    PointType p_knots = curve_knots.evaluate(t);
    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_NEAR(p[i], p_knots[i], 1e-13);
    }
  }

  // Inserting knots shouldn't increase its multiplicty
  //  greater than the degree
  const int num_knots = curve_knots.getNumKnots();
  curve_knots.insertKnot(0.5, 3);
  curve_knots.insertKnot(0.0, 1);
  curve_knots.insertKnot(1.0, 1);
  EXPECT_EQ(curve_knots.getNumKnots(), num_knots);

  curve_knots.insertKnot(0.4, 2);
  EXPECT_EQ(curve_knots.getNumKnots(), num_knots + 1);

  for(double t = 0.0; t <= 1.0; t += 0.1)
  {
    PointType p = curve.evaluate(t);
    PointType p_knots = curve_knots.evaluate(t);
    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_NEAR(p[i], p_knots[i], 1e-13);
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbscurve, curve_splitting)
{
  SLIC_INFO("Testing NURBS curve splitting");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSCurveType = primal::NURBSCurve<CoordType, DIM>;

  PointType data[4] = {PointType {0.6, 1.2, 1.0},
                       PointType {1.3, 1.6, 1.8},
                       PointType {2.9, 2.4, 2.3},
                       PointType {3.2, 3.5, 3.0}};

  double weights[4] = {1.0, 2.0, 3.0, 4.0};

  NURBSCurveType curve(data, weights, 4, 3);

  // Do some knot insertion to make it interesting
  curve.insertKnot(0.3, 2);
  curve.insertKnot(0.6, 1);
  curve.insertKnot(0.8, 1);

  NURBSCurveType curve1, curve2, curve3;
  curve.split(0.3, curve1, curve2);
  curve2.split((0.6 - 0.3) / (1 - 0.3), curve2, curve3);

  for(double t = 0.0; t < 1.0; t += 0.05)
  {
    PointType p = curve.evaluate(t);
    
    if( t <= 0.3 )
    {
      PointType p1 = curve1.evaluate( t / 0.3 );
      for(int i = 0; i < DIM; ++i)
      {
        EXPECT_NEAR(p[i], p1[i], 1e-13);
      }
    }

    if( t >= 0.3 && t <= 0.6 )
    {
      PointType p2 = curve2.evaluate( (t - 0.3) / (0.6 - 0.3) );
      for(int i = 0; i < DIM; ++i)
      {
        EXPECT_NEAR(p[i], p2[i], 1e-13);
      }
    }

    if( t >= 0.6 )
    {
      PointType p3 = curve3.evaluate( (t - 0.6) / (1 - 0.6) );
      for(int i = 0; i < DIM; ++i)
      {
        EXPECT_NEAR(p[i], p3[i], 1e-13);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbscurve, bezier_extraction)
{
  return; 

  SLIC_INFO("Testing NURBS Bezier extraction");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSCurveType = primal::NURBSCurve<CoordType, DIM>;

  PointType data[4] = {PointType {0.6, 1.2, 1.0},
                       PointType {1.3, 1.6, 1.8},
                       PointType {2.9, 2.4, 2.3},
                       PointType {3.2, 3.5, 3.0}};

  double weights[4] = {1.0, 2.0, 3.0, 4.0};

  NURBSCurveType curve(data, weights, 4, 3);

  // Do knot insertion, which determines where the Bezier
  //  splitting happens
  curve.insertKnot(0.3, 3);
  curve.insertKnot(0.6, 1);
  curve.insertKnot(0.8, 2);

  NURBSCurveType bezier_curve;
  auto bezier_list = curve.extractBezier();

  EXPECT_EQ(bezier_list.size(), 4);

  for(double t = 0.0; t < 1.0; t += 0.05)
  {
    PointType p = curve.evaluate(t);
    
    if( t <= 0.3 )
    {
      PointType p1 = bezier_list[0].evaluate( t / 0.3 );
      for(int i = 0; i < DIM; ++i)
      {
        EXPECT_NEAR(p[i], p1[i], 1e-13);
      }
    }

    if( t >= 0.3 && t <= 0.6 )
    {
      PointType p2 = bezier_list[1].evaluate( (t - 0.3) / (0.6 - 0.3) );
      for(int i = 0; i < DIM; ++i)
      {
        EXPECT_NEAR(p[i], p2[i], 1e-13);
      }
    }

    if( t >= 0.6 && t <= 0.8 )
    {
      PointType p3 = bezier_list[2].evaluate( (t - 0.6) / (0.8 - 0.6) );
      for(int i = 0; i < DIM; ++i)
      {
        EXPECT_NEAR(p[i], p3[i], 1e-13);
      }
    }

    if( t >= 0.8 )
    {
      PointType p4 = bezier_list[3].evaluate( (t - 0.8) / (1 - 0.8) );
      for(int i = 0; i < DIM; ++i)
      {
        EXPECT_NEAR(p[i], p4[i], 1e-13);
      }
    }
  }
}

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}

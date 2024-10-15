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
TEST(primal_nurbscurve, constructor)
{
  const int DIM = 3;
  using CoordType = double;
  using NURBSCurveType = primal::NURBSCurve<CoordType, DIM>;

  {
    SLIC_INFO("Testing default NURBSCurve constructor");
    NURBSCurveType nCurve;

    int expDegree = -1;
    EXPECT_EQ(expDegree, nCurve.getDegree());
    EXPECT_EQ(expDegree + 1, nCurve.getOrder());
    EXPECT_EQ(expDegree + 1, nCurve.getControlPoints().size());
    EXPECT_EQ(expDegree + 1, nCurve.getKnots().size());
    EXPECT_FALSE(nCurve.isRational());
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbscurve, set_degree)
{
  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSCurveType = primal::NURBSCurve<CoordType, DIM>;

  SLIC_INFO("Test adding control points to empty Bezier curve");

  NURBSCurveType nCurve;
  EXPECT_EQ(-1, nCurve.getDegree());

  const int degree = 1;
  const int npts = 3;
  PointType controlPoints[3] = {PointType {0.6, 1.2, 1.0},
                                PointType {0.0, 1.6, 1.8},
                                PointType {0.2, 1.4, 2.0}};

  nCurve.setNumControlPoints(npts);
  nCurve.setDegree(degree + 1);
  nCurve.setDegree(degree);

  EXPECT_EQ(nCurve.getDegree(), degree);
  EXPECT_EQ(nCurve.getNumControlPoints(), npts);
  EXPECT_EQ(nCurve.getNumKnots(), npts + degree + 1);

  nCurve[0] = controlPoints[0];
  nCurve[1] = controlPoints[1];
  nCurve[2] = controlPoints[2];

  for(int p = 0; p < npts; ++p)
  {
    auto& pt = nCurve[p];
    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_DOUBLE_EQ(controlPoints[p][i], pt[i]);
    }
  }

  nCurve.clear();
  EXPECT_EQ(-1, nCurve.getDegree());
  EXPECT_FALSE(nCurve.isRational());

  nCurve.setParameters(npts, degree);
  nCurve.makeRational();
  EXPECT_TRUE(nCurve.isRational());

  nCurve.setWeight(0, 2.0);
  EXPECT_EQ(nCurve.getWeight(0), 2.0);
}

//------------------------------------------------------------------------------
TEST(primal_nurbscurve, point_array_constructor)
{
  SLIC_INFO("Testing point array constructor");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSCurveType = primal::NURBSCurve<CoordType, DIM>;

  const int npts = 2;
  const int degree = 1;
  PointType controlPoints[npts] = {PointType {0.6, 1.2, 1.0},
                                   PointType {0.0, 1.6, 1.8}};

  double weights[npts] = {1.0, 2.0};

  NURBSCurveType nCurve(controlPoints, npts, degree);
  NURBSCurveType wCurve(controlPoints, weights, npts, degree);

  EXPECT_EQ(1, nCurve.getDegree());
  for(int p = 0; p <= nCurve.getDegree(); ++p)
  {
    auto& pt1 = nCurve[p];
    auto& pt2 = wCurve[p];
    double w = wCurve.getWeight(p);

    EXPECT_DOUBLE_EQ(weights[p], w);
    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_DOUBLE_EQ(controlPoints[p][i], pt1[i]);
      EXPECT_DOUBLE_EQ(controlPoints[p][i], pt2[i]);
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbscurve, axom_array_constructor)
{
  SLIC_INFO("Testing point array constructor");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSCurveType = primal::NURBSCurve<CoordType, DIM>;

  const int npts = 2;
  const int degree = 1;
  axom::Array<PointType> controlPoints {PointType {0.6, 1.2, 1.0},
                                        PointType {0.0, 1.6, 1.8}};

  axom::Array<double> weights {1.0, 2.0};

  NURBSCurveType nCurve(controlPoints, degree);
  NURBSCurveType wCurve(controlPoints, weights, degree);

  EXPECT_EQ(1, nCurve.getDegree());
  for(int p = 0; p <= nCurve.getDegree(); ++p)
  {
    auto& pt1 = nCurve[p];
    auto& pt2 = wCurve[p];
    double w = wCurve.getWeight(p);

    EXPECT_DOUBLE_EQ(weights[p], w);
    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_DOUBLE_EQ(controlPoints[p][i], pt1[i]);
      EXPECT_DOUBLE_EQ(controlPoints[p][i], pt2[i]);
    }
  }
}

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

  double weights[4] = {1.0, 2.0, 3.0, 4.0};

  // clang-format off
  PointType exp_start_vals[4][4] =  // degree 0                  degree 1                   degree 2                   degree 3
           /* 1 pt */ {{PointType {0.6, 1.2, 1.0}, PointType {-999, -999, -999}, PointType {-999, -999, -999}, PointType {-999, -999, -999}},
           /* 2 pt */  {PointType {0.6, 1.2, 1.0}, PointType { 0.6,  1.2,  1.0}, PointType {-999, -999, -999}, PointType {-999, -999, -999}},
           /* 3 pt */  {PointType {0.6, 1.2, 1.0}, PointType { 0.6,  1.2,  1.0}, PointType { 0.6,  1.2,  1.0}, PointType {-999, -999, -999}},
           /* 4 pt */  {PointType {0.6, 1.2, 1.0}, PointType { 0.6,  1.2,  1.0}, PointType { 0.6,  1.2,  1.0}, PointType { 0.6,  1.2,  1.0}}};

  PointType exp_mid_vals[4][4] =  // degree 0                           degree 1                                    degree 2                         degree 3
           /* 1 pt */ {{PointType {0.6, 1.2, 1.0}, PointType {   -999,     -999,     -999}, PointType {  -999,  -999,   -999}, PointType { -999, -999,  -999}},
           /* 2 pt */  {PointType {1.3, 1.6, 1.8}, PointType {16./15.,  22./15.,  23./15.}, PointType {  -999,  -999,   -999}, PointType { -999, -999,  -999}},
           /* 3 pt */  {PointType {1.3, 1.6, 1.8}, PointType {    1.3,      1.6,      1.8}, PointType {1.8125,  1.85, 1.8875}, PointType { -999, -999,  -999}},
           /* 4 pt */  {PointType {2.9, 2.4, 2.3}, PointType {   2.26,     2.08,      2.1}, PointType {  2.26,  2.08,    2.1}, PointType {2.365, 2.32, 2.225}}};

  PointType exp_end_vals[4][4] =     // degree 0                   degree 1                     degree 2                      degree 3
           /* 1 pt */ {{PointType {0.6, 1.2, 1.0}, PointType {-999, -999, -999}, PointType {-999, -999, -999}, PointType {-999, -999, -999}},  
           /* 2 pt */  {PointType {1.3, 1.6, 1.8}, PointType { 1.3,  1.6,  1.8}, PointType {-999, -999, -999}, PointType {-999, -999, -999}},
           /* 3 pt */  {PointType {2.9, 2.4, 2.3}, PointType { 2.9,  2.4,  2.3}, PointType { 2.9,  2.4,  2.3}, PointType {-999, -999, -999}},
           /* 4 pt */  {PointType {3.2, 3.5, 3.0}, PointType { 3.2,  3.5,  3.0}, PointType { 3.2,  3.5,  3.0}, PointType { 3.2,  3.5,  3.0}}};
  // clang-format on

  // Test evaluation at various degrees and various number of control
  //  points in `data`.
  for(int npts = 1; npts <= 4; ++npts)
  {
    for(int deg = 0; deg <= npts - 1; ++deg)
    {
      NURBSCurveType curve(data, weights, npts, deg);

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

  double weights[4] = {1.0, 2.0, 3.0, 4.0};

  // clang-format off
  VectorType exp_start_vals[4][4] =  // degree 0                  degree 1                      degree 2                   degree 3
           /* 1 pt */ {{VectorType {0.0, 0.0, 0.0}, VectorType {-999, -999, -999}, VectorType {-999, -999, -999}, VectorType {-999, -999, -999}},
           /* 2 pt */  {VectorType {0.0, 0.0, 0.0}, VectorType { 1.4,  0.8,  1.6}, VectorType {-999, -999, -999}, VectorType {-999, -999, -999}},
           /* 3 pt */  {VectorType {0.0, 0.0, 0.0}, VectorType { 2.8,  1.6,  3.2}, VectorType { 2.8,  1.6,  3.2}, VectorType {-999, -999, -999}},
           /* 4 pt */  {VectorType {0.0, 0.0, 0.0}, VectorType { 4.2,  2.4,  4.8}, VectorType { 5.6,  3.2,  6.4}, VectorType { 4.2,  2.4,  4.8}}};

  VectorType exp_mid_vals[4][4] =  // degree 0                    degree 1                    degree 2                         degree 3
           /* 1 pt */ {{VectorType {0.0, 0.0, 0.0}, VectorType {   -999,    -999,    -999}, VectorType {  -999,  -999,   -999}, VectorType { -999,  -999, -999}},
           /* 2 pt */  {VectorType {0.0, 0.0, 0.0}, VectorType {28./45., 16./45., 32./45.}, VectorType {  -999,  -999,   -999}, VectorType { -999,  -999, -999}},
           /* 3 pt */  {VectorType {0.0, 0.0, 0.0}, VectorType {    4.8,     2.4,     1.5}, VectorType {2.2375,  1.15, 1.0625}, VectorType { -999,  -999, -999}},
           /* 4 pt */  {VectorType {0.0, 0.0, 0.0}, VectorType {  4.608,   2.304,    1.44}, VectorType { 3.072, 1.536,   0.96}, VectorType {2.652, 2.256, 1.62}}};

  VectorType exp_end_vals[4][4] =  // degree 0                  degree 1                   degree 2                   degree 3
           /* 1 pt */ {{VectorType {0.0, 0.0, 0.0}, VectorType {   -999,    -999,  -999}, VectorType {   -999,    -999,  -999}, VectorType { -999,  -999,  -999}},
           /* 2 pt */  {VectorType {0.0, 0.0, 0.0}, VectorType {   0.35,     0.2,   0.4}, VectorType {   -999,    -999,  -999}, VectorType { -999,  -999,  -999}},
           /* 3 pt */  {VectorType {0.0, 0.0, 0.0}, VectorType {32./15., 16./15., 2./3.}, VectorType {32./15., 16./15., 2./3.}, VectorType { -999,  -999,  -999}},
           /* 4 pt */  {VectorType {0.0, 0.0, 0.0}, VectorType {  0.675,   2.475, 1.575}, VectorType {    0.9,     3.3,   2.1}, VectorType {0.675, 2.475, 1.575}}};
  // clang-format on

  // Test evaluation at various degrees and various number of control
  //  points in `data`.
  for(int npts = 1; npts <= 4; ++npts)
  {
    for(int deg = 0; deg <= npts - 1; ++deg)
    {
      NURBSCurveType curve(data, weights, npts, deg);

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

  double weights[4] = {1.0, 2.0, 3.0, 4.0};

  // clang-format off
  VectorType exp_start_vals[4][4] =  // degree 0                  degree 1                   degree 2                   degree 3
           /* 1 pt */ {{VectorType {0.0, 0.0, 0.0},  VectorType { -999,  -999,  -999},  VectorType { -999,  -999,  -999},  VectorType {-999, -999,  -999}},
           /* 2 pt */  {VectorType {0.0, 0.0, 0.0},  VectorType { -2.8,  -1.6,  -3.2},  VectorType { -999,  -999,  -999},  VectorType {-999, -999,  -999}},
           /* 3 pt */  {VectorType {0.0, 0.0, 0.0},  VectorType {-11.2,  -6.4, -12.8},  VectorType { -3.0,  -2.4, -11.4},  VectorType {-999, -999,  -999}},
           /* 4 pt */  {VectorType {0.0, 0.0, 0.0},  VectorType {-25.2, -14.4, -28.8},  VectorType {-34.0, -20.8, -54.8},  VectorType {-0.6, -2.4, -24.6}}};

  VectorType exp_mid_vals[4][4] =  // degree 0                     degree 1                    degree 2                         degree 3
           /* 1 pt */ {{VectorType {0.0, 0.0, 0.0}, VectorType {      -999,       -999,       -999}, VectorType {   -999,   -999,   -999}, VectorType {   -999,   -999,   -999}},
           /* 2 pt */  {VectorType {0.0, 0.0, 0.0}, VectorType {-112./135., -64.0/135., -128./135.}, VectorType {   -999,   -999,   -999}, VectorType {   -999,   -999,   -999}},
           /* 3 pt */  {VectorType {0.0, 0.0, 0.0}, VectorType {      -9.6,       -4.8,       -3.0}, VectorType { -0.375,   -0.3, -1.425}, VectorType {   -999,   -999,   -999}},
           /* 4 pt */  {VectorType {0.0, 0.0, 0.0}, VectorType {  -11.0592,    -5.5296,     -3.456}, VectorType {-5.1712, 9.5744,  6.144}, VectorType {-3.8448, 0.3456, -0.888}}};

  VectorType exp_end_vals[4][4] =  // degree 0                             degree 1                               degree 2                              degree 3
           /* 1 pt */ {{VectorType {0.0, 0.0, 0.0}, VectorType {     -999,     -999,    -999}, VectorType {  -999,    -999,     -999}, VectorType {   -999,   -999,   -999}},
           /* 2 pt */  {VectorType {0.0, 0.0, 0.0}, VectorType {    -0.35,     -0.2,    -0.4}, VectorType {  -999,    -999,     -999}, VectorType {   -999,   -999,   -999}},
           /* 3 pt */  {VectorType {0.0, 0.0, 0.0}, VectorType {-128./45., -64./45.,  -8./9.}, VectorType {-1./9., -8./90., -19./45.}, VectorType {   -999,   -999,   -999}},
           /* 4 pt */  {VectorType {0.0, 0.0, 0.0}, VectorType {  -1.0125,  -3.7125, -2.3625}, VectorType {  -2.9,    -0.5,     -0.3}, VectorType {-4.0125, 0.4875, 0.3375}}};
  // clang-format on

  // Test evaluation at various degrees and various number of control
  //  points in `data`.
  for(int npts = 1; npts <= 4; ++npts)
  {
    for(int deg = 0; deg <= npts - 1; ++deg)
    {
      NURBSCurveType curve(data, weights, npts, deg);

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

  for(int deg = 1; deg <= 3; ++deg)
  {
    NURBSCurveType curve(data, weights, 4, deg);

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

      if(t <= 0.3)
      {
        PointType p1 = curve1.evaluate(t / 0.3);
        for(int i = 0; i < DIM; ++i)
        {
          EXPECT_NEAR(p[i], p1[i], 1e-13);
        }
      }

      if(t >= 0.3 && t <= 0.6)
      {
        PointType p2 = curve2.evaluate((t - 0.3) / (0.6 - 0.3));
        for(int i = 0; i < DIM; ++i)
        {
          EXPECT_NEAR(p[i], p2[i], 1e-13);
        }
      }

      if(t >= 0.6)
      {
        PointType p3 = curve3.evaluate((t - 0.6) / (1 - 0.6));
        for(int i = 0; i < DIM; ++i)
        {
          EXPECT_NEAR(p[i], p3[i], 1e-13);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_nurbscurve, bezier_extraction)
{
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
  curve.insertKnot(0.33, 3);
  curve.insertKnot(0.66, 1);
  curve.insertKnot(0.77, 2);

  NURBSCurveType bezier_curve;
  auto bezier_list = curve.extractBezier();

  EXPECT_EQ(bezier_list.size(), 4);

  for(double t = 0.0; t < 1.0; t += 0.05)
  {
    PointType p = curve.evaluate(t);

    if(t <= 0.33)
    {
      PointType p1 = bezier_list[0].evaluate(t / 0.33);
      for(int i = 0; i < DIM; ++i)
      {
        EXPECT_NEAR(p[i], p1[i], 1e-13);
      }
    }

    if(t >= 0.33 && t <= 0.66)
    {
      PointType p2 = bezier_list[1].evaluate((t - 0.33) / (0.66 - 0.33));
      for(int i = 0; i < DIM; ++i)
      {
        EXPECT_NEAR(p[i], p2[i], 1e-13);
      }
    }

    if(t >= 0.66 && t <= 0.77)
    {
      PointType p3 = bezier_list[2].evaluate((t - 0.66) / (0.77 - 0.66));
      for(int i = 0; i < DIM; ++i)
      {
        EXPECT_NEAR(p[i], p3[i], 1e-13);
      }
    }

    if(t >= 0.77)
    {
      PointType p4 = bezier_list[3].evaluate((t - 0.77) / (1 - 0.77));
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

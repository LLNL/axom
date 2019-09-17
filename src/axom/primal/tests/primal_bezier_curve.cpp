// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/* 
 * /file primal_bezier_curve.cpp
 * /brief This file tests primal's Bezier curve functionality
 */

#include "gtest/gtest.h"

#include "axom/slic.hpp"

#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/operators/squared_distance.hpp"

namespace primal = axom::primal;

//------------------------------------------------------------------------------
TEST(primal_beziercurve, constructor)
{
  const int DIM = 3;
  using CoordType = double;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;
  using CoordsVec = BezierCurveType::CoordsVec;

  {
    SLIC_INFO("Testing default BezierCurve constructor ");
    BezierCurveType bCurve;

    int expOrder = -1;
    EXPECT_EQ(expOrder, bCurve.getOrder());
    EXPECT_EQ(expOrder + 1, bCurve.getControlPoints().size());
    EXPECT_EQ(CoordsVec(), bCurve.getControlPoints());
  }

  {
    SLIC_INFO("Testing BezierCurve order constructor ");

    BezierCurveType bCurve(1);
    int expOrder = 1;
    EXPECT_EQ(expOrder, bCurve.getOrder());
    EXPECT_EQ(expOrder + 1, static_cast<int>(bCurve.getControlPoints().size()));
  }
}

//----------------------------------------------------------------------------------
TEST(primal_beziercurve, set_order)
{
  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  SLIC_INFO("Test adding control points to empty Bezier curve");

  BezierCurveType bCurve;
  EXPECT_EQ(-1, bCurve.getOrder());

  const int order = 1;
  PointType controlPoints[2] = {PointType::make_point(0.6, 1.2, 1.0),
                                PointType::make_point(0.0, 1.6, 1.8)};

  bCurve.setOrder(order);
  EXPECT_EQ(order, bCurve.getOrder());

  bCurve[0] = controlPoints[0];
  bCurve[1] = controlPoints[1];

  for(int p = 0; p <= bCurve.getOrder(); ++p)
  {
    auto& pt = bCurve[p];
    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_DOUBLE_EQ(controlPoints[p][i], pt[i]);
    }
  }
}

//----------------------------------------------------------------------------------
TEST(primal_beziercurve, point_array_constructor)
{
  SLIC_INFO("Testing point array constructor");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  PointType controlPoints[2] = {PointType::make_point(0.6, 1.2, 1.0),
                                PointType::make_point(0.0, 1.6, 1.8)};

  BezierCurveType bCurve(controlPoints, 1);

  EXPECT_EQ(1, bCurve.getOrder());
  for(int p = 0; p <= bCurve.getOrder(); ++p)
  {
    auto& pt = bCurve[p];
    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_DOUBLE_EQ(controlPoints[p][i], pt[i]);
    }
  }
}

//----------------------------------------------------------------------------------
TEST(primal_beziercurve, coordinate_array_constructor)
{
  SLIC_INFO("Testing coordinate array constructor");

  const int DIM = 3;
  using CoordType = double;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  // clang-format off
  // Note: Order of coordinates is by dimension
  CoordType coords[6] = {0.6, 0.0,  // x-coords for control points
                         1.2, 1.6,  // y-coords for control points
                         1.0, 1.8}; // z-coords for control points
  // clang-format on

  BezierCurveType bCurve(coords, 1);
  EXPECT_EQ(1, bCurve.getOrder());

  EXPECT_DOUBLE_EQ(coords[0], bCurve[0][0]);
  EXPECT_DOUBLE_EQ(coords[2], bCurve[0][1]);
  EXPECT_DOUBLE_EQ(coords[4], bCurve[0][2]);

  EXPECT_DOUBLE_EQ(coords[1], bCurve[1][0]);
  EXPECT_DOUBLE_EQ(coords[3], bCurve[1][1]);
  EXPECT_DOUBLE_EQ(coords[5], bCurve[1][2]);
}

//------------------------------------------------------------------------------
TEST(primal_beziercurve, evaluate)
{
  SLIC_INFO("Testing Bezier evaluation");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  const int order = 3;
  PointType data[order + 1] = {PointType::make_point(0.6, 1.2, 1.0),
                               PointType::make_point(1.3, 1.6, 1.8),
                               PointType::make_point(2.9, 2.4, 2.3),
                               PointType::make_point(3.2, 3.5, 3.0)};

  BezierCurveType b2Curve(data, order);

  PointType midtval = PointType::make_point(2.05, 2.0875, 2.0375);

  // Evaluate the curve at several parameter values
  // Curve should interpolate endpoints
  PointType eval0 = b2Curve.evaluate(0.0);
  PointType eval1 = b2Curve.evaluate(1.0);
  PointType evalMid = b2Curve.evaluate(0.5);

  for(int i = 0; i < DIM; ++i)
  {
    EXPECT_DOUBLE_EQ(b2Curve[0][i], eval0[i]);
    EXPECT_DOUBLE_EQ(b2Curve[order][i], eval1[i]);
    EXPECT_DOUBLE_EQ(midtval[i], evalMid[i]);
  }
}

//------------------------------------------------------------------------------
TEST(primal_beziercurve_, tangent)
{
  SLIC_INFO("Testing Bezier tangent calculation");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  const int order = 3;
  PointType data[order + 1] = {PointType::make_point(0.6, 1.2, 1.0),
                               PointType::make_point(1.3, 1.6, 1.8),
                               PointType::make_point(2.9, 2.4, 2.3),
                               PointType::make_point(3.2, 3.5, 3.0)};

  BezierCurveType b2Curve(data, order);

  PointType midtval = PointType::make_point(3.15, 2.325, 1.875);
  PointType starttval = PointType::make_point(2.1, 1.2, 2.4);
  PointType endtval = PointType::make_point(.9, 3.3, 2.1);

  // Evaluate the curve at several parameter values
  // Curve should be tangent to control net at endpoints
  PointType eval0 = b2Curve.dt(0.0);
  PointType eval1 = b2Curve.dt(1.0);
  PointType evalMid = b2Curve.dt(0.5);

  for(int i = 0; i < DIM; ++i)
  {
    EXPECT_NEAR(starttval[i], eval0[i], 1e-15);
    EXPECT_NEAR(endtval[i], eval1[i], 1e-15);
    EXPECT_NEAR(midtval[i], evalMid[i], 1e-15);
  }
}

//------------------------------------------------------------------------------
TEST(primal_beziercurve, sector_area_cubic)
{
  const int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  {
    SLIC_INFO("Testing Bezier sector area calculation for a cubic");
    const int order = 3;
    PointType data[order + 1] = {PointType::make_point(0.6, 1.2),
                                 PointType::make_point(1.3, 1.6),
                                 PointType::make_point(2.9, 2.4),
                                 PointType::make_point(3.2, 3.5)};

    BezierCurveType bCurve(data, order);
    EXPECT_TRUE(axom::utilities::isNearlyEqual(bCurve.sectorArea(), .1455));
  }
}

//------------------------------------------------------------------------------
TEST(primal_beziercurve, sector_moment_cubic)
{
  const int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  {
    SLIC_INFO("Testing Bezier sector moment calculation for a cubic");
    const int order = 3;
    PointType data[order + 1] = {PointType::make_point(0.6, 1.2),
                                 PointType::make_point(1.3, 1.6),
                                 PointType::make_point(2.9, 2.4),
                                 PointType::make_point(3.2, 3.5)};

    BezierCurveType bCurve(data, order);
    PointType M = bCurve.sectorCentroid();
    EXPECT_NEAR(M[0], -.429321428571429, 2e-15);
    EXPECT_NEAR(M[1], -.354010714285715, 2e-15);
  }
}

//------------------------------------------------------------------------------
TEST(primal_beziercurve, sector_area_point)
{
  const int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  {
    SLIC_INFO("Testing Bezier sector area calculation for a point");
    const int order = 0;
    PointType data[order + 1] = {PointType::make_point(0.6, 1.2)};

    BezierCurveType bCurve(data, order);
    EXPECT_DOUBLE_EQ(bCurve.sectorArea(), 0.0);
  }
}

//------------------------------------------------------------------------------
TEST(primal_beziercurve, sector_moment_point)
{
  const int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  {
    SLIC_INFO("Testing Bezier sector moment calculation for a point");
    const int order = 0;
    PointType data[order + 1] = {PointType::make_point(0.6, 1.2)};

    BezierCurveType bCurve(data, order);
    PointType M = bCurve.sectorCentroid();
    EXPECT_DOUBLE_EQ(M[0], 0.0);
    EXPECT_DOUBLE_EQ(M[1], 0.0);
  }
}

//------------------------------------------------------------------------------
TEST(primal_beziercurve, split_cubic)
{
  SLIC_INFO("Testing Bezier splitting of a cubic");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  const int order = 3;
  PointType data[order + 1] = {PointType::make_point(0.6, 1.2, 1.0),
                               PointType::make_point(1.3, 1.6, 1.8),
                               PointType::make_point(2.9, 2.4, 2.3),
                               PointType::make_point(3.2, 3.5, 3.0)};
  BezierCurveType b2Curve(data, order);

  BezierCurveType b3Curve(order);  // Checks split with order constructor
  BezierCurveType b4Curve;         // Checks split with default constructor
  b2Curve.split(.5, b3Curve, b4Curve);

  // clang-format off
  CoordType b3Coords[12] = {0.6, .95, 1.525, 2.05,
                            1.2, 1.4, 1.7,   2.0875,
                            1.0, 1.4, 1.725, 2.0375};
  CoordType b4Coords[12] = {2.05,   2.575, 3.05, 3.2,
                            2.0875, 2.475, 2.95, 3.5,
                            2.0375, 2.35,  2.65, 3.0};
  // clang-format on

  BezierCurveType b3True(b3Coords, 3);
  BezierCurveType b4True(b4Coords, 3);
  for(int i = 0; i < DIM; ++i)
  {
    for(int p = 0; p <= order; ++p)
    {
      EXPECT_DOUBLE_EQ(b3Curve[p][i], b3True[p][i]);
      EXPECT_DOUBLE_EQ(b4Curve[p][i], b4True[p][i]);
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_beziercurve, split_point)
{
  SLIC_INFO("Testing Bezier splitting for order 0");

  const int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  // Test order-0
  {
    const int order = 0;
    PointType data[order + 1] = {PointType::make_point(0.6, 1.2)};
    BezierCurveType b(data, order);

    BezierCurveType c1, c2;
    b.split(0.5, c1, c2);

    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_DOUBLE_EQ(data[0][i], c1[0][i]);
      EXPECT_DOUBLE_EQ(data[0][i], c2[0][i]);
    }
  }
}

TEST(primal_beziercurve, split_linear)
{
  SLIC_INFO("Testing Bezier splitting for order 1");

  const int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  const int order = 1;
  PointType data[order + 1] = {PointType::make_point(-1, -5),
                               PointType::make_point(1, 5)};
  BezierCurveType b(data, order);

  {
    BezierCurveType c1, c2;
    b.split(0.5, c1, c2);

    EXPECT_DOUBLE_EQ(-1., c1[0][0]);
    EXPECT_DOUBLE_EQ(-5., c1[0][1]);
    EXPECT_DOUBLE_EQ(0., c1[1][0]);
    EXPECT_DOUBLE_EQ(0., c1[1][1]);

    EXPECT_DOUBLE_EQ(0., c2[0][0]);
    EXPECT_DOUBLE_EQ(0., c2[0][1]);
    EXPECT_DOUBLE_EQ(1., c2[1][0]);
    EXPECT_DOUBLE_EQ(5., c2[1][1]);
  }

  {
    BezierCurveType c1, c2;
    const double t = 0.25;
    b.split(0.25, c1, c2);

    PointType interp = PointType::lerp(data[0], data[1], t);

    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_DOUBLE_EQ(data[0][i], c1[0][i]);
      EXPECT_DOUBLE_EQ(interp[i], c1[1][i]);

      EXPECT_DOUBLE_EQ(interp[i], c2[0][i]);
      EXPECT_DOUBLE_EQ(data[1][i], c2[1][i]);
    }
  }
}

TEST(primal_beziercurve, split_quadratic)
{
  SLIC_INFO("Testing Bezier splitting for order 2");

  const int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  const double t = .42;
  const int order = 2;

  // Control points for the three levels of the quadratic de Casteljau algorithm
  PointType lev0[3] = {PointType::make_point(1.1, 1.1),
                       PointType::make_point(5.5, 5.5),
                       PointType::make_point(9.9, 2.2)};

  PointType lev1[2] = {PointType::lerp(lev0[0], lev0[1], t),
                       PointType::lerp(lev0[1], lev0[2], t)};

  PointType lev2[1] = {PointType::lerp(lev1[0], lev1[1], t)};

  BezierCurveType b(lev0, order);

  // Define expected control points for curves 1 and 2
  BezierCurveType expC1(order);
  expC1[0] = lev0[0];
  expC1[1] = lev1[0];
  expC1[2] = lev2[0];

  BezierCurveType expC2(order);
  expC2[0] = lev2[0];
  expC2[1] = lev1[1];
  expC2[2] = lev0[2];

  // Split the curve
  BezierCurveType c1, c2;
  b.split(t, c1, c2);

  SLIC_INFO(""                                          //
            << "Original quadratic: " << b              //
            << "\nCurves after splitting at t = " << t  //
            << "\n\t c1: " << c1                        //
            << "\n\t c2: " << c2);

  // Check values
  for(int p = 0; p <= order; ++p)
  {
    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_DOUBLE_EQ(expC1[p][i], c1[p][i]);
      EXPECT_DOUBLE_EQ(expC2[p][i], c2[p][i]);
    }
  }
}

TEST(primal_beziercurve, isLinear)
{
  SLIC_INFO("Testing isLinear() on Bezier curves");

  const int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using SegmentType = primal::Segment<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  // order 0 -- always true
  {
    const int order = 0;
    auto curve = BezierCurveType(order);
    EXPECT_TRUE(curve.isLinear());

    curve[0] = PointType::make_point(1., 1.);
    EXPECT_TRUE(curve.isLinear());
  }

  // order 1 -- always true
  {
    const int order = 1;
    auto curve = BezierCurveType(order);
    EXPECT_TRUE(curve.isLinear());

    curve[0] = PointType::make_point(1., 1.8);
    curve[1] = PointType::make_point(-12., 3.5);
    EXPECT_TRUE(curve.isLinear());
  }

  // order 2
  {
    const int order = 2;
    auto curve = BezierCurveType(order);
    EXPECT_TRUE(curve.isLinear());

    // straight line
    curve[0] = PointType::make_point(1, 1);
    curve[1] = PointType::make_point(2, 2);
    curve[2] = PointType::make_point(3, 3);
    EXPECT_TRUE(curve.isLinear());

    // move middle point and check linearity with different tolerances
    VectorType v(curve[2], curve[0]);
    auto normal = VectorType::make_vector(-v[1], v[0]);
    curve[1].array() += 0.005 * normal.array();
    SLIC_INFO("Updated curve: " << curve);

    SegmentType s(curve[2], curve[0]);
    SLIC_INFO("Sq dist: " << primal::squared_distance(curve[1], s));

    // linear for a coarse tolerance
    {
      const double tol = 0.1;
      const double tol_sq = tol * tol;
      EXPECT_TRUE(curve.isLinear(tol_sq));
    }

    // non-linear for a finer tolerance
    {
      const double tol = 0.01;
      const double tol_sq = tol * tol;
      EXPECT_FALSE(curve.isLinear(tol_sq));
    }
  }
}

TEST(primal_beziercurve, sector_weights)
{
  SLIC_INFO("Testing weights for BezierCurve::sectorArea()");

  // NOTE: Expected weights are provided in the reference paper [Ueda99]
  // See doxygen comment for BezierCurve::sectorArea()

  using CoordType = double;

  primal::internal::MemoizedSectorAreaWeights<CoordType> memoizedSectorWeights;

  // order 1
  {
    const int ord = 1;
    auto weights = memoizedSectorWeights.getWeights(ord);

    double binomInv = 1. / axom::utilities::binomialCoefficient(2, 1);
    axom::numerics::Matrix<CoordType> exp(ord + 1, ord + 1);
    // clang-format off
    exp(0,0) =  0; exp(0,1) =  1;
    exp(1,0) = -1; exp(1,1) =  0;
    // clang-format on

    for(int i = 0; i <= ord; ++i)
    {
      for(int j = 0; j <= ord; ++j)
      {
        EXPECT_DOUBLE_EQ(exp(i, j) * binomInv, weights(i, j));
      }
    }
  }

  // order 2
  {
    const int ord = 2;
    auto weights = memoizedSectorWeights.getWeights(ord);

    double binomInv = 1. / axom::utilities::binomialCoefficient(4, 2);
    axom::numerics::Matrix<CoordType> exp(ord + 1, ord + 1);
    // clang-format off
    exp(0,0) =  0; exp(0,1) =  2; exp(0,2) =  1;
    exp(1,0) = -2; exp(1,1) =  0; exp(1,2) =  2;
    exp(2,0) = -1; exp(2,1) = -2; exp(2,2) =  0;
    // clang-format on

    for(int i = 0; i <= ord; ++i)
    {
      for(int j = 0; j <= ord; ++j)
      {
        EXPECT_DOUBLE_EQ(exp(i, j) * binomInv, weights(i, j));
      }
    }
  }

  // order 3
  {
    const int ord = 3;
    auto weights = memoizedSectorWeights.getWeights(ord);

    double binomInv = 1. / axom::utilities::binomialCoefficient(6, 3);
    axom::numerics::Matrix<CoordType> exp(ord + 1, ord + 1);
    // clang-format off
    exp(0,0) =  0; exp(0,1) =  6; exp(0,2) =  3; exp(0,3) =  1;
    exp(1,0) = -6; exp(1,1) =  0; exp(1,2) =  3; exp(1,3) =  3;
    exp(2,0) = -3; exp(2,1) = -3; exp(2,2) =  0; exp(2,3) =  6;
    exp(3,0) = -1; exp(3,1) = -3; exp(3,2) = -6; exp(3,3) =  0;
    // clang-format on

    for(int i = 0; i <= ord; ++i)
    {
      for(int j = 0; j <= ord; ++j)
      {
        EXPECT_DOUBLE_EQ(exp(i, j) * binomInv, weights(i, j));
      }
    }
  }

  // order 4
  {
    const int ord = 4;
    auto weights = memoizedSectorWeights.getWeights(ord);

    double binomInv = 1. / axom::utilities::binomialCoefficient(8, 4);
    axom::numerics::Matrix<CoordType> exp(ord + 1, ord + 1);
    // clang-format off
    exp(0,0) =  0; exp(0,1) = 20; exp(0,2) = 10; exp(0,3) =  4; exp(0,4) =  1;
    exp(1,0) =-20; exp(1,1) =  0; exp(1,2) =  8; exp(1,3) =  8; exp(1,4) =  4;
    exp(2,0) =-10; exp(2,1) = -8; exp(2,2) =  0; exp(2,3) =  8; exp(2,4) = 10;
    exp(3,0) = -4; exp(3,1) = -8; exp(3,2) = -8; exp(3,3) =  0; exp(3,4) = 20;
    exp(4,0) = -1; exp(4,1) = -4; exp(4,2) =-10; exp(4,3) =-20; exp(4,4) =  0;
    // clang-format on

    for(int i = 0; i <= ord; ++i)
    {
      for(int j = 0; j <= ord; ++j)
      {
        EXPECT_DOUBLE_EQ(exp(i, j) * binomInv, weights(i, j));
      }
    }
  }

  // order 5
  {
    const int ord = 5;
    auto weights = memoizedSectorWeights.getWeights(ord);

    double binomInv = 1. / axom::utilities::binomialCoefficient(10, 5);
    axom::numerics::Matrix<CoordType> exp(ord + 1, ord + 1);
    // clang-format off
    exp(0,0) =  0; exp(0,1) = 70; exp(0,2) = 35; exp(0,3) = 15; exp(0,4) =  5; exp(0,5) =  1;
    exp(1,0) =-70; exp(1,1) =  0; exp(1,2) = 25; exp(1,3) = 25; exp(1,4) = 15; exp(1,5) =  5;
    exp(2,0) =-35; exp(2,1) =-25; exp(2,2) =  0; exp(2,3) = 20; exp(2,4) = 25; exp(2,5) = 15;
    exp(3,0) =-15; exp(3,1) =-25; exp(3,2) =-20; exp(3,3) =  0; exp(3,4) = 25; exp(3,5) = 35;
    exp(4,0) = -5; exp(4,1) =-15; exp(4,2) =-25; exp(4,3) =-25; exp(4,4) =  0; exp(4,5) = 70;
    exp(5,0) = -1; exp(5,1) = -5; exp(5,2) =-15; exp(5,3) =-35; exp(5,4) =-70; exp(5,5) =  0;
    // clang-format on

    for(int i = 0; i <= ord; ++i)
    {
      for(int j = 0; j <= ord; ++j)
      {
        EXPECT_DOUBLE_EQ(exp(i, j) * binomInv, weights(i, j));
      }
    }
  }
}

//------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}

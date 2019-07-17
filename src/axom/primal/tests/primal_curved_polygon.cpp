// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/* /file bezier_test.cpp
 * /brief This file tests the BezierCurve.hpp and eval_bezier.hpp files
 */

#include "gtest/gtest.h"

#include "axom/slic.hpp"
#include "axom/primal/geometry/CurvedPolygon.hpp"

namespace primal = axom::primal;

std::map<int, std::vector<double>> primal::UedaAreaMats;  //TODO: put this in the BezierCurve.cpp file
//------------------------------------------------------------------------------
TEST(primal_curvedpolygon, constructor)
{
  const int DIM = 3;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  {
    SLIC_INFO("Testing default CurvedPolygon constructor ");
    CurvedPolygonType bPolygon;

    int expNumEdges = 0;
    EXPECT_EQ(expNumEdges, bPolygon.numEdges());
    EXPECT_EQ(expNumEdges, bPolygon.getEdges().size());
    EXPECT_EQ(std::vector<BezierCurveType>(), bPolygon.getEdges());
  }

  {
    SLIC_INFO("Testing CurvedPolygon order constructor ");

    CurvedPolygonType bPolygon(1);
    int expNumEdges = 1;
    EXPECT_EQ(expNumEdges, bPolygon.numEdges());
    EXPECT_EQ(expNumEdges, static_cast<int>(bPolygon.getEdges().size()));
  }
}

//----------------------------------------------------------------------------------
TEST(primal_curvedpolygon, add_edges)
{
  const int DIM = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  SLIC_INFO("Test adding edges to empty CurvedPolygon");

  CurvedPolygonType bPolygon;
  EXPECT_EQ(0, bPolygon.numEdges());

  PointType controlPoints[2] = {PointType::make_point(0.6, 1.2),
                                PointType::make_point(0.0, 1.6)};

  BezierCurveType bCurve(controlPoints, 1);

  bPolygon.addEdge(bCurve);
  bPolygon.addEdge(bCurve);

  EXPECT_EQ(2, bPolygon.numEdges());
  for(int p = 0; p < bPolygon.numEdges(); ++p)
  {
    BezierCurveType& bc = bPolygon[p];
    for(int sz = 0; sz <= bc.getOrder(); ++sz)
    {
      auto& pt = bc[sz];
      for(int i = 0; i < DIM; ++i)
      {
        EXPECT_DOUBLE_EQ(controlPoints[sz][i], pt[i]);
      }
    }
  }
}

//----------------------------------------------------------------------------------
TEST(primal_curvedpolygon, is_Valid)
{
  const int DIM = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  SLIC_INFO("Test checking if CurvedPolygon is closed.");

  CurvedPolygonType bPolygon;
  EXPECT_EQ(0, bPolygon.numEdges());

  PointType controlPoints[2] = {PointType::make_point(0.6, 1.2),
                                PointType::make_point(0.0, 1.6)};

  PointType controlPoints2[2] = {PointType::make_point(0.0, 1.6),
                                 PointType::make_point(0.3, 2.0)};

  PointType controlPoints3[2] = {PointType::make_point(0.3, 2.0),
                                 PointType::make_point(0.6, 1.2)};

  BezierCurveType bCurve(controlPoints, 1);
  bPolygon.addEdge(bCurve);

  BezierCurveType bCurve2(controlPoints2, 1);
  bPolygon.addEdge(bCurve2);

  EXPECT_EQ(2, bPolygon.numEdges());
  EXPECT_EQ(false, bPolygon.isClosed());

  BezierCurveType bCurve3(controlPoints3, 1);
  bPolygon.addEdge(bCurve3);

  EXPECT_EQ(3, bPolygon.numEdges());
  EXPECT_EQ(true, bPolygon.isClosed());
}

//----------------------------------------------------------------------------------
TEST(primal_beziercurve, area)
{
  const int DIM = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  SLIC_INFO("Test checking if CurvedPolygon is closed.");

  CurvedPolygonType bPolygon;
  EXPECT_EQ(0, bPolygon.numEdges());

  PointType controlPoints[2] = {PointType::make_point(0.6, 1.2),
                                PointType::make_point(0.0, 1.6)};

  PointType controlPoints2[2] = {PointType::make_point(0.0, 1.6),
                                 PointType::make_point(0.3, 2.0)};

  PointType controlPoints3[2] = {PointType::make_point(0.3, 2.0),
                                 PointType::make_point(0.6, 1.2)};

  BezierCurveType bCurve(controlPoints, 1);
  bPolygon.addEdge(bCurve);

  BezierCurveType bCurve2(controlPoints2, 1);
  bPolygon.addEdge(bCurve2);

  BezierCurveType bCurve3(controlPoints3, 1);
  bPolygon.addEdge(bCurve3);

  CoordType A = bPolygon.area();
  CoordType trueA = .18;

  EXPECT_TRUE(axom::utilities::isNearlyEqual(trueA, A));
}

/*
//----------------------------------------------------------------------------------
TEST( primal_beziercurve, coordinate_array_constructor )
{
  SLIC_INFO("Testing coordinate array constructor");

  const int DIM = 3;
  using CoordType = double;
  using BezierCurveType = primal::BezierCurve< CoordType, DIM >;

  // Note: Order of coordinates is by dimension
  CoordType coords[6] = {0.6, 0.0,  // x-coords for control points
                         1.2, 1.6,  // y-coords for control points
                         1.0, 1.8}; // z-coords for control points

  BezierCurveType bCurve(coords,1);
  EXPECT_EQ(1, bCurve.getOrder());


  EXPECT_DOUBLE_EQ(coords[0], bCurve[0][0]);
  EXPECT_DOUBLE_EQ(coords[2], bCurve[0][1]);
  EXPECT_DOUBLE_EQ(coords[4], bCurve[0][2]);

  EXPECT_DOUBLE_EQ(coords[1], bCurve[1][0]);
  EXPECT_DOUBLE_EQ(coords[3], bCurve[1][1]);
  EXPECT_DOUBLE_EQ(coords[5], bCurve[1][2]);
}

//------------------------------------------------------------------------------
TEST( primal_beziercurve, evaluate)
{
  SLIC_INFO("Testing Bezier evaluation");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point< CoordType, DIM >;
  using BezierCurveType = primal::BezierCurve< CoordType, DIM >;

  const int order = 3;
  PointType data[order+1] =  { PointType::make_point(0.6, 1.2, 1.0),
                               PointType::make_point(1.3, 1.6, 1.8),
                               PointType::make_point(2.9, 2.4, 2.3),
                               PointType::make_point(3.2, 3.5, 3.0) };

  BezierCurveType b2Curve(data, order);

  PointType midtval = PointType::make_point(2.05,2.0875,2.0375);

  // Evaluate the curve at several parameter values
  // Curve should interpolate endpoints
  PointType eval0 = b2Curve.evaluate(0.0);
  PointType eval1 = b2Curve.evaluate(1.0);
  PointType evalMid = b2Curve.evaluate(0.5);

  for ( int i=0 ; i<DIM ; ++i)
  {
    EXPECT_DOUBLE_EQ(b2Curve[0][i],     eval0[i]);
    EXPECT_DOUBLE_EQ(b2Curve[order][i], eval1[i]);
    EXPECT_DOUBLE_EQ(midtval[i],        evalMid[i]);
  }
}

//------------------------------------------------------------------------------
TEST( primal_beziercurve, split_cubic )
{
  SLIC_INFO("Testing Bezier splitting of a cubic");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point< CoordType, DIM >;
  using BezierCurveType = primal::BezierCurve< CoordType, DIM >;

  const int order = 3;
  PointType data[order+1] = { PointType::make_point(0.6, 1.2, 1.0),
                              PointType::make_point(1.3, 1.6, 1.8),
                              PointType::make_point(2.9, 2.4, 2.3),
                              PointType::make_point(3.2, 3.5, 3.0) };
  BezierCurveType b2Curve(data, order);

  BezierCurveType b3Curve(order); // Checks split with order constructor
  BezierCurveType b4Curve;        // Checks split with default constructor
  b2Curve.split(.5, b3Curve, b4Curve);

  CoordType b3Coords[12] = {0.6, .95, 1.525, 2.05,
                            1.2, 1.4, 1.7,   2.0875,
                            1.0, 1.4, 1.725, 2.0375};
  CoordType b4Coords[12] = {2.05,   2.575, 3.05, 3.2,
                            2.0875, 2.475, 2.95, 3.5,
                            2.0375, 2.35,  2.65, 3.0};
  BezierCurveType b3True(b3Coords,3);
  BezierCurveType b4True(b4Coords,3);
  for ( int i=0 ; i<DIM ; ++i)
  {
    for ( int p=0 ; p<= order ; ++p)
    {
      EXPECT_DOUBLE_EQ(b3Curve[p][i], b3True[p][i]);
      EXPECT_DOUBLE_EQ(b4Curve[p][i], b4True[p][i]);
    }
  }
}

//------------------------------------------------------------------------------
TEST( primal_beziercurve, split_point )
{
  SLIC_INFO("Testing Bezier splitting for order 0");

  const int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point< CoordType, DIM >;
  using BezierCurveType = primal::BezierCurve< CoordType, DIM >;

  // Test order-0
  {
    const int order = 0;
    PointType data[order+1] = { PointType::make_point(0.6, 1.2)};
    BezierCurveType b(data, order);

    BezierCurveType c1,c2;
    b.split(0.5, c1, c2);

    for ( int i=0 ; i<DIM ; ++i)
    {
      EXPECT_DOUBLE_EQ(data[0][i], c1[0][i]);
      EXPECT_DOUBLE_EQ(data[0][i], c2[0][i]);
    }
  }
}

TEST( primal_beziercurve, split_linear )
{
  SLIC_INFO("Testing Bezier splitting for order 1");

  const int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point< CoordType, DIM >;
  using BezierCurveType = primal::BezierCurve< CoordType, DIM >;

  const int order = 1;
  PointType data[order+1] = { PointType::make_point(-1, -5),
                              PointType::make_point( 1,  5) };
  BezierCurveType b(data, order);

  {
    BezierCurveType c1,c2;
    b.split(0.5, c1, c2);

    EXPECT_DOUBLE_EQ(-1., c1[0][0]);
    EXPECT_DOUBLE_EQ(-5., c1[0][1]);
    EXPECT_DOUBLE_EQ( 0., c1[1][0]);
    EXPECT_DOUBLE_EQ( 0., c1[1][1]);

    EXPECT_DOUBLE_EQ( 0., c2[0][0]);
    EXPECT_DOUBLE_EQ( 0., c2[0][1]);
    EXPECT_DOUBLE_EQ( 1., c2[1][0]);
    EXPECT_DOUBLE_EQ( 5., c2[1][1]);

  }

  {
    BezierCurveType c1,c2;
    const double t = 0.25;
    b.split(0.25, c1, c2);

    PointType interp = PointType::lerp(data[0], data[1], t);

    for ( int i=0 ; i<DIM ; ++i)
    {
      EXPECT_DOUBLE_EQ(data[0][i], c1[0][i]);
      EXPECT_DOUBLE_EQ(interp[i],  c1[1][i]);

      EXPECT_DOUBLE_EQ(interp[i],  c2[0][i]);
      EXPECT_DOUBLE_EQ(data[1][i], c2[1][i]);
    }
  }
}

TEST( primal_beziercurve, split_quadratic)
{
  SLIC_INFO("Testing Bezier splitting for order 2");

  const int DIM = 2;
  using CoordType = double;
  using PointType = primal::Point< CoordType, DIM >;
  using BezierCurveType = primal::BezierCurve< CoordType, DIM >;

  const double t = .42;
  const int order = 2;

  // Control points for the three levels of the quadratic de Casteljau algorithm
  PointType lev0[3] = { PointType::make_point(1.1, 1.1),
                        PointType::make_point( 5.5, 5.5),
                        PointType::make_point( 9.9, 2.2) };

  PointType lev1[2] = { PointType::lerp(lev0[0], lev0[1], t),
                        PointType::lerp(lev0[1], lev0[2], t) };

  PointType lev2[1] = { PointType::lerp(lev1[0], lev1[1], t) };

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
  BezierCurveType c1,c2;
  b.split(t, c1, c2);

  SLIC_INFO(""
            <<"Original quadratic: "<< b
            << "\nCurves after splitting at t = "<< t
            << "\n\t c1: " << c1
            << "\n\t c2: " << c2);

  // Check values
  for(int p=0 ; p <= order ; ++p)
  {
    for ( int i=0 ; i<DIM ; ++i)
    {
      EXPECT_DOUBLE_EQ(expC1[p][i], c1[p][i]);
      EXPECT_DOUBLE_EQ(expC2[p][i], c2[p][i]);
    }
  }
}


//------------------------------------------------------------------------------
*/
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}

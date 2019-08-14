// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/* /file bezier_test.cpp
 * /brief This file tests the BezierCurve.hpp and eval_bezier.hpp files
 */

#include "gtest/gtest.h"

#include "axom/slic.hpp"
#include "axom/primal/geometry/CurvedPolygon.hpp"
#include "axom/primal/operators/intersect_curved_poly.hpp"

namespace primal = axom::primal;

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
    SLIC_INFO("Testing CurvedPolygon numEdges constructor ");

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
TEST(primal_curvedpolygon, isClosed)
{
  const int DIM = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  SLIC_INFO("Test checking if CurvedPolygon is closed.");

  CurvedPolygonType bPolygon;
  EXPECT_EQ(0, bPolygon.numEdges());
  EXPECT_EQ(false, bPolygon.isClosed());

  PointType controlPoints[2] = {PointType::make_point(0.6, 1.2),
                                PointType::make_point(0.3, 2.0)};

  PointType controlPoints2[2] = {PointType::make_point(0.3, 2.0),
                                 PointType::make_point(0.0, 1.6)};

  PointType controlPoints3[2] = {PointType::make_point(0.0, 1.6),
                                 PointType::make_point(0.6, 1.2)};

  BezierCurveType bCurve(controlPoints, 1);
  bPolygon.addEdge(bCurve);
  EXPECT_EQ(false, bPolygon.isClosed());

  BezierCurveType bCurve2(controlPoints2, 1);
  bPolygon.addEdge(bCurve2);

  EXPECT_EQ(2, bPolygon.numEdges());
  EXPECT_EQ(false, bPolygon.isClosed());

  BezierCurveType bCurve3(controlPoints3, 1);
  bPolygon.addEdge(bCurve3);

  EXPECT_EQ(3, bPolygon.numEdges());
  EXPECT_EQ(true, bPolygon.isClosed());

  bPolygon[2][1][0] -= 2e-15;
  EXPECT_EQ(false, bPolygon.isClosed(1e-15));
}

//----------------------------------------------------------------------------------
TEST(primal_curvedpolygon, split_edge)
{
  const int DIM = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  SLIC_INFO("Test checking CurvedPolygon edge split.");

  CurvedPolygonType bPolygon;
  EXPECT_EQ(0, bPolygon.numEdges());

  PointType controlPoints[2] = {PointType::make_point(0.6, 1.2),
                                PointType::make_point(0.3, 2.0)};

  PointType controlPoints2[2] = {PointType::make_point(0.3, 2.0),
                                 PointType::make_point(0.0, 1.6)};

  PointType controlPoints3[2] = {PointType::make_point(0.0, 1.6),
                                 PointType::make_point(0.6, 1.2)};

  BezierCurveType bCurve(controlPoints, 1);
  bPolygon.addEdge(bCurve);

  BezierCurveType bCurve2(controlPoints2, 1);
  bPolygon.addEdge(bCurve2);

  BezierCurveType bCurve3(controlPoints3, 1);
  bPolygon.addEdge(bCurve3);

  bPolygon.splitEdge(0, .5);
  bCurve.split(.5, bCurve2, bCurve3);

  EXPECT_EQ(bPolygon.numEdges(), 4);
  for(int i = 0; i < bPolygon[0].getOrder(); ++i)
  {
    for(int dimi = 0; dimi < DIM; ++dimi)
    {
      EXPECT_EQ(bPolygon[0][i][dimi], bCurve2[i][dimi]);
      EXPECT_EQ(bPolygon[1][i][dimi], bCurve3[i][dimi]);
    }
  }
}

//----------------------------------------------------------------------------------
TEST(primal_curvedpolygon, area_triangle_degenerate)
{
  const int DIM = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  SLIC_INFO(
    "Test checking CurvedPolygon degenerate triangle area computation.");

  CurvedPolygonType bPolygon;
  EXPECT_EQ(0, bPolygon.numEdges());
  EXPECT_EQ(0.0, bPolygon.area());

  PointType controlPoints[2] = {PointType::make_point(0.6, 1.2),
                                PointType::make_point(0.3, 2.0)};

  PointType controlPoints2[2] = {PointType::make_point(0.3, 2.0),
                                 PointType::make_point(0.0, 1.6)};

  PointType controlPoints3[2] = {PointType::make_point(0.0, 1.6),
                                 PointType::make_point(0.6, 1.2)};

  BezierCurveType bCurve(controlPoints, 1);
  bPolygon.addEdge(bCurve);
  EXPECT_EQ(0.0, bPolygon.area());

  BezierCurveType bCurve2(controlPoints2, 1);
  bPolygon.addEdge(bCurve2);
  EXPECT_EQ(0.0, bPolygon.area());

  BezierCurveType bCurve3(controlPoints3, 1);
  bPolygon.addEdge(bCurve3);

  bPolygon[2][1][0] -= 2e-15;
  EXPECT_EQ(0.0, bPolygon.area());
}

//----------------------------------------------------------------------------------
TEST(primal_curvedpolygon, area_triangle_linear)
{
  const int DIM = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  SLIC_INFO("Test checking CurvedPolygon linear triangle area computation.");

  CurvedPolygonType bPolygon;
  EXPECT_EQ(0, bPolygon.numEdges());

  PointType controlPoints[2] = {PointType::make_point(0.6, 1.2),
                                PointType::make_point(0.3, 2.0)};

  PointType controlPoints2[2] = {PointType::make_point(0.3, 2.0),
                                 PointType::make_point(0.0, 1.6)};

  PointType controlPoints3[2] = {PointType::make_point(0.0, 1.6),
                                 PointType::make_point(0.6, 1.2)};

  BezierCurveType bCurve(controlPoints, 1);
  bPolygon.addEdge(bCurve);

  BezierCurveType bCurve2(controlPoints2, 1);
  bPolygon.addEdge(bCurve2);

  BezierCurveType bCurve3(controlPoints3, 1);
  bPolygon.addEdge(bCurve3);

  CoordType A = bPolygon.area();
  CoordType trueA = -.18;

  EXPECT_DOUBLE_EQ(trueA, A);
}

//----------------------------------------------------------------------------------
TEST(primal_curvedpolygon, area_triangle_quadratic)
{
  const int DIM = 2;
  const int order = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  SLIC_INFO("Test checking CurvedPolygon quadratic triangle area computation.");

  CurvedPolygonType bPolygon;
  EXPECT_EQ(0, bPolygon.numEdges());

  PointType controlPoints[order + 1] = {PointType::make_point(0.6, 1.2),
                                        PointType::make_point(0.4, 1.3),
                                        PointType::make_point(0.3, 2.0)};

  PointType controlPoints2[order + 1] = {PointType::make_point(0.3, 2.0),
                                         PointType::make_point(0.27, 1.5),
                                         PointType::make_point(0.0, 1.6)};

  PointType controlPoints3[order + 1] = {PointType::make_point(0.0, 1.6),
                                         PointType::make_point(0.1, 1.5),
                                         PointType::make_point(0.6, 1.2)};

  BezierCurveType bCurve(controlPoints, order);
  bPolygon.addEdge(bCurve);

  BezierCurveType bCurve2(controlPoints2, order);
  bPolygon.addEdge(bCurve2);

  BezierCurveType bCurve3(controlPoints3, order);
  bPolygon.addEdge(bCurve3);

  CoordType A = bPolygon.area();

  CoordType trueA = -.09733333333333333333;
  EXPECT_DOUBLE_EQ(trueA, A);
}

//----------------------------------------------------------------------------------
TEST(primal_curvedpolygon, area_triangle_mixed_order)
{
  const int DIM = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  SLIC_INFO(
    "Test checking CurvedPolygon mixed order triangle area computation.");

  CurvedPolygonType bPolygon;
  EXPECT_EQ(0, bPolygon.numEdges());

  PointType controlPoints[3] = {PointType::make_point(0.6, 1.2),
                                PointType::make_point(0.4, 1.3),
                                PointType::make_point(0.3, 2.0)};

  PointType controlPoints2[3] = {PointType::make_point(0.3, 2.0),
                                 PointType::make_point(0.27, 1.5),
                                 PointType::make_point(0.0, 1.6)};

  PointType controlPoints3[2] = {PointType::make_point(0.0, 1.6),
                                 PointType::make_point(0.6, 1.2)};

  BezierCurveType bCurve(controlPoints, 2);
  bPolygon.addEdge(bCurve);

  BezierCurveType bCurve2(controlPoints2, 2);
  bPolygon.addEdge(bCurve2);

  BezierCurveType bCurve3(controlPoints3, 1);
  bPolygon.addEdge(bCurve3);

  CoordType A = bPolygon.area();

  CoordType trueA = -.0906666666666666666666;
  EXPECT_DOUBLE_EQ(trueA, A);
}

//----------------------------------------------------------------------------------
TEST(primal_curvedpolygon, moment_triangle_linear)
{
  const int DIM = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  SLIC_INFO("Test checking CurvedPolygon linear triangle moment computation.");

  CurvedPolygonType bPolygon;
  EXPECT_EQ(0, bPolygon.numEdges());

  PointType controlPoints[2] = {PointType::make_point(0.6, 1.2),
                                PointType::make_point(0.3, 2.0)};

  PointType controlPoints2[2] = {PointType::make_point(0.3, 2.0),
                                 PointType::make_point(0.0, 1.6)};
  PointType controlPoints3[2] = {PointType::make_point(0.0, 1.6),
                                 PointType::make_point(0.6, 1.2)};

  BezierCurveType bCurve(controlPoints, 1);
  bPolygon.addEdge(bCurve);

  BezierCurveType bCurve2(controlPoints2, 1);
  bPolygon.addEdge(bCurve2);

  BezierCurveType bCurve3(controlPoints3, 1);
  bPolygon.addEdge(bCurve3);

  PointType M = bPolygon.moment();
  CoordType trueM1 = 0.3;
  CoordType trueM2 = 1.6;

  EXPECT_DOUBLE_EQ(trueM1, M[0]);
  EXPECT_DOUBLE_EQ(trueM2, M[1]);
}

//----------------------------------------------------------------------------------
TEST(primal_curvedpolygon, moment_triangle_mixed_order)
{
  const int DIM = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  SLIC_INFO(
    "Test checking CurvedPolygon mixed order triangle area computation.");

  CurvedPolygonType bPolygon;
  EXPECT_EQ(0, bPolygon.numEdges());

  PointType controlPoints[3] = {PointType::make_point(0.6, 1.2),
                                PointType::make_point(0.4, 1.3),
                                PointType::make_point(0.3, 2.0)};

  PointType controlPoints2[3] = {PointType::make_point(0.3, 2.0),
                                 PointType::make_point(0.27, 1.5),
                                 PointType::make_point(0.0, 1.6)};

  PointType controlPoints3[2] = {PointType::make_point(0.0, 1.6),
                                 PointType::make_point(0.6, 1.2)};

  BezierCurveType bCurve(controlPoints, 2);
  bPolygon.addEdge(bCurve);

  BezierCurveType bCurve2(controlPoints2, 2);
  bPolygon.addEdge(bCurve2);

  BezierCurveType bCurve3(controlPoints3, 1);
  bPolygon.addEdge(bCurve3);

  PointType M = bPolygon.moment();
  CoordType trueM2 = 1.55764705882353;
  CoordType trueM1 = .2970147058823527;

  EXPECT_DOUBLE_EQ(trueM1, M[0]);
  EXPECT_DOUBLE_EQ(trueM2, M[1]);
}

//----------------------------------------------------------------------------------
TEST(primal_curvedpolygon, intersection_triangle_linear)
{
  const int DIM = 2;
  const int order = 1;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  SLIC_INFO(
    "Test intersecting two linear triangular CurvedPolygons (single region).");

  CurvedPolygonType bPolygon;
  EXPECT_EQ(0, bPolygon.numEdges());

  PointType controlPoints[order + 1] = {PointType::make_point(0.6, 1.2),
                                        PointType::make_point(0.3, 2.0)};

  PointType controlPoints2[order + 1] = {PointType::make_point(0.3, 2.0),
                                         PointType::make_point(0.0, 1.6)};

  PointType controlPoints3[order + 1] = {PointType::make_point(0.0, 1.6),
                                         PointType::make_point(0.6, 1.2)};

  BezierCurveType bCurve(controlPoints, order);
  bPolygon.addEdge(bCurve);

  BezierCurveType bCurve2(controlPoints2, order);
  bPolygon.addEdge(bCurve2);

  BezierCurveType bCurve3(controlPoints3, order);
  bPolygon.addEdge(bCurve3);
  CurvedPolygonType bPolygon2 = bPolygon;
  for(int i = 0; i < DIM; ++i)
  {
    for(int j = 0; j <= order; ++j)
    {
      for(int k = 0; k < bPolygon2.numEdges(); ++k)
      {
        bPolygon2[k][j][i] += .11;
      }
    }
  }

  PointType expcontrolPoints[order + 1] = {
    PointType::make_point(0.3091666666666666666666, 1.9755555555555555555),
    PointType::make_point(0.11, 1.71)};

  PointType expcontrolPoints2[order + 1] = {
    PointType::make_point(0.11, 1.71),
    PointType::make_point(0.5083333333333333333, 1.44444444444444444444)};

  PointType expcontrolPoints3[order + 1] = {
    PointType::make_point(0.5083333333333333333, 1.44444444444444444444),
    PointType::make_point(0.3091666666666666666666, 1.9755555555555555555)};

  CurvedPolygonType expbPolygon;

  BezierCurveType expbCurve(expcontrolPoints, order);
  expbPolygon.addEdge(expbCurve);
  BezierCurveType expbCurve2(expcontrolPoints2, order);
  expbPolygon.addEdge(expbCurve2);
  BezierCurveType expbCurve3(expcontrolPoints3, order);
  expbPolygon.addEdge(expbCurve3);

  std::vector<CurvedPolygonType> bPolygons3;
  bool didIntersect = intersect_polygon(bPolygon, bPolygon2, bPolygons3);
  EXPECT_TRUE(didIntersect);

  for(int i = 0; i < DIM; ++i)
  {
    for(int j = 0; j <= order; ++j)
    {
      for(int idxcurve = 0; idxcurve < static_cast<int>(bPolygons3.size());
          ++idxcurve)
      {
        for(int k = 0; k < bPolygons3[idxcurve].numEdges(); ++k)
        {
          EXPECT_TRUE(axom::utilities::isNearlyEqual(expbPolygon[k][j][i],
                                                     bPolygons3[idxcurve][k][j][i],
                                                     1e-14));
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------
TEST(primal_curvedpolygon, intersection_triangle_quadratic)
{
  const int DIM = 2;
  const int order = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  SLIC_INFO(
    "Test intersecting two quadratic triangular CurvedPolygons (single "
    "region).");

  CurvedPolygonType bPolygon;
  EXPECT_EQ(0, bPolygon.numEdges());

  PointType controlPoints[order + 1] = {PointType::make_point(0.6, 1.2),
                                        PointType::make_point(0.4, 1.3),
                                        PointType::make_point(0.3, 2.0)};

  PointType controlPoints2[order + 1] = {PointType::make_point(0.3, 2.0),
                                         PointType::make_point(0.27, 1.5),
                                         PointType::make_point(0.0, 1.6)};

  PointType controlPoints3[order + 1] = {PointType::make_point(0.0, 1.6),
                                         PointType::make_point(0.1, 1.5),
                                         PointType::make_point(0.6, 1.2)};

  BezierCurveType bCurve(controlPoints, order);
  bPolygon.addEdge(bCurve);

  BezierCurveType bCurve2(controlPoints2, order);
  bPolygon.addEdge(bCurve2);

  BezierCurveType bCurve3(controlPoints3, order);
  bPolygon.addEdge(bCurve3);
  CurvedPolygonType bPolygon2 = bPolygon;
  for(int i = 0; i < DIM; ++i)
  {
    for(int j = 0; j <= order; ++j)
    {
      for(int k = 0; k < bPolygon2.numEdges(); ++k)
      {
        bPolygon2[k][j][i] += .11;
      }
    }
  }

  PointType expcontrolPoints[order + 1] = {
    PointType::make_point(0.335956890729522, 1.784126953773395),
    PointType::make_point(0.297344765794753, 1.718171485335525),
    PointType::make_point(0.239567753301698, 1.700128235793372)};

  PointType expcontrolPoints2[order + 1] = {
    PointType::make_point(0.2395677533016981, 1.700128235793371),
    PointType::make_point(0.221884203146682, 1.662410644580941),
    PointType::make_point(0.199328465398189, 1.636873522352205)};

  PointType expcontrolPoints3[order + 1] = {
    PointType::make_point(0.199328465398188, 1.636873522352206),
    PointType::make_point(0.277429214338182, 1.579562422716502),
    PointType::make_point(0.408882616650578, 1.495574996394597)};

  PointType expcontrolPoints4[order + 1] = {
    PointType::make_point(0.408882616650588, 1.495574996394586),
    PointType::make_point(0.368520120719339, 1.616453177259694),
    PointType::make_point(0.335956890729522, 1.784126953773394)};

  CurvedPolygonType expbPolygon;

  BezierCurveType expbCurve(expcontrolPoints, order);
  expbPolygon.addEdge(expbCurve);
  BezierCurveType expbCurve2(expcontrolPoints2, order);
  expbPolygon.addEdge(expbCurve2);
  BezierCurveType expbCurve3(expcontrolPoints3, order);
  expbPolygon.addEdge(expbCurve3);
  BezierCurveType expbCurve4(expcontrolPoints4, order);
  expbPolygon.addEdge(expbCurve4);

  std::vector<CurvedPolygonType> bPolygons3;
  bool didIntersect = intersect_polygon(bPolygon, bPolygon2, bPolygons3);
  EXPECT_TRUE(didIntersect);

  for(int i = 0; i < DIM; ++i)
  {
    for(int j = 0; j <= order; ++j)
    {
      for(int idxcurve = 0; idxcurve < static_cast<int>(bPolygons3.size());
          ++idxcurve)
      {
        for(int k = 0; k < bPolygons3[idxcurve].numEdges(); ++k)
        {
          EXPECT_TRUE(axom::utilities::isNearlyEqual(expbPolygon[k][j][i],
                                                     bPolygons3[idxcurve][k][j][i],
                                                     1e-14));
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------
TEST(primal_curvedpolygon, area_intersection_triangle_linear)
{
  const int DIM = 2;
  const int order = 1;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  SLIC_INFO(
    "Test finding area of intersection two linear triangular CurvedPolygons "
    "(single region).");

  CurvedPolygonType bPolygon;
  EXPECT_EQ(0, bPolygon.numEdges());

  PointType controlPoints[order + 1] = {PointType::make_point(0.6, 1.2),
                                        PointType::make_point(0.3, 2.0)};

  PointType controlPoints2[order + 1] = {PointType::make_point(0.3, 2.0),
                                         PointType::make_point(0.0, 1.6)};

  PointType controlPoints3[order + 1] = {PointType::make_point(0.0, 1.6),
                                         PointType::make_point(0.6, 1.2)};

  BezierCurveType bCurve(controlPoints, order);
  bPolygon.addEdge(bCurve);

  BezierCurveType bCurve2(controlPoints2, order);
  bPolygon.addEdge(bCurve2);

  BezierCurveType bCurve3(controlPoints3, order);
  bPolygon.addEdge(bCurve3);
  CurvedPolygonType bPolygon2 = bPolygon;
  for(int i = 0; i < DIM; ++i)
  {
    for(int j = 0; j <= order; ++j)
    {
      for(int k = 0; k < bPolygon2.numEdges(); ++k)
      {
        bPolygon2[k][j][i] += .11;
      }
    }
  }

  std::vector<CurvedPolygonType> bPolygons3;
  bool didIntersect = intersect_polygon(bPolygon, bPolygon2, bPolygons3);
  EXPECT_TRUE(didIntersect);

  CoordType A = 0.0;
  for(int i = 0; i < static_cast<int>(bPolygons3.size()); ++i)
  {
    A += bPolygons3[i].area(1e-14);
  }
  CoordType expA = -0.0793347222222222222;
  EXPECT_NEAR(A, expA, 1e-14);
}

//----------------------------------------------------------------------------------
TEST(primal_curvedpolygon, area_intersection_triangle_quadratic)
{
  const int DIM = 2;
  const int order = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  SLIC_INFO(
    "Test intersecting two quadratic triangular CurvedPolygons (single "
    "region).");

  CurvedPolygonType bPolygon;
  EXPECT_EQ(0, bPolygon.numEdges());

  PointType controlPoints[order + 1] = {PointType::make_point(0.6, 1.2),
                                        PointType::make_point(0.4, 1.3),
                                        PointType::make_point(0.3, 2.0)};

  PointType controlPoints2[order + 1] = {PointType::make_point(0.3, 2.0),
                                         PointType::make_point(0.27, 1.5),
                                         PointType::make_point(0.0, 1.6)};

  PointType controlPoints3[order + 1] = {PointType::make_point(0.0, 1.6),
                                         PointType::make_point(0.1, 1.5),
                                         PointType::make_point(0.6, 1.2)};

  BezierCurveType bCurve(controlPoints, order);
  bPolygon.addEdge(bCurve);

  BezierCurveType bCurve2(controlPoints2, order);
  bPolygon.addEdge(bCurve2);

  BezierCurveType bCurve3(controlPoints3, order);
  bPolygon.addEdge(bCurve3);
  CurvedPolygonType bPolygon2 = bPolygon;
  for(int i = 0; i < DIM; ++i)
  {
    for(int j = 0; j <= order; ++j)
    {
      for(int k = 0; k < bPolygon2.numEdges(); ++k)
      {
        bPolygon2[k][j][i] += .11;
      }
    }
  }

  std::vector<CurvedPolygonType> bPolygons3;
  bool didIntersect = intersect_polygon(bPolygon, bPolygon2, bPolygons3);
  EXPECT_TRUE(didIntersect);

  CoordType A = 0.0;
  for(int i = 0; i < static_cast<int>(bPolygons3.size()); ++i)
  {
    A += bPolygons3[i].area(1e-14);
  }
  CoordType expA = -0.024649833203616;
  EXPECT_NEAR(A, expA, 1e-14);
}

TEST(primal_curvedpolygon, area_intersection_triangle_quadratic_two_regions)
{
  const int DIM = 2;
  const int order = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  SLIC_INFO(
    "Test intersecting two quadratic triangular CurvedPolygons (two regions).");

  CurvedPolygonType bPolygon;
  EXPECT_EQ(0, bPolygon.numEdges());

  PointType controlPoints[order + 1] = {PointType::make_point(0.6, 1.2),
                                        PointType::make_point(0.4, 1.3),
                                        PointType::make_point(0.3, 2.0)};

  PointType controlPoints2[order + 1] = {PointType::make_point(0.3, 2.0),
                                         PointType::make_point(0.27, 1.5),
                                         PointType::make_point(0.0, 1.6)};

  PointType controlPoints3[order + 1] = {PointType::make_point(0.0, 1.6),
                                         PointType::make_point(0.1, 1.5),
                                         PointType::make_point(0.6, 1.2)};

  PointType controlPoints4[order + 1] = {PointType::make_point(1.0205, 1.6699),
                                         PointType::make_point(0.8339, 1.5467),
                                         PointType::make_point(0.1777, 1.8101)};

  PointType controlPoints5[order + 1] = {PointType::make_point(0.1777, 1.8101),
                                         PointType::make_point(0.5957, 1.5341),
                                         PointType::make_point(0.3741, 1.3503)};

  PointType controlPoints6[order + 1] = {PointType::make_point(0.3741, 1.3503),
                                         PointType::make_point(0.5107, 1.3869),
                                         PointType::make_point(1.0205, 1.6699)};
  CurvedPolygonType bPolygon2;

  BezierCurveType bCurve(controlPoints, order);
  BezierCurveType bCurve4(controlPoints4, order);
  bPolygon2.addEdge(bCurve4);
  bPolygon.addEdge(bCurve);

  BezierCurveType bCurve2(controlPoints2, order);
  BezierCurveType bCurve5(controlPoints5, order);
  bPolygon2.addEdge(bCurve5);
  bPolygon.addEdge(bCurve2);

  BezierCurveType bCurve3(controlPoints3, order);
  BezierCurveType bCurve6(controlPoints6, order);
  bPolygon2.addEdge(bCurve6);
  bPolygon.addEdge(bCurve3);

  std::vector<CurvedPolygonType> bPolygons3;
  bool didIntersect = intersect_polygon(bPolygon, bPolygon2, bPolygons3);
  EXPECT_TRUE(didIntersect);
  EXPECT_EQ(bPolygons3.size(), 2);
}

//----------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}

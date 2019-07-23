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
TEST(primal_beziercurve, area_triangle_linear)
{
  const int DIM = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  SLIC_INFO("Test checking CurvedPolygon linear area triangle computation.");

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

  EXPECT_DOUBLE_EQ(trueA, A);
}

//----------------------------------------------------------------------------------
TEST(primal_beziercurve, area_triangle_quadratic)
{
  const int DIM = 2;
  const int order = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  SLIC_INFO("Test checking CurvedPolygon linear area triangle computation.");

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
TEST(primal_beziercurve, area_triangle_mixed_order)
{
  const int DIM = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  SLIC_INFO("Test checking CurvedPolygon linear area triangle computation.");

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
TEST(primal_beziercurve, split_edge)
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

  bPolygon.splitEdge(0, .5);
  bCurve.split(.5, bCurve2, bCurve3);

  CurvedPolygonType bPolygon2 = bPolygon;
  std::vector<CurvedPolygonType> bPolygon3;

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

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}

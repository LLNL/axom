// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/* /file primal_curved_polygon.cpp
 * /brief This file tests the CurvedPolygon class
 */

#include "gtest/gtest.h"

#include "axom/slic.hpp"
#include "axom/primal/geometry/CurvedPolygon.hpp"

namespace primal = axom::primal;

/**
 * Helper function to compute the area and centroid of a curved polygon and to check that they match expectations, stored in \a expArea and \a expCentroid. Areas and Moments are computed within tolerance \a eps and checks use \a test_eps.
 */
template <typename CoordType, int DIM>
void checkMoments(const primal::CurvedPolygon<CoordType, DIM>& bPolygon,
                  const CoordType expArea,
                  const primal::Point<CoordType, DIM>& expMoment,
                  double eps,
                  double test_eps)
{
  using Array = std::vector<CoordType>;

  EXPECT_DOUBLE_EQ(expArea, bPolygon.area(eps));
  for(int i = 0; i < DIM; ++i)
  {
    EXPECT_DOUBLE_EQ(expMoment[i], bPolygon.centroid(eps)[i]);
  }
}

template <typename CoordType, int DIM>
primal::CurvedPolygon<CoordType, DIM> createPolygon(
  std::vector<primal::Point<CoordType, DIM>> ControlPoints,
  std::vector<int> orders)
{
  using PointType = primal::Point<CoordType, DIM>;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  const int num_edges = orders.size();
  const int num_unique_control_points = ControlPoints.size();

  //std::cout << num_edges << ", " << num_unique_control_points << std::endl;
  //std::cout << ControlPoints << std::endl;
  for(int i = 0; i < num_edges; ++i)
  {
    std::cout << orders[i] << std::endl;
  }

  //checks if the orders and control points given will give a valid polygon
  EXPECT_EQ(accumulate(orders.begin(), orders.end(), 0) + 1,
            num_unique_control_points);

  CurvedPolygonType bPolygon;
  int iter = 0;
  for(int j = 0; j < num_edges; ++j)
  {
    std::vector<PointType> subCP;
    subCP.assign(ControlPoints.begin() + iter,
                 ControlPoints.begin() + iter + orders[j] + 1);
    BezierCurveType addCurve(subCP, orders[j]);
    bPolygon.addEdge(addCurve);
    iter += (orders[j]);
  }
  std::cout << bPolygon << std::endl;
  return bPolygon;
}

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

  SLIC_INFO("Test checking if CurvedPolygon is closed.");

  CurvedPolygonType bPolygon;
  EXPECT_EQ(0, bPolygon.numEdges());
  EXPECT_EQ(false, bPolygon.isClosed());

  std::vector<PointType> CP = {PointType::make_point(0.6, 1.2),
                               PointType::make_point(0.3, 2.0),
                               PointType::make_point(0.0, 1.6),
                               PointType::make_point(0.6, 1.2)};
  std::vector<int> orders = {1, 1, 1};

  std::vector<PointType> subCP = {PointType::make_point(0.6, 1.2),
                                  PointType::make_point(0.3, 2.0)};
  std::vector<int> suborders = {1};
  CurvedPolygonType subPolygon = createPolygon(subCP, suborders);
  EXPECT_EQ(false, subPolygon.isClosed());

  bPolygon = createPolygon(CP, orders);

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

  std::vector<PointType> CP = {PointType::make_point(0.6, 1.2),
                               PointType::make_point(0.3, 2.0),
                               PointType::make_point(0.0, 1.6),
                               PointType::make_point(0.6, 1.2)};

  std::vector<int> orders32 = {1, 1, 1};
  CurvedPolygonType bPolygon32 = createPolygon(CP, orders32);
  std::cout << "Got here!! " << std::endl;
  std::vector<PointType> subCP;

  subCP.assign(CP.begin(), CP.begin() + 2);
  BezierCurveType bCurve(subCP, 1);
  bPolygon32.splitEdge(0, .5);

  BezierCurveType bCurve2;
  BezierCurveType bCurve3;
  bCurve.split(.5, bCurve2, bCurve3);

  EXPECT_EQ(bPolygon32.numEdges(), 4);
  for(int i = 0; i < bPolygon32[0].getOrder(); ++i)
  {
    for(int dimi = 0; dimi < DIM; ++dimi)
    {
      EXPECT_EQ(bPolygon32[0][i][dimi], bCurve2[i][dimi]);
      EXPECT_EQ(bPolygon32[1][i][dimi], bCurve3[i][dimi]);
    }
  }
}

//----------------------------------------------------------------------------------
TEST(primal_curvedpolygon, moments_triangle_degenerate)
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
  PointType origin = PointType::make_point(0.0, 0.0);

  PointType controlPoints[2] = {PointType::make_point(0.6, 1.2),
                                PointType::make_point(0.3, 2.0)};

  PointType controlPoints2[2] = {PointType::make_point(0.3, 2.0),
                                 PointType::make_point(0.0, 1.6)};

  PointType controlPoints3[2] = {PointType::make_point(0.0, 1.6),
                                 PointType::make_point(0.6, 1.2)};

  BezierCurveType bCurve(controlPoints, 1);
  bPolygon.addEdge(bCurve);
  checkMoments(bPolygon, 0.0, origin, 1e-14, 1e-15);

  BezierCurveType bCurve2(controlPoints2, 1);
  bPolygon.addEdge(bCurve2);
  checkMoments(bPolygon, 0.0, origin, 1e-14, 1e-15);

  BezierCurveType bCurve3(controlPoints3, 1);
  bPolygon.addEdge(bCurve3);

  bPolygon[2][1][0] -= 1e-11;
  checkMoments(bPolygon, 0.0, origin, 1e-14, 1e-15);
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

  PointType M = bPolygon.centroid();
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

  PointType M = bPolygon.centroid();
  CoordType trueM2 = 1.55764705882353;
  CoordType trueM1 = .2970147058823527;

  EXPECT_DOUBLE_EQ(trueM1, M[0]);
  EXPECT_DOUBLE_EQ(trueM2, M[1]);
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

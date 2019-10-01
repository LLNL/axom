// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/* /file primal_curved_polygon.cpp
 * /brief This file tests the CurvedPolygon class and associated intersections
 */

#include "gtest/gtest.h"

#include "axom/slic.hpp"
#include "axom/primal/geometry/CurvedPolygon.hpp"
#include "axom/primal/operators/intersect.hpp"

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
  EXPECT_NEAR(expArea, bPolygon.area(eps), test_eps);
  for(int i = 0; i < DIM; ++i)
  {
    EXPECT_NEAR(expMoment[i], bPolygon.centroid(eps)[i], test_eps);
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
TEST(primal_curvedpolygon, moments_triangle_linear)
{
  const int DIM = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;

  SLIC_INFO("Test checking CurvedPolygon linear triangle moment computation.");
  std::vector<PointType> CP = {PointType::make_point(0.6, 1.2),
                               PointType::make_point(0.3, 2.0),
                               PointType::make_point(0.0, 1.6),
                               PointType::make_point(0.6, 1.2)};

  std::vector<int> orders = {1, 1, 1};
  CurvedPolygonType bPolygon = createPolygon(CP, orders);

  CoordType trueA = -.18;
  PointType trueC = PointType::make_point(0.3, 1.6);

  checkMoments(bPolygon, trueA, trueC, 1e-14, 1e-15);
}

//----------------------------------------------------------------------------------
TEST(primal_curvedpolygon, moments_triangle_quadratic)
{
  const int DIM = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;

  SLIC_INFO("Test checking CurvedPolygon quadratic triangle area computation.");

  SLIC_INFO(
    "Test checking CurvedPolygon mixed order triangle moment computation.");

  std::vector<PointType> CP = {PointType::make_point(0.6, 1.2),
                               PointType::make_point(0.4, 1.3),
                               PointType::make_point(0.3, 2.0),
                               PointType::make_point(0.27, 1.5),
                               PointType::make_point(0.0, 1.6),
                               PointType::make_point(0.1, 1.5),
                               PointType::make_point(0.6, 1.2)};

  std::vector<int> orders = {2, 2, 2};
  CurvedPolygonType bPolygon = createPolygon(CP, orders);

  CoordType trueA = -0.097333333333333;
  PointType trueC = PointType::make_point(.294479452054794, 1.548219178082190);

  checkMoments(bPolygon, trueA, trueC, 1e-15, 1e-14);
}

//----------------------------------------------------------------------------------
TEST(primal_curvedpolygon, moments_triangle_mixed_order)
{
  const int DIM = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;

  SLIC_INFO(
    "Test checking CurvedPolygon mixed order triangle moment computation.");

  std::vector<PointType> CP = {PointType::make_point(0.6, 1.2),
                               PointType::make_point(0.4, 1.3),
                               PointType::make_point(0.3, 2.0),
                               PointType::make_point(0.27, 1.5),
                               PointType::make_point(0.0, 1.6),
                               PointType::make_point(0.6, 1.2)};

  std::vector<int> orders = {2, 2, 1};
  CurvedPolygonType bPolygon = createPolygon(CP, orders);

  CoordType trueA = -.0906666666666666666666;
  PointType trueC = PointType::make_point(.2970147058823527, 1.55764705882353);

  checkMoments(bPolygon, trueA, trueC, 1e-14, 1e-15);
}

//----------------------------------------------------------------------------------
TEST(primal_curvedpolygon, moments_quad_all_orders)
{
  const int DIM = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;

  SLIC_INFO("Test checking CurvedPolygon linear triangle moment computation.");
  std::vector<PointType> CPorig = {PointType::make_point(0.0, 0.0),
                                   PointType::make_point(0.0, 1.0),
                                   PointType::make_point(1.0, 1.0),
                                   PointType::make_point(1.0, 0.0),
                                   PointType::make_point(0.0, 0.0)};

  std::vector<int> orders = {1, 1, 1, 1};
  CurvedPolygonType bPolygon = createPolygon(CPorig, orders);

  CoordType trueA = 1.0;
  PointType trueC = PointType::make_point(0.5, 0.5);

  checkMoments(bPolygon, trueA, trueC, 1e-14, 1e-15);
  for(int p = 2; p < 11; ++p)
  {
    std::vector<PointType> CP = CPorig;
    for(int side = 0; side < 4; ++side)
    {
      for(int i = 1; i < p; ++i)
      {
        switch(side)
        {
        case 0:
          CP.insert(CP.begin() + i + (side * p),
                    (PointType::make_point(0.0, 1.0 * i / p)));
          break;
        case 1:
          CP.insert(CP.begin() + i + (side * p),
                    (PointType::make_point(1.0 * i / p, 1.0)));
          break;
        case 2:
          CP.insert(CP.begin() + i + (side * p),
                    (PointType::make_point(1.0, 1.0 - (1.0 * i / p))));
          break;
        case 3:
          CP.insert(CP.begin() + i + (side * p),
                    (PointType::make_point(1.0 - (1.0 * i / p), 0.0)));
          break;
        }
      }
      orders[side] += 1;
    }
    /*for (int i=0; i<CP.size(); ++i)
  {
        std::cout << CP[i] << std::endl;
      }*/
    bPolygon = createPolygon(CP, orders);
    checkMoments(bPolygon, trueA, trueC, 1e-14, 1e-15);
  }
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

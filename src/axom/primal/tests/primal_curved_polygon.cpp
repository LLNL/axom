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

#include <numeric>

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

/* Helper function to create a CurvedPolygon from a list of control points and a list of orders of component curves. Control points should be given as a list of Points in order of orientation with no duplicates except that the first control point should also be the last control point (if the polygon is closed). Orders should be given as a list of ints in order of orientation, representing the orders of the component curves.
 */
template <typename CoordType, int DIM>
primal::CurvedPolygon<CoordType, DIM> createPolygon(
  const std::vector<primal::Point<CoordType, DIM>> ControlPoints,
  const std::vector<int> orders)
{
  using PointType = primal::Point<CoordType, DIM>;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  const int num_edges = orders.size();
  const int num_unique_control_points = ControlPoints.size();

  //checks if the orders and control points given will yield a valid curved polygon
  {
    const int sum_of_orders = std::accumulate(orders.begin(), orders.end(), 0);
    EXPECT_EQ(sum_of_orders + 1, num_unique_control_points);
  }

  //Converts the control points to BezierCurves of specified orders and stores them in a CurvedPolygon object.
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

  SLIC_INFO(
    "Test checking CurvedPolygon quadratic triangle moment computation.");

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
  bool didIntersect = intersect(bPolygon, bPolygon2, bPolygons3);
  EXPECT_TRUE(didIntersect);
  std::cout << bPolygons3.size() << std::endl;
  // int nEd= bPolygons3.size();
  for(int i = 0; i < static_cast<int>(bPolygons3.size()); ++i)
  {
    std::cout << bPolygons3[i] << std::endl;
  }
  for(int i = 0; i < DIM; ++i)
  {
    for(int j = 0; j <= order; ++j)
    {
      for(int idxcurve = 0; idxcurve < static_cast<int>(bPolygons3.size());
          ++idxcurve)
      {
        for(int k = 0; k < bPolygons3[idxcurve].numEdges(); ++k)
        {
          EXPECT_NEAR(expbPolygon[k][j][i], bPolygons3[idxcurve][k][j][i], 1e-10);
          /*EXPECT_TRUE(axom::utilities::isNearlyEqual(expbPolygon[k][j][i],bPolygons3[idxcurve][k][j][i],2.));*/
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
  bool didIntersect = intersect(bPolygon, bPolygon2, bPolygons3);
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
                                                     2.));
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
  bool didIntersect = intersect(bPolygon, bPolygon2, bPolygons3);
  EXPECT_TRUE(didIntersect);

  CoordType A = 0.0;
  for(int i = 0; i < static_cast<int>(bPolygons3.size()); ++i)
  {
    A += bPolygons3[i].area(1e-14);
  }
  CoordType expA = -0.0793347222222222222;
  EXPECT_NEAR(A, expA, 1e-10);
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
  bool didIntersect = intersect(bPolygon, bPolygon2, bPolygons3, 1e-14);
  EXPECT_TRUE(didIntersect);

  CoordType A = 0.0;
  for(int i = 0; i < static_cast<int>(bPolygons3.size()); ++i)
  {
    A += bPolygons3[i].area(1e-8);
  }
  CoordType expA = -0.024649833203616;
  EXPECT_NEAR(A, expA, 1e-10);
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
  bool didIntersect = intersect(bPolygon, bPolygon2, bPolygons3);
  EXPECT_TRUE(didIntersect);
  EXPECT_EQ(bPolygons3.size(), 2);
}

TEST(primal_curvedpolygon, area_intersection_triangle_inclusion)
{
  const int DIM = 2;
  const int order = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  SLIC_INFO(
    "Test intersecting two quadratic triangular CurvedPolygons (inclusion).");

  CurvedPolygonType bPolygon;
  EXPECT_EQ(0, bPolygon.numEdges());

  PointType controlPoints[order + 1] = {PointType::make_point(0.0, 0.0),
                                        PointType::make_point(0.5, 0.0),
                                        PointType::make_point(1.0, 0.0)};

  PointType controlPoints2[order + 1] = {PointType::make_point(1.0, 0.0),
                                         PointType::make_point(0.5, 0.5),
                                         PointType::make_point(0.0, 1.0)};

  PointType controlPoints3[order + 1] = {PointType::make_point(0.0, 1.0),
                                         PointType::make_point(0.0, 0.5),
                                         PointType::make_point(0.0, 0.0)};

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
        bPolygon2[k][j][i] = bPolygon2[k][j][i] * .5 + .05;
      }
    }
  }

  std::vector<CurvedPolygonType> bPolygons3;
  bool didIntersect = intersect(bPolygon, bPolygon2, bPolygons3);
  bPolygons3.clear();
  bool didIntersect2 = intersect(bPolygon2, bPolygon, bPolygons3);
  EXPECT_TRUE(didIntersect);
  EXPECT_TRUE(didIntersect2);
  EXPECT_EQ(bPolygons3.size(), 1);
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

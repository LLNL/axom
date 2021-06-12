// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file primal_curved_polygon.cpp
 * \brief This file tests the CurvedPolygon class
 */

#include "gtest/gtest.h"

#include "axom/slic.hpp"
#include "axom/primal/geometry/CurvedPolygon.hpp"
#include "axom/primal/operators/intersect.hpp"

#include <numeric>

namespace primal = axom::primal;

/*!
 * Helper function to compute the area and centroid of a curved polygon and to check that they match expectations, 
 * stored in \a expArea and \a expCentroid. Areas and Moments are computed within tolerance \a eps and checks use \a test_eps.
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

/*!
 * Helper function to create a CurvedPolygon from a list of control points and a list of orders of component curves. 
 * Control points should be given as a list of Points in order of orientation with no duplicates except that 
 * the first control point should also be the last control point (if the polygon is closed). 
 * Orders should be given as a list of ints in order of orientation, representing the orders of the component curves.
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

  PointType controlPoints[2] = {PointType {0.6, 1.2}, PointType {0.0, 1.6}};

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

  std::vector<PointType> CP = {PointType {0.6, 1.2},
                               PointType {0.3, 2.0},
                               PointType {0.0, 1.6},
                               PointType {0.6, 1.2}};
  std::vector<int> orders = {1, 1, 1};

  std::vector<PointType> subCP = {PointType {0.6, 1.2}, PointType {0.3, 2.0}};
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

  std::vector<PointType> CP = {PointType {0.6, 1.2},
                               PointType {0.3, 2.0},
                               PointType {0.0, 1.6},
                               PointType {0.6, 1.2}};

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

  SLIC_INFO("Testing area computation of degenerate triangles");

  CurvedPolygonType bPolygon;
  EXPECT_EQ(0, bPolygon.numEdges());
  EXPECT_EQ(0.0, bPolygon.area());
  PointType origin = PointType::make_point(0.0, 0.0);

  PointType controlPoints[2] = {PointType {0.6, 1.2}, PointType {0.3, 2.0}};

  PointType controlPoints2[2] = {PointType {0.3, 2.0}, PointType {0.0, 1.6}};
  PointType controlPoints3[2] = {PointType {0.0, 1.6}, PointType {0.6, 1.2}};

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

  SLIC_INFO("Test moment computation of a linear triangle");
  std::vector<PointType> CP = {PointType {0.6, 1.2},
                               PointType {0.3, 2.0},
                               PointType {0.0, 1.6},
                               PointType {0.6, 1.2}};

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

  SLIC_INFO("Test moment computation of quadratic triangle");

  std::vector<PointType> CP = {PointType {0.6, 1.2},
                               PointType {0.4, 1.3},
                               PointType {0.3, 2.0},
                               PointType {0.27, 1.5},
                               PointType {0.0, 1.6},
                               PointType {0.1, 1.5},
                               PointType {0.6, 1.2}};

  std::vector<int> orders = {2, 2, 2};
  CurvedPolygonType bPolygon = createPolygon(CP, orders);

  CoordType trueA = -0.097333333333333;
  PointType trueC {.294479452054794, 1.548219178082190};

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
    "Test moment computation for curved triangle with mixed order edges");

  std::vector<PointType> CP = {PointType {0.6, 1.2},
                               PointType {0.4, 1.3},
                               PointType {0.3, 2.0},
                               PointType {0.27, 1.5},
                               PointType {0.0, 1.6},
                               PointType {0.6, 1.2}};

  std::vector<int> orders = {2, 2, 1};
  CurvedPolygonType bPolygon = createPolygon(CP, orders);

  CoordType trueA = -.0906666666666666666666;
  PointType trueC {.2970147058823527, 1.55764705882353};

  checkMoments(bPolygon, trueA, trueC, 1e-14, 1e-15);
}

//----------------------------------------------------------------------------------
TEST(primal_curvedpolygon, moments_quad_all_orders)
{
  const int DIM = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;

  SLIC_INFO("Test moment computation for quads of different orders");
  std::vector<PointType> CPorig = {PointType {0.0, 0.0},
                                   PointType {0.0, 1.0},
                                   PointType {1.0, 1.0},
                                   PointType {1.0, 0.0},
                                   PointType {0.0, 0.0}};

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
        const int offset = i + (side * p);
        const double t = static_cast<double>(i) / p;
        switch(side)
        {
        case 0:
          CP.insert(CP.begin() + offset, PointType {0., t});
          break;
        case 1:
          CP.insert(CP.begin() + offset, PointType {t, 1.});
          break;
        case 2:
          CP.insert(CP.begin() + offset, PointType {1., 1. - t});
          break;
        case 3:
          CP.insert(CP.begin() + offset, PointType {1. - t, 0.});
          break;
        }
      }
      orders[side] += 1;
    }
    // for (int i=0; i<CP.size(); ++i)
    // {
    //   std::cout << CP[i] << std::endl;
    // }

    bPolygon = createPolygon(CP, orders);
    checkMoments(bPolygon, trueA, trueC, 1e-14, 1e-15);
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

  SLIC_INFO("Test intersection area of two linear triangles (single region)");

  CurvedPolygonType bPolygon;
  EXPECT_EQ(0, bPolygon.numEdges());

  PointType controlPoints[order + 1] = {PointType {0.6, 1.2},
                                        PointType {0.3, 2.0}};

  PointType controlPoints2[order + 1] = {PointType {0.3, 2.0},
                                         PointType {0.0, 1.6}};

  PointType controlPoints3[order + 1] = {PointType {0.0, 1.6},
                                         PointType {0.6, 1.2}};

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
    "Test intersection area of two quadratic triangles (single region)");

  CurvedPolygonType bPolygon;
  EXPECT_EQ(0, bPolygon.numEdges());

  PointType controlPoints[order + 1] = {PointType {0.6, 1.2},
                                        PointType {0.4, 1.3},
                                        PointType {0.3, 2.0}};

  PointType controlPoints2[order + 1] = {PointType {0.3, 2.0},
                                         PointType {0.27, 1.5},
                                         PointType {0.0, 1.6}};

  PointType controlPoints3[order + 1] = {PointType {0.0, 1.6},
                                         PointType {0.1, 1.5},
                                         PointType {0.6, 1.2}};

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

//----------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}

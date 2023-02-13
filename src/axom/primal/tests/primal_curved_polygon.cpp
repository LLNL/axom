// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file primal_curved_polygon.cpp
 * \brief This file tests the CurvedPolygon class
 */

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/slic.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/CurvedPolygon.hpp"
#include "axom/primal/geometry/OrientationResult.hpp"
#include "axom/primal/operators/intersect.hpp"
#include "axom/primal/operators/compute_moments.hpp"
#include "axom/primal/operators/orientation.hpp"

#include <numeric>
#include <math.h>

namespace primal = axom::primal;

/*!
 * Helper function to compute the area and centroid of a curved polygon and to check that they match expectations,
 * stored in \a expArea and \a expCentroid. Areas and Moments are computed within tolerance \a eps and checks use \a test_eps.
 */
template <typename CoordType>
void checkMoments(const primal::CurvedPolygon<CoordType, 2>& bPolygon,
                  const CoordType expArea,
                  const primal::Point<CoordType, 2>& expMoment,
                  double eps,
                  double test_eps)
{
  EXPECT_NEAR(expArea, primal::area(bPolygon, eps), test_eps);

  const auto centroid = primal::centroid(bPolygon, eps);
  for(int i = 0; i < 2; ++i)
  {
    EXPECT_NEAR(expMoment[i], centroid[i], test_eps);
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
  const axom::Array<primal::Point<CoordType, DIM>> ControlPoints,
  const axom::Array<int> orders)
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
    axom::Array<PointType> subCP(orders[j] + 1);
    for(int i = 0; i < orders[j] + 1; i++) subCP[i] = ControlPoints[i + iter];

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
    EXPECT_EQ(axom::Array<BezierCurveType>(), bPolygon.getEdges());
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

  {
    CurvedPolygonType bPolygon;
    EXPECT_EQ(0, bPolygon.numEdges());
    EXPECT_FALSE(bPolygon.isClosed());
  }

  axom::Array<PointType> CP = {PointType {0.6, 1.2},
                               PointType {0.3, 2.0},
                               PointType {0.0, 1.6},
                               PointType {0.6, 1.2}};
  axom::Array<int> orders = {1, 1, 1};

  {
    axom::Array<PointType> subCP = {PointType {0.6, 1.2}, PointType {0.3, 2.0}};
    axom::Array<int> suborders = {1};
    CurvedPolygonType subPolygon = createPolygon(subCP, suborders);
    EXPECT_FALSE(subPolygon.isClosed());
  }

  {
    CurvedPolygonType bPolygon = createPolygon(CP, orders);
    EXPECT_EQ(3, bPolygon.numEdges());
    EXPECT_TRUE(bPolygon.isClosed());

    bPolygon[2][1][0] -= 2e-15;
    EXPECT_FALSE(bPolygon.isClosed(1e-15));
  }

  {
    CurvedPolygonType bPolygon = createPolygon(CP, orders);

    bPolygon[1][0][0] = 5;
    EXPECT_FALSE(bPolygon.isClosed(1e-15));
  }
}

//----------------------------------------------------------------------------------
TEST(primal_curvedpolygon, isClosed_BiGon)
{
  const int DIM = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;

  SLIC_INFO("Test checking if CurvedPolygon is closed for a Bi-Gon");

  CurvedPolygonType bPolygon;
  EXPECT_EQ(0, bPolygon.numEdges());
  EXPECT_FALSE(bPolygon.isClosed());

  // Bi-gon defined by a quadratic edge and a straight line
  axom::Array<PointType> CP = {PointType {0.8, .25},
                               PointType {2.0, .50},
                               PointType {0.8, .75},
                               PointType {0.8, .25}};
  axom::Array<int> orders = {2, 1};

  CurvedPolygonType poly = createPolygon(CP, orders);
  EXPECT_TRUE(poly.isClosed());

  // modify a vertex of the quadratic and check again
  CurvedPolygonType poly2 = poly;
  poly2[0][2] = PointType {0.8, 1.0};
  EXPECT_FALSE(poly2.isClosed());
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

  axom::Array<PointType> CP = {PointType {0.6, 1.2},
                               PointType {0.3, 2.0},
                               PointType {0.0, 1.6},
                               PointType {0.6, 1.2}};

  axom::Array<int> orders32 = {1, 1, 1};
  CurvedPolygonType bPolygon32 = createPolygon(CP, orders32);

  // Needs to be std::vector to use .assign
  axom::Array<PointType> subCP = {PointType {0.6, 1.2}, PointType {0.3, 2.0}};

  BezierCurveType bCurve(subCP, 1);
  //std::cout << "Got here!! " << std::endl;
  //std::cout << bPolygon32 << std::endl;
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
  EXPECT_EQ(0.0, primal::area(bPolygon));
  PointType origin {0.0, 0.0};

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
  axom::Array<PointType> CP = {PointType {0.6, 1.2},
                               PointType {0.3, 2.0},
                               PointType {0.0, 1.6},
                               PointType {0.6, 1.2}};

  axom::Array<int> orders = {1, 1, 1};
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

  axom::Array<PointType> CP = {PointType {0.6, 1.2},
                               PointType {0.4, 1.3},
                               PointType {0.3, 2.0},
                               PointType {0.27, 1.5},
                               PointType {0.0, 1.6},
                               PointType {0.1, 1.5},
                               PointType {0.6, 1.2}};

  axom::Array<int> orders = {2, 2, 2};
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

  axom::Array<PointType> CP = {PointType {0.6, 1.2},
                               PointType {0.4, 1.3},
                               PointType {0.3, 2.0},
                               PointType {0.27, 1.5},
                               PointType {0.0, 1.6},
                               PointType {0.6, 1.2}};

  axom::Array<int> orders = {2, 2, 1};
  CurvedPolygonType bPolygon = createPolygon(CP, orders);

  CoordType trueA = -.0906666666666666666666;
  PointType trueC {.2970147058823527, 1.55764705882353};

  checkMoments(bPolygon, trueA, trueC, 1e-14, 1e-14);
}

//----------------------------------------------------------------------------------
TEST(primal_curvedpolygon, moments_quad_all_orders)
{
  const int DIM = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;

  SLIC_INFO("Test moment computation for quads of different orders");
  axom::Array<PointType> CPorig = {PointType {0.0, 0.0},
                                   PointType {0.0, 1.0},
                                   PointType {1.0, 1.0},
                                   PointType {1.0, 0.0},
                                   PointType {0.0, 0.0}};

  axom::Array<int> orders = {1, 1, 1, 1};
  CurvedPolygonType bPolygon = createPolygon(CPorig, orders);

  CoordType trueA = 1.0;
  PointType trueC = PointType::make_point(0.5, 0.5);

  checkMoments(bPolygon, trueA, trueC, 1e-14, 1e-15);
  for(int p = 2; p < 11; ++p)
  {
    axom::Array<PointType> CP = CPorig;
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
TEST(primal_curvedpolygon, reverseOrientation)
{
  const int DIM = 2;
  const int order = 1;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;
  using SegmentType = primal::Segment<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;

  // Test several n-gons discretizing the unit circle
  const int MAX_SEG = 10;
  const PointType origin;
  for(int nseg = 3; nseg < MAX_SEG; ++nseg)
  {
    // Create an n-gon with line segments going CCW along the unit circle
    CurvedPolygonType poly(nseg);
    axom::Array<PointType> pts(nseg + 1);
    for(int i = 0; i < nseg; ++i)
    {
      const double theta = i / static_cast<double>(nseg) * (2. * M_PI);
      pts[i] = PointType {cos(theta), sin(theta)};
    }
    pts[nseg] = pts[0];

    for(int i = 0; i < nseg; ++i)
    {
      poly[i] = BezierCurveType(&pts[i], order);
    }
    EXPECT_TRUE(poly.isClosed());

    // Perform some checks on the polygon
    for(int i = 0; i < nseg; ++i)
    {
      // check that the end point of each segment is equal to the start of the next
      auto& currentEnd = poly[i][order];
      auto& nextStart = poly[(i + 1) % nseg][0];
      EXPECT_EQ(currentEnd, nextStart);

      // check that the orientation of segment midpoints goes in the same direction
      SegmentType seg(poly[i].evaluate(0.5), poly[(i + 1) % nseg].evaluate(0.5));
      EXPECT_EQ(primal::ON_NEGATIVE_SIDE, primal::orientation(origin, seg));
    }

    // Create a polygon with reversed orientation
    CurvedPolygonType reversed = poly;
    reversed.reverseOrientation();

    EXPECT_EQ(poly.numEdges(), reversed.numEdges());

    // Perform some checks on the reversed polygon
    for(int i = 0; i < nseg; ++i)
    {
      // check that order of each segment stayed the same
      EXPECT_EQ(poly[i].getOrder(), reversed[i].getOrder());

      // check that the end point of each segment is equal to the start of the next
      auto& currentEnd = reversed[i][order];
      auto& nextStart = reversed[(i + 1) % nseg][0];
      EXPECT_EQ(currentEnd, nextStart);

      // check that segment midpoints are oriented in same direction (opposite of origin)
      SegmentType seg(reversed[i].evaluate(0.5),
                      reversed[(i + 1) % nseg].evaluate(0.5));
      EXPECT_EQ(primal::ON_POSITIVE_SIDE, primal::orientation(origin, seg));
    }

    // Check that reversing twice yields the original
    CurvedPolygonType reversedAgain = reversed;
    reversedAgain.reverseOrientation();
    EXPECT_EQ(poly, reversedAgain);
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

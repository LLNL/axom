// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/slic.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Triangle.hpp"

#include "fmt/fmt.hpp"

#include <cmath>

using namespace axom;

TEST(primal_triangle, triangle_area_2D)
{
  static const int DIM = 2;
  static const double EPS = 1e-12;
  typedef double CoordType;
  typedef primal::Point<CoordType, DIM> QPoint;
  typedef primal::Triangle<CoordType, DIM> QTri;

  QPoint pt[3] = {
    QPoint::make_point(0, 0),
    QPoint::make_point(0, 1),
    QPoint::make_point(1, 0),
  };

  QTri tri(pt[0], pt[1], pt[2]);
  EXPECT_NEAR(tri.area(), 0.5, EPS);

  tri = QTri(pt[1], pt[2], pt[0]);
  EXPECT_NEAR(tri.area(), 0.5, EPS);

  tri = QTri(pt[2], pt[1], pt[0]);
  EXPECT_NEAR(tri.area(), 0.5, EPS);

  tri = QTri(pt[0], pt[2], pt[1]);
  EXPECT_NEAR(tri.area(), 0.5, EPS);
}

//------------------------------------------------------------------------------
TEST(primal_triangle, triangle_area_3D)
{
  static const int DIM = 3;
  static const double EPS = 1e-12;
  typedef double CoordType;
  typedef primal::Point<CoordType, DIM> QPoint;
  typedef primal::Triangle<CoordType, DIM> QTri;

  QPoint pt[4] = {
    QPoint::make_point(0, 0, 0),
    QPoint::make_point(1, 0, 0),
    QPoint::make_point(0, 1, 0),
    QPoint::make_point(0, 0, 1),
  };

  QTri tri(pt[0], pt[1], pt[2]);
  EXPECT_NEAR(tri.area(), 0.5, EPS);

  tri = QTri(pt[0], pt[2], pt[3]);
  EXPECT_NEAR(tri.area(), 0.5, EPS);

  tri = QTri(pt[0], pt[1], pt[3]);
  EXPECT_NEAR(tri.area(), 0.5, EPS);

  tri = QTri(pt[1], pt[2], pt[3]);
  EXPECT_NEAR(tri.area(), std::sqrt(3.) / 2., EPS);
}

//------------------------------------------------------------------------------
TEST(primal_triangle, triangle_physical_to_bary)
{
  static const int DIM = 3;
  static const double EPS = 1e-12;
  typedef double CoordType;
  typedef primal::Point<CoordType, DIM> QPoint;
  typedef primal::Triangle<CoordType, DIM> QTri;

  QPoint pt[3] = {
    QPoint::make_point(1, 0, 0),
    QPoint::make_point(0, 1, 0),
    QPoint::make_point(0, 0, 1),
  };

  QTri tri(pt[0], pt[1], pt[2]);

  typedef std::vector<std::pair<QPoint, QPoint>> TestVec;
  TestVec testData;

  // Test the three vertices
  testData.push_back(std::make_pair(pt[0], QPoint::make_point(1., 0., 0.)));
  testData.push_back(std::make_pair(pt[1], QPoint::make_point(0., 1., 0.)));
  testData.push_back(std::make_pair(pt[2], QPoint::make_point(0., 0., 1.)));

  // Test the three edge midpoints
  testData.push_back(std::make_pair(QPoint(0.5 * (pt[0].array() + pt[1].array())),
                                    QPoint::make_point(0.5, 0.5, 0.)));
  testData.push_back(std::make_pair(QPoint(0.5 * (pt[0].array() + pt[2].array())),
                                    QPoint::make_point(0.5, 0., 0.5)));
  testData.push_back(std::make_pair(QPoint(0.5 * (pt[1].array() + pt[2].array())),
                                    QPoint::make_point(0., 0.5, 0.5)));

  // Test the triangle midpoint
  testData.push_back(std::make_pair(
    QPoint(1. / 3. * (pt[0].array() + pt[1].array() + pt[2].array())),
    QPoint::make_point(1. / 3., 1. / 3., 1. / 3.)));

  // Test a point outside the triangle
  testData.push_back(std::make_pair(
    QPoint(-0.4 * pt[0].array() + 1.2 * pt[1].array() + 0.2 * pt[2].array()),
    QPoint::make_point(-0.4, 1.2, 0.2)));

  // Now run the actual tests
  for(TestVec::const_iterator it = testData.begin(); it != testData.end(); ++it)
  {
    const QPoint& query = it->first;
    const QPoint& expBary = it->second;
    QPoint bary = tri.physToBarycentric(query);
    QPoint phys = tri.baryToPhysical(bary);

    SLIC_DEBUG(fmt::format(
      "Computed barycentric coordinates for triangle {} and point {} are {}",
      tri,
      query,
      bary));
    for(int i = 0; i < 3; ++i)
    {
      EXPECT_NEAR(bary[i], expBary[i], EPS);
      EXPECT_NEAR(phys[i], query[i], EPS);
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_triangle, triangle_bary_to_physical)
{
  static const int DIM = 3;
  static const double EPS = 1e-12;
  typedef double CoordType;
  typedef primal::Point<CoordType, DIM> QPoint;
  typedef primal::Triangle<CoordType, DIM> QTri;

  QPoint pt[3] = {
    QPoint::make_point(1, 0, 0),
    QPoint::make_point(0, 1, 0),
    QPoint::make_point(0, 0, 1),
  };

  QTri tri(pt[0], pt[1], pt[2]);

  typedef std::vector<std::pair<QPoint, QPoint>> TestVec;
  TestVec testData;

  // Test the three vertices
  testData.push_back(std::make_pair(QPoint::make_point(1., 0., 0.), pt[0]));
  testData.push_back(std::make_pair(QPoint::make_point(0., 1., 0.), pt[1]));
  testData.push_back(std::make_pair(QPoint::make_point(0., 0., 1.), pt[2]));

  // Test the three edge midpoints
  testData.push_back(
    std::make_pair(QPoint::make_point(0.5, 0.5, 0.),
                   QPoint(0.5 * (pt[0].array() + pt[1].array()))));
  testData.push_back(
    std::make_pair(QPoint::make_point(0.5, 0., 0.5),
                   QPoint(0.5 * (pt[0].array() + pt[2].array()))));
  testData.push_back(
    std::make_pair(QPoint::make_point(0., 0.5, 0.5),
                   QPoint(0.5 * (pt[1].array() + pt[2].array()))));

  // Test the triangle midpoint
  testData.push_back(std::make_pair(
    QPoint::make_point(1. / 3., 1. / 3., 1. / 3.),
    QPoint(1. / 3. * (pt[0].array() + pt[1].array() + pt[2].array()))));

  // Test a point outside the triangle
  testData.push_back(std::make_pair(
    QPoint::make_point(-0.4, 1.2, 0.2),
    QPoint(-0.4 * pt[0].array() + 1.2 * pt[1].array() + 0.2 * pt[2].array())));

  // Now run the actual tests
  for(TestVec::const_iterator it = testData.begin(); it != testData.end(); ++it)
  {
    const QPoint& query = it->first;
    const QPoint& expWorld = it->second;
    QPoint phys = tri.baryToPhysical(query);
    QPoint bary = tri.physToBarycentric(phys);

    SLIC_DEBUG(fmt::format(
      "Computed physical coordinates for triangle {} at barycentric {} are {}",
      tri,
      query,
      phys));

    for(int i = 0; i < 3; ++i)
    {
      EXPECT_NEAR(phys[i], expWorld[i], EPS);
      EXPECT_NEAR(bary[i], query[i], EPS);
    }
  }
}

//-----------------------------------------------------------------------------
TEST(primal_triangle, triangle_2D_point_containment)
{
  static const int DIM = 2;
  static const double EPS = 1e-12;
  typedef double CoordType;
  typedef primal::Point<CoordType, DIM> QPoint;
  typedef primal::Triangle<CoordType, DIM> QTri;

  // Test triangle
  QPoint pt[3] = {
    QPoint::make_point(1, 0),
    QPoint::make_point(1, 1),
    QPoint::make_point(0, 0),
  };

  QTri tri(pt[0], pt[1], pt[2]);

  typedef std::vector<QPoint> TestVec;
  TestVec successes, failures;

  // Tests that should succeed:
  // Test the three vertices
  successes.push_back(pt[0]);
  successes.push_back(pt[1]);
  successes.push_back(pt[2]);
  // Test points on the edges
  successes.push_back(QPoint::make_point(0.3, 0.3));
  successes.push_back(QPoint::make_point(0.5, 0.0));
  successes.push_back(QPoint::make_point(1.0, 0.7));
  // Test some points in the interior
  successes.push_back(QPoint::make_point(0.2, 0.15));
  successes.push_back(QPoint::make_point(0.6, 0.3));

  // Tests that should fail:
  // Point not coplanar with tri (only applicable in 3D)
  // Points outside triangle boundaries
  failures.push_back(QPoint::make_point(1, 1.01));
  failures.push_back(QPoint::make_point(50, 1000));
  // Points very close to vertices
  failures.push_back(QPoint::make_point(1.00001, 1.000001));

  // Actually run the tests
  for(TestVec::const_iterator it = successes.begin(); it != successes.end(); ++it)
  {
    EXPECT_TRUE(tri.checkInTriangle(*it, EPS));
  }
  for(TestVec::const_iterator it = failures.begin(); it != failures.end(); ++it)
  {
    EXPECT_FALSE(tri.checkInTriangle(*it, EPS));
  }
}

//------------------------------------------------------------------------------
TEST(primal_triangle, triangle_3D_point_containment)
{
  static const int DIM = 3;
  static const double EPS = 1e-12;
  typedef double CoordType;
  typedef primal::Point<CoordType, DIM> QPoint;
  typedef primal::Triangle<CoordType, DIM> QTri;

  // Test triangle
  QPoint pt[3] = {
    QPoint::make_point(1, 0, 0),
    QPoint::make_point(1, 1, 0),
    QPoint::make_point(0, 0, 0),
  };

  QTri tri(pt[0], pt[1], pt[2]);

  typedef std::vector<QPoint> TestVec;
  TestVec successes, failures;

  // Tests that should succeed:
  // Test the three vertices
  successes.push_back(pt[0]);
  successes.push_back(pt[1]);
  successes.push_back(pt[2]);
  // Test points on the edges
  successes.push_back(QPoint::make_point(0.3, 0.3, 0));
  successes.push_back(QPoint::make_point(0.5, 0.0, 0));
  successes.push_back(QPoint::make_point(1.0, 0.7, 0));
  // Test some points in the interior
  successes.push_back(QPoint::make_point(0.2, 0.15, 0));
  successes.push_back(QPoint::make_point(0.6, 0.3, 0));

  // Tests that should fail:
  // Point not coplanar with tri (only applicable in 3D)
  failures.push_back(QPoint::make_point(0.2, 0.15, 0.00001));
  failures.push_back(QPoint::make_point(0.6, 0.3, 0.1));
  failures.push_back(QPoint::make_point(0.9999, 0.99, -0.0000001));
  // Points outside triangle boundaries
  failures.push_back(QPoint::make_point(1, 1.01, 0));
  failures.push_back(QPoint::make_point(50, 1000, 0));
  // Points very close to vertices
  failures.push_back(QPoint::make_point(1.00001, 1.000001, 0));

  // Actually run the tests
  for(TestVec::const_iterator it = successes.begin(); it != successes.end(); ++it)
  {
    EXPECT_TRUE(tri.checkInTriangle(*it, EPS));
  }
  for(TestVec::const_iterator it = failures.begin(); it != failures.end(); ++it)
  {
    EXPECT_FALSE(tri.checkInTriangle(*it, EPS));
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
using axom::slic::SimpleLogger;

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);

  SimpleLogger logger;  // create & initialize test logger,
  axom::slic::setLoggingMsgLevel(axom::slic::message::Info);

  int result = RUN_ALL_TESTS();
  return result;
}

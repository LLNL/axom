// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <limits>
#include <algorithm>

#include "gtest/gtest.h"
#include "axom/primal.hpp"

namespace primal = axom::primal;

//------------------------------------------------------------------------------
TEST(primal_closest_point, obb_test_closest_point_interior)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QVector = primal::Vector<CoordType, DIM>;
  using QOBBox = primal::OrientedBoundingBox<CoordType, DIM>;

  constexpr double ONE_OVER_SQRT_TWO = 0.7071;
  QPoint pt1;      // origin
  QVector u[DIM];  // make axes
  u[0][0] = ONE_OVER_SQRT_TWO;
  u[0][1] = ONE_OVER_SQRT_TWO;
  u[1][0] = ONE_OVER_SQRT_TWO;
  u[1][1] = -ONE_OVER_SQRT_TWO;
  u[2][2] = 1.;

  QVector e = QVector(1.);

  QOBBox obbox1(pt1, u, e);

  // interior point
  EXPECT_TRUE(primal::closest_point(pt1, obbox1) == pt1);
}

//------------------------------------------------------------------------------
TEST(primal_closest_point, obb_test_closest_point_vertex)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QVector = primal::Vector<CoordType, DIM>;
  using QOBBox = primal::OrientedBoundingBox<CoordType, DIM>;

  constexpr double ONE_OVER_SQRT_TWO = 0.7071;
  QPoint pt1;      // origin
  QVector u[DIM];  // make axes
  u[0][0] = ONE_OVER_SQRT_TWO;
  u[0][1] = ONE_OVER_SQRT_TWO;
  u[1][0] = ONE_OVER_SQRT_TWO;
  u[1][1] = -ONE_OVER_SQRT_TWO;
  u[2][2] = 1.;

  QVector e = QVector(1.);
  QOBBox obbox1(pt1, u, e);
  std::vector<QPoint> verts = obbox1.vertices();

  QPoint pt2 = verts[0];
  QPoint pt3(10. * pt2.array());

  // closest point is a vertex
  EXPECT_TRUE(primal::closest_point(pt3, obbox1) == pt2);
}

//------------------------------------------------------------------------------
TEST(primal_closest_point, obb_test_closest_point_face)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QVector = primal::Vector<CoordType, DIM>;
  using QOBBox = primal::OrientedBoundingBox<CoordType, DIM>;

  constexpr double ONE_OVER_SQRT_TWO = 0.7071;
  constexpr double EPS = 0.01;
  QPoint pt1;      // origin
  QVector u[DIM];  // make standard axes
  u[0][0] = ONE_OVER_SQRT_TWO;
  u[0][1] = ONE_OVER_SQRT_TWO;
  u[1][0] = ONE_OVER_SQRT_TWO;
  u[1][1] = -ONE_OVER_SQRT_TWO;
  u[2][2] = 1.;

  QVector e = QVector(1.);
  QOBBox obbox1(pt1, u, e);

  QPoint pt2(u[0].array());
  QPoint pt3(10. * pt2.array());

  QVector found(primal::closest_point(pt3, obbox1));
  QVector expected(pt2);

  // closest point is on a face
  EXPECT_TRUE((found - expected).squared_norm() < EPS);
}

//------------------------------------------------------------------------------
TEST(primal_closest_point, obb_test_closest_point_edge)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QVector = primal::Vector<CoordType, DIM>;
  using QOBBox = primal::OrientedBoundingBox<CoordType, DIM>;

  constexpr double ONE_OVER_SQRT_TWO = 0.7071;
  constexpr double EPS = 0.01;
  QPoint pt1;      // origin
  QVector u[DIM];  // make standard axes
  u[0][0] = ONE_OVER_SQRT_TWO;
  u[0][1] = ONE_OVER_SQRT_TWO;
  u[1][0] = ONE_OVER_SQRT_TWO;
  u[1][1] = -ONE_OVER_SQRT_TWO;
  u[2][2] = 1.;

  QVector e = QVector(1.);
  QOBBox obbox1(pt1, u, e);

  QPoint pt2(u[0].array() + u[1].array());
  QPoint pt3(10. * pt2.array());

  QVector found(primal::closest_point(pt3, obbox1));
  QVector expected(pt2);

  // closest point is on an edge
  EXPECT_TRUE((found - expected).squared_norm() < EPS);
}

//------------------------------------------------------------------------------
TEST(primal_closest_point, sphere_point_2D)
{
  constexpr int DIM = 2;
  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QVector = primal::Vector<CoordType, DIM>;
  using QSphere = primal::Sphere<CoordType, DIM>;

  constexpr double EPS = 1e-12;

  // define several test points
  axom::Array<QPoint> test_points = {
    QPoint {1, 0},
    QPoint {0, 1},
    QPoint {-1, 0},
    QPoint {0, -1},
    QPoint {1, 1},
    QPoint {2, 2},
    QPoint {5.234, -6.432},
  };

  // check test points against unit sphere
  for(const auto& pt : test_points)
  {
    auto unit_sphere = QSphere(1.);
    auto cp = primal::closest_point(pt, unit_sphere);
    auto exp_cp = unit_sphere.getCenter() + QVector(pt).unitVector();

    EXPECT_NEAR(exp_cp[0], cp[0], EPS);
    EXPECT_NEAR(exp_cp[1], cp[1], EPS);
    EXPECT_NEAR(0., primal::squared_distance(exp_cp, cp), EPS);
  }

  // define several spheres of different centers and radii
  axom::Array<double> radii = {1., .12345, 543.21, .4, 0.001, 87, 0.};
  axom::Array<QPoint> centers = {QPoint {0, 0},
                                 QPoint {1, 0},
                                 QPoint {7.5, 8.25},
                                 QPoint {-1.3, 2.3},
                                 QPoint {-20, -40},
                                 QPoint {4.3, 0},
                                 QPoint {0, -3.4}};

  // add sphere centers to test_points
  for(const auto& pt : centers)
  {
    test_points.push_back(pt);
  }

  // check that closest_point lies on the sphere
  for(const auto& pt : test_points)
  {
    for(const auto r : radii)
    {
      for(const auto& ctr : centers)
      {
        auto sphere = QSphere(ctr, r);
        auto cp = primal::closest_point(pt, sphere);
        EXPECT_NEAR(sphere.getRadius(),
                    sqrt(primal::squared_distance(sphere.getCenter(), cp)),
                    EPS);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_closest_point, sphere_point_3D)
{
  constexpr int DIM = 3;
  using CoordType = double;
  using QPoint = primal::Point<CoordType, DIM>;
  using QVector = primal::Vector<CoordType, DIM>;
  using QSphere = primal::Sphere<CoordType, DIM>;

  constexpr double EPS = 1e-12;

  // define several test points
  axom::Array<QPoint> test_points = {
    QPoint {1, 0, 0},
    QPoint {0, 1, 0},
    QPoint {0, 0, 1},
    QPoint {-1, 0, 1},
    QPoint {0, -1, -1},
    QPoint {1, 1, 1},
    QPoint {2, 4, 8},
    QPoint {5.234, -6.432, 7.345},
  };

  // check test points against unit sphere
  for(const auto& pt : test_points)
  {
    auto unit_sphere = QSphere(1.);
    auto cp = primal::closest_point(pt, unit_sphere);
    auto exp_cp = unit_sphere.getCenter() + QVector(pt).unitVector();

    EXPECT_NEAR(exp_cp[0], cp[0], EPS);
    EXPECT_NEAR(exp_cp[1], cp[1], EPS);
    EXPECT_NEAR(exp_cp[2], cp[2], EPS);
    EXPECT_NEAR(0., primal::squared_distance(exp_cp, cp), EPS);
  }

  // define several spheres of different centers and radii
  axom::Array<double> radii = {1., .12345, 543.21, .4, 0.001, 87, 0.};
  axom::Array<QPoint> centers = {QPoint {0, 0, 0},
                                 QPoint {1, 0, 0},
                                 QPoint {7.5, 8.25, 9.7125},
                                 QPoint {-1.3, 2.3, -3.45},
                                 QPoint {-20, -40, -80},
                                 QPoint {0, -3.4, -1}};

  // add sphere centers to test_points
  for(const auto& pt : centers)
  {
    test_points.push_back(pt);
  }

  // check that closest_point lies on the sphere
  for(const auto& pt : test_points)
  {
    for(const auto r : radii)
    {
      for(const auto& ctr : centers)
      {
        auto sphere = QSphere(ctr, r);
        auto cp = primal::closest_point(pt, sphere);
        EXPECT_NEAR(sphere.getRadius(),
                    sqrt(primal::squared_distance(sphere.getCenter(), cp)),
                    EPS);
      }
    }
  }
}
/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#include <limits>
#include <algorithm>

#include "gtest/gtest.h"

#include "primal/NumericArray.hpp"
#include "primal/Point.hpp"
#include "primal/Vector.hpp"
#include "primal/OrientedBoundingBox.hpp"
#include "primal/closest_point.hpp"

using namespace axom;

//------------------------------------------------------------------------------
TEST( primal_closest_point, obb_test_closest_point_interior )
{
  static const int DIM = 3;
  typedef double CoordType;
  typedef primal::Point< CoordType, DIM > QPoint;
  typedef primal::Vector< CoordType, DIM > QVector;
  typedef primal::OrientedBoundingBox< CoordType, DIM > QOBBox;

  const double ONE_OVER_SQRT_TWO = 0.7071;
  QPoint pt1;  // origin
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
TEST( primal_closest_point, obb_test_closest_point_vertex )
{
  static const int DIM = 3;
  typedef double CoordType;
  typedef primal::Point< CoordType, DIM > QPoint;
  typedef primal::Vector< CoordType, DIM > QVector;
  typedef primal::OrientedBoundingBox< CoordType, DIM > QOBBox;

  const double ONE_OVER_SQRT_TWO = 0.7071;
  QPoint pt1;  // origin
  QVector u[DIM];  // make axes
  u[0][0] = ONE_OVER_SQRT_TWO;
  u[0][1] = ONE_OVER_SQRT_TWO;
  u[1][0] = ONE_OVER_SQRT_TWO;
  u[1][1] = -ONE_OVER_SQRT_TWO;
  u[2][2] = 1.;

  QVector e = QVector(1.);
  QOBBox obbox1(pt1, u, e);
  std::vector< QPoint > verts = obbox1.vertices();

  QPoint pt2 = verts[0];
  QPoint pt3(10.*pt2.array());

  // closest point is a vertex
  EXPECT_TRUE(primal::closest_point(pt3, obbox1) == pt2);
}

//------------------------------------------------------------------------------
TEST( primal_closest_point, obb_test_closest_point_face )
{
  static const int DIM = 3;
  typedef double CoordType;
  typedef primal::Point< CoordType, DIM > QPoint;
  typedef primal::Vector< CoordType, DIM > QVector;
  typedef primal::OrientedBoundingBox< CoordType, DIM > QOBBox;

  const double ONE_OVER_SQRT_TWO = 0.7071;
  const double EPS = 0.01;
  QPoint pt1;  // origin
  QVector u[DIM];  // make standard axes
  u[0][0] = ONE_OVER_SQRT_TWO;
  u[0][1] = ONE_OVER_SQRT_TWO;
  u[1][0] = ONE_OVER_SQRT_TWO;
  u[1][1] = -ONE_OVER_SQRT_TWO;
  u[2][2] = 1.;

  QVector e = QVector(1.);
  QOBBox obbox1(pt1, u, e);

  QPoint pt2(u[0].array());
  QPoint pt3(10.*pt2.array());

  QVector found(primal::closest_point(pt3, obbox1));
  QVector expected(pt2);

  // closest point is on a face
  EXPECT_TRUE((found - expected).squared_norm() < EPS);
}

//------------------------------------------------------------------------------
TEST( primal_closest_point, obb_test_closest_point_edge )
{
  static const int DIM = 3;
  typedef double CoordType;
  typedef primal::Point< CoordType, DIM > QPoint;
  typedef primal::Vector< CoordType, DIM > QVector;
  typedef primal::OrientedBoundingBox< CoordType, DIM > QOBBox;

  const double ONE_OVER_SQRT_TWO = 0.7071;
  const double EPS = 0.01;
  QPoint pt1;  // origin
  QVector u[DIM];  // make standard axes
  u[0][0] = ONE_OVER_SQRT_TWO;
  u[0][1] = ONE_OVER_SQRT_TWO;
  u[1][0] = ONE_OVER_SQRT_TWO;
  u[1][1] = -ONE_OVER_SQRT_TWO;
  u[2][2] = 1.;

  QVector e = QVector(1.);
  QOBBox obbox1(pt1, u, e);

  QPoint pt2(u[0].array() + u[1].array());
  QPoint pt3(10.*pt2.array());

  QVector found(primal::closest_point(pt3, obbox1));
  QVector expected(pt2);

  // closest point is on an edge
  EXPECT_TRUE((found - expected).squared_norm() < EPS);
}

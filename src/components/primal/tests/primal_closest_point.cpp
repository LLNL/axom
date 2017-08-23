/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
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

  QPoint pt1;  // origin
  QVector u[DIM];  // make standard axes
  for (int i = 0; i < DIM; i++) {
    u[i] = QVector();
    u[i][i] = 1.;
  }

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

  QPoint pt1;  // origin
  QVector u[DIM];  // make standard axes
  for (int i = 0; i < DIM; i++) {
    u[i] = QVector();
    u[i][i] = 1.;
  }

  QVector e = QVector(1.);
  QPoint pt2 = QPoint(1.);
  QPoint pt3 = QPoint(2.);;

  QOBBox obbox1(pt1, u, e);

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

  QPoint pt1;  // origin
  QVector u[DIM];  // make standard axes
  for (int i = 0; i < DIM; i++) {
    u[i] = QVector();
    u[i][i] = 1.;
  }

  QVector e = QVector(1.);
  QPoint pt2;
  QPoint pt3;

  QOBBox obbox1(pt1, u, e);

  pt2.array()[2] = 2.;
  pt3.array()[2] = 1.;

  // closest point is on a face
  EXPECT_TRUE(primal::closest_point(pt2, obbox1) == pt3);
}

//------------------------------------------------------------------------------
TEST( primal_closest_point, obb_test_closest_point_edge )
{
  static const int DIM = 3;
  typedef double CoordType;
  typedef primal::Point< CoordType, DIM > QPoint;
  typedef primal::Vector< CoordType, DIM > QVector;
  typedef primal::OrientedBoundingBox< CoordType, DIM > QOBBox;

  QPoint pt1;  // origin
  QVector u[DIM];  // make standard axes
  for (int i = 0; i < DIM; i++) {
    u[i] = QVector();
    u[i][i] = 1.;
  }

  QVector e = QVector(1.);
  QPoint pt2(2.);
  QPoint pt3(1.);

  QOBBox obbox1(pt1, u, e);

  pt2.array()[0] = 0.;
  pt3.array()[0] = 0.;

  // closest point is on an edge
  EXPECT_TRUE(primal::closest_point(pt2, obbox1) == pt3);
}

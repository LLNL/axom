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
TEST( primal_closest_point, obb_test_closest_point )
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

  EXPECT_TRUE(primal::closest_point(pt1, obbox1) == pt1);
  EXPECT_TRUE(primal::closest_point(pt3, obbox1) == pt2);

  QPoint pt4;
  QPoint pt5;
  pt4.array()[2] = 2.;
  pt5.array()[2] = 1.;
  EXPECT_TRUE(primal::closest_point(pt4, obbox1) == pt5);

  pt4.array()[1] = 2.;
  pt5.array()[1] = 1.;
  EXPECT_TRUE(primal::closest_point(pt4, obbox1) == pt5);
}

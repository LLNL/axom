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
#include "primal/compute_bounding_box.hpp"

using namespace axom;

TEST( primal_compute_bounding_box, OBB_from_points_test)
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
  QVector e(1.);

  // test basic cube
  QOBBox obbox1(pt1, u, e);

  std::vector< QPoint > v = obbox1.vertices();
  QPoint verts[8];
  for (int i = 0; i < 8; i++) verts[i] = v[i];

  QOBBox obbox2 = primal::OBB_from_points< CoordType, DIM >(verts, 8);

  EXPECT_TRUE(obbox2.contains(obbox1));

  // now test line of points
  QPoint pt2(1.);
  QPoint pt3(2.);
  verts[0] = pt1;
  verts[1] = pt2;
  verts[2] = pt3;

  QOBBox obbox3 = primal::OBB_from_points< CoordType, DIM >(verts, 3);

  EXPECT_TRUE(obbox3.contains(pt1));
  EXPECT_TRUE(obbox3.contains(pt2));
  EXPECT_TRUE(obbox3.contains(pt3));

  EXPECT_TRUE(obbox3.extents()[1] < 1E-8);
  EXPECT_TRUE(obbox3.extents()[2] < 1E-8);
}

// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <limits>
#include <algorithm>

#include "gtest/gtest.h"

#include "axom/primal/geometry/NumericArray.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/OrientedBoundingBox.hpp"
#include "axom/primal/operators/compute_bounding_box.hpp"

using namespace axom;

TEST(primal_compute_bounding_box, compute_oriented_box_test)
{
  static const int DIM = 3;
  typedef double CoordType;
  typedef primal::Point<CoordType, DIM> QPoint;
  typedef primal::Vector<CoordType, DIM> QVector;
  typedef primal::OrientedBoundingBox<CoordType, DIM> QOBBox;

  QPoint pt1;      // origin
  QVector u[DIM];  // make standard axes
  for(int i = 0; i < DIM; i++)
  {
    u[i] = QVector();
    u[i][i] = 1.;
  }
  QVector e(1.);

  // test basic cube
  QOBBox obbox1(pt1, u, e);

  std::vector<QPoint> v = obbox1.vertices();
  QPoint verts[8];
  for(int i = 0; i < 8; i++)
  {
    verts[i] = v[i];
  }

  QOBBox obbox2 = primal::compute_oriented_bounding_box<CoordType, DIM>(verts, 8);

  EXPECT_TRUE(obbox2.contains(obbox1));

  // now test line of points
  QPoint pt2(1.);
  QPoint pt3(2.);
  verts[0] = pt1;
  verts[1] = pt2;
  verts[2] = pt3;

  QOBBox obbox3 = primal::compute_oriented_bounding_box<CoordType, DIM>(verts, 3);

  EXPECT_TRUE(obbox3.contains(pt1));
  EXPECT_TRUE(obbox3.contains(pt2));
  EXPECT_TRUE(obbox3.contains(pt3));

  EXPECT_TRUE(obbox3.getExtents()[1] < 1E-8);
  EXPECT_TRUE(obbox3.getExtents()[2] < 1E-8);
}

TEST(primal_compute_bounding_box, merge_oriented_box_test)
{
  static const int DIM = 3;
  typedef double CoordType;
  typedef primal::Point<CoordType, DIM> QPoint;
  typedef primal::Vector<CoordType, DIM> QVector;
  typedef primal::OrientedBoundingBox<CoordType, DIM> QOBBox;

  const double ONE_OVER_SQRT_TWO = 0.7071;
  QPoint pt1;      // origin
  QVector u[DIM];  // make axes
  u[0][0] = ONE_OVER_SQRT_TWO;
  u[0][1] = ONE_OVER_SQRT_TWO;
  u[1][0] = ONE_OVER_SQRT_TWO;
  u[1][1] = -ONE_OVER_SQRT_TWO;
  u[2][2] = 1.;

  QVector e = QVector(1.);
  QOBBox obbox1(pt1, u, e);

  // test when they're identical
  EXPECT_TRUE((obbox1 == primal::merge_boxes<CoordType, DIM>(obbox1, obbox1)));

  // test when they're completely separate
  QPoint pt2(10.);

  QOBBox obbox2(pt2, u, e);
  QOBBox obbox3 = primal::merge_boxes<CoordType, DIM>(obbox1, obbox2);
  EXPECT_TRUE(obbox3.contains(obbox1));
  EXPECT_TRUE(obbox3.contains(obbox2));

  // test when they only partially intersect
  QPoint pt3(0.5);
  QOBBox obbox4(pt3, u, e);
  QOBBox obbox5 = primal::merge_boxes<CoordType, DIM>(obbox1, obbox4);
  EXPECT_TRUE(obbox5.contains(obbox1));
  EXPECT_TRUE(obbox5.contains(obbox4));
}

TEST(primal_compute_bounding_box, merge_aligned_box_test)
{
  static const int DIM = 3;
  typedef double CoordType;
  typedef primal::Point<CoordType, DIM> QPoint;
  typedef primal::Vector<CoordType, DIM> QVector;
  typedef primal::BoundingBox<CoordType, DIM> QBBox;

  QPoint pt1(0.);
  QPoint pt2(1.);
  QBBox bbox1(pt1);
  QBBox bbox2(pt2);
  QBBox bbox3(pt1, pt2);

  EXPECT_TRUE((bbox3 == primal::merge_boxes<CoordType, DIM>(bbox1, bbox2)));

  QBBox bbox4(bbox3);
  QVector s(10.);
  bbox3.shift(s);

  QBBox bbox5 = primal::merge_boxes<CoordType, DIM>(bbox3, bbox4);

  EXPECT_TRUE(bbox5.contains(bbox3));
  EXPECT_TRUE(bbox5.contains(bbox4));
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);

#ifdef COMPUTE_BOUNDING_BOX_TESTER_SHOULD_SEED
  std::srand(std::time(0));
#else
  std::srand(105);
#endif

  int result = RUN_ALL_TESTS();
  return result;
}

// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/primal/geometry/Plane.hpp"
#include "axom/core/numerics/matvecops.hpp"

#include "gtest/gtest.h"

// C/C++ includes
#include <cmath>

namespace primal = axom::primal;
namespace numerics = axom::numerics;

//------------------------------------------------------------------------------
// INTERNAL HELPER METHODS
//------------------------------------------------------------------------------
namespace
{
using PointType3 = primal::Point<double, 3>;
using VectorType3 = primal::Vector<double, 3>;
using PointType2 = primal::Point<double, 2>;
using VectorType2 = primal::Vector<double, 2>;

template <typename T, int NDIMS>
void ensure_unit_norm(const primal::Vector<T, NDIMS>& v)
{
  const double norm = v.norm();
  EXPECT_DOUBLE_EQ(norm, 1.0);
}

//------------------------------------------------------------------------------
template <int NDIMS>
void check_copy_constructor()
{
  const double MAGIC_NUM = 42;
  primal::Vector<double, NDIMS> n1(MAGIC_NUM);
  double offset = MAGIC_NUM;

  primal::Plane<double, NDIMS> p2(n1, offset);
  primal::Plane<double, NDIMS> p1(p2);

  numerics::normalize(n1.data(), NDIMS);

  for(int i = 0; i < NDIMS; ++i)
  {
    EXPECT_DOUBLE_EQ(p1.getNormal()[i], p2.getNormal()[i]);
    EXPECT_DOUBLE_EQ(p1.getNormal()[i], n1[i]);
  }

  EXPECT_DOUBLE_EQ(p1.getOffset(), p2.getOffset());
  EXPECT_DOUBLE_EQ(p1.getOffset(), MAGIC_NUM);
}

//------------------------------------------------------------------------------
template <int NDIMS>
void check_assignment_operator()
{
  const double MAGIC_NUM = 42;
  primal::Point<double, NDIMS> x(0.0);
  primal::Vector<double, NDIMS> n0(0.0);
  n0[0] = 1.;
  primal::Vector<double, NDIMS> n1(MAGIC_NUM);
  double offset = MAGIC_NUM;

  primal::Plane<double, NDIMS> p1(n0, x);
  primal::Plane<double, NDIMS> p2(n1, offset);

  /* test assignment */
  p1 = p2;

  numerics::normalize(n1.data(), NDIMS);

  for(int i = 0; i < NDIMS; ++i)
  {
    EXPECT_DOUBLE_EQ(p1.getNormal()[i], p2.getNormal()[i]);
    EXPECT_DOUBLE_EQ(p1.getNormal()[i], n1[i]);
  }

  EXPECT_DOUBLE_EQ(p1.getOffset(), p2.getOffset());
  EXPECT_DOUBLE_EQ(p1.getOffset(), MAGIC_NUM);
}
} /* end anonymous namespace */

//------------------------------------------------------------------------------
// UNIT TEST
//------------------------------------------------------------------------------
TEST(primal_plane, construct_from_normal_and_point)
{
  VectorType3 normal {0.0, 0.0, 10.0};
  PointType3 x {0.0, 0.0, 2.0};

  // test 3D
  primal::Plane<double, 3> P(normal, x);
  EXPECT_DOUBLE_EQ(P.getOffset(), 2.0);
  EXPECT_DOUBLE_EQ(P.getNormal()[0], 0.0);
  EXPECT_DOUBLE_EQ(P.getNormal()[1], 0.0);
  EXPECT_DOUBLE_EQ(P.getNormal()[2], 1.0);
  EXPECT_EQ(P.getDimension(), 3);

  // test 2D
  VectorType2 normal2 {1.0, 2.0};
  PointType2 x2 {1.0, 2.0};

  primal::Plane<double, 2> P2(normal2, x2);
  ensure_unit_norm(P2.getNormal());
  EXPECT_DOUBLE_EQ(P2.getOffset(), std::sqrt(5.0));
  EXPECT_EQ(P2.getDimension(), 2);
}

//------------------------------------------------------------------------------
TEST(primal_plane, construct_from_normal_and_offset)
{
  VectorType3 normal {0.0, 0.0, 1.0};
  double offset = 2.0;

  // test 3D
  primal::Plane<double, 3> P(normal, offset);
  EXPECT_DOUBLE_EQ(P.getNormal()[0], 0.0);
  EXPECT_DOUBLE_EQ(P.getNormal()[1], 0.0);
  EXPECT_DOUBLE_EQ(P.getNormal()[2], 1.0);
  EXPECT_DOUBLE_EQ(P.getOffset(), offset);

  // test 2D
  VectorType2 normal2 {1.0, 2.0};
  offset = std::sqrt(5.0);
  primal::Plane<double, 2> P2(normal2, offset);
  ensure_unit_norm(P2.getNormal());
  EXPECT_DOUBLE_EQ(P2.getOffset(), offset);
}

//------------------------------------------------------------------------------
TEST(primal_plane, construct_from_points)
{
  primal::Point<double, 3> x1 {1.0, 1.0, 3.0};
  primal::Point<double, 3> x2 {2.0, 2.0, 3.0};
  primal::Point<double, 3> x3 {1.0, 3.0, 3.0};

  // test 3D
  primal::Plane<double, 3> P = primal::make_plane(x1, x2, x3);
  ensure_unit_norm(P.getNormal());
  EXPECT_DOUBLE_EQ(P.getOffset(), 3.0);

  // test 2D
  primal::Point<double, 2> a {2.0, -1.0};
  primal::Point<double, 2> b {2.0, 2.0};
  primal::Plane<double, 2> P2 = primal::make_plane(a, b);
  ensure_unit_norm(P2.getNormal());
  EXPECT_DOUBLE_EQ(P2.getOffset(), -2.0);
}

//------------------------------------------------------------------------------
TEST(primal_plane, copy_constructor)
{
  check_copy_constructor<2>();
  check_copy_constructor<3>();
}

//------------------------------------------------------------------------------
TEST(primal_plane, assignment_operator)
{
  check_assignment_operator<2>();
  check_assignment_operator<3>();
}

//------------------------------------------------------------------------------
TEST(primal_plane, signed_distance_and_orientation)
{
  primal::Point<double, 3> x1 {1.0, 1.0, 3.0};
  primal::Point<double, 3> x2 {2.0, 2.0, 3.0};
  primal::Point<double, 3> x3 {1.0, 3.0, 3.0};

  double signed_distance = 0.0;  // stores computed signed distance.
  PointType3 q {0.0, 0.0, 0.0};  // test query point

  // STEP 0: test 3D
  primal::Plane<double, 3> P = primal::make_plane(x1, x2, x3);

  // (a) test point below plane
  signed_distance = P.signedDistance(q);
  EXPECT_DOUBLE_EQ(signed_distance, -3.0);
  EXPECT_EQ(P.getOrientation(q), primal::ON_NEGATIVE_SIDE);

  // (b) test point above plane
  q[2] = 6.0;
  signed_distance = P.signedDistance(q);
  EXPECT_DOUBLE_EQ(signed_distance, 3.0);
  EXPECT_EQ(P.getOrientation(q), primal::ON_POSITIVE_SIDE);

  // (c) test point on plane
  q[2] = 3.0;
  signed_distance = P.signedDistance(q);
  EXPECT_DOUBLE_EQ(signed_distance, 0.0);
  EXPECT_EQ(P.getOrientation(q), primal::ON_BOUNDARY);

  // STEP 1: test 2D
  primal::Point<double, 2> a {2.0, -1.0};
  primal::Point<double, 2> b {2.0, 2.0};
  primal::Plane<double, 2> P2 = primal::make_plane(a, b);
  PointType2 q2 {0.0, 0.0};  // test query point

  // (a) test point above plane
  signed_distance = P2.signedDistance(q2);
  EXPECT_DOUBLE_EQ(signed_distance, 2.0);
  EXPECT_EQ(P2.getOrientation(q2), primal::ON_POSITIVE_SIDE);

  // (b) test point below plane
  q2[0] = 4.0;
  signed_distance = P2.signedDistance(q2);
  EXPECT_DOUBLE_EQ(signed_distance, -2.0);
  EXPECT_EQ(P2.getOrientation(q2), primal::ON_NEGATIVE_SIDE);

  // (c) test point on plane
  q2[0] = 2.0;
  signed_distance = P2.signedDistance(q2);
  EXPECT_DOUBLE_EQ(signed_distance, 0.0);
  EXPECT_EQ(P2.getOrientation(q2), primal::ON_BOUNDARY);
}

//------------------------------------------------------------------------------
TEST(primal_plane, project_point)
{
  primal::Point<double, 3> x1 {1.0, 1.0, 3.0};
  primal::Point<double, 3> x2 {2.0, 2.0, 3.0};
  primal::Point<double, 3> x3 {1.0, 3.0, 3.0};

  primal::Point<double, 3> q {0.0, 0.0, 0.0};
  primal::Point<double, 3> qproj {0.0, 0.0, 0.0};

  // STEP 0: test 3D
  primal::Plane<double, 3> P = primal::make_plane(x1, x2, x3);

  // (a) test project point below plane
  qproj = P.projectPoint(q);
  EXPECT_EQ(P.getOrientation(qproj), primal::ON_BOUNDARY);
  EXPECT_DOUBLE_EQ(qproj[0], 0.0);
  EXPECT_DOUBLE_EQ(qproj[1], 0.0);
  EXPECT_DOUBLE_EQ(qproj[2], 3.0);

  // (b) test project point above plane
  q[2] = 6.0;
  qproj[0] = qproj[1] = qproj[2] = 0.0;
  qproj = P.projectPoint(q);
  EXPECT_EQ(P.getOrientation(qproj), primal::ON_BOUNDARY);
  EXPECT_DOUBLE_EQ(qproj[0], 0.0);
  EXPECT_DOUBLE_EQ(qproj[1], 0.0);
  EXPECT_DOUBLE_EQ(qproj[2], 3.0);

  // (c) test project point (already) on plane
  q[2] = 3.0;
  qproj[0] = qproj[1] = qproj[2] = 0.0;
  qproj = P.projectPoint(q);
  EXPECT_EQ(P.getOrientation(qproj), primal::ON_BOUNDARY);
  EXPECT_DOUBLE_EQ(qproj[0], q[0]);
  EXPECT_DOUBLE_EQ(qproj[1], q[1]);
  EXPECT_DOUBLE_EQ(qproj[2], q[2]);

  // STEP 1: test 2D
  primal::Point<double, 2> a {2.0, -1.0};
  primal::Point<double, 2> b {2.0, 2.0};
  primal::Plane<double, 2> P2 = primal::make_plane(a, b);
  primal::Point<double, 2> q2 {0.0, 0.0};
  primal::Point<double, 2> qproj2 {0.0, 0.0};

  // (a) test project point below plane
  q2[0] = 4.0;
  qproj2[0] = qproj2[1] = 0.0;
  qproj2 = P2.projectPoint(q2);
  EXPECT_EQ(P2.getOrientation(qproj2), primal::ON_BOUNDARY);
  EXPECT_DOUBLE_EQ(qproj2[0], 2.0);
  EXPECT_DOUBLE_EQ(qproj2[1], 0.0);

  // (b) test project point above plane
  q2[0] = 0.0;
  qproj2[0] = qproj2[1] = 0.0;
  qproj2 = P2.projectPoint(q2);
  EXPECT_EQ(P2.getOrientation(qproj2), primal::ON_BOUNDARY);
  EXPECT_DOUBLE_EQ(qproj2[0], 2.0);
  EXPECT_DOUBLE_EQ(qproj2[1], 0.0);

  // (c) test project point (already) on plane
  q2[0] = 2.0;
  qproj2[0] = qproj2[1] = 0.0;
  qproj2 = P2.projectPoint(q2);
  EXPECT_EQ(P2.getOrientation(qproj2), primal::ON_BOUNDARY);
  EXPECT_DOUBLE_EQ(qproj2[0], q2[0]);
  EXPECT_DOUBLE_EQ(qproj2[1], q2[1]);
}

//------------------------------------------------------------------------------
TEST(primal_plane, flip)
{
  primal::Point<double, 3> x1 {1.0, 1.0, 3.0};
  primal::Point<double, 3> x2 {2.0, 2.0, 3.0};
  primal::Point<double, 3> x3 {1.0, 3.0, 3.0};

  PointType3 q {0.0, 0.0, 0.0};

  // STEP 0: test 3D
  primal::Plane<double, 3> P = primal::make_plane(x1, x2, x3);
  EXPECT_EQ(P.getOrientation(q), primal::ON_NEGATIVE_SIDE);
  P.flip();
  EXPECT_EQ(P.getOrientation(q), primal::ON_POSITIVE_SIDE);

  // STEP 1: test 2D
  primal::Point<double, 2> a {2.0, -1.0};
  primal::Point<double, 2> b {2.0, 2.0};
  primal::Plane<double, 2> P2 = primal::make_plane(a, b);
  PointType2 q2 {0.0, 0.0};
  EXPECT_EQ(P2.getOrientation(q2), primal::ON_POSITIVE_SIDE);
  P2.flip();
  EXPECT_EQ(P2.getOrientation(q2), primal::ON_NEGATIVE_SIDE);
}

//------------------------------------------------------------------------------
#include "axom/slic/core/SimpleLogger.hpp"
using axom::slic::SimpleLogger;

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  SimpleLogger logger;  // create & initialize test logger,

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}

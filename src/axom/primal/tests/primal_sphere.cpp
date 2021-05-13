// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/primal/geometry/Sphere.hpp"

#include "gtest/gtest.h"

namespace primal = axom::primal;

//------------------------------------------------------------------------------
// INTERNAL HELPER METHODS
//------------------------------------------------------------------------------
namespace
{
void check_array(const double* a, const double* b, int size)
{
  for(int i = 0; i < size; ++i)
  {
    EXPECT_DOUBLE_EQ(a[i], b[i]);
  }
}

//------------------------------------------------------------------------------
template <int NDIMS>
void check_constructor()
{
  const double DEFAULT_RADIUS = 1.0;
  const double DEFAULT_CENTER[] = {0.0, 0.0, 0.0};
  const double PRESCRIBED_RADIUS = 5.0;
  const double PRESCRIBED_CENTER[] = {1.0, 1.0, 1.0};

  primal::Sphere<double, NDIMS> S0;
  EXPECT_DOUBLE_EQ(S0.getRadius(), DEFAULT_RADIUS);
  check_array(S0.getCenter(), DEFAULT_CENTER, NDIMS);

  primal::Sphere<double, NDIMS> S1(PRESCRIBED_RADIUS);
  EXPECT_DOUBLE_EQ(S1.getRadius(), PRESCRIBED_RADIUS);
  check_array(S1.getCenter(), DEFAULT_CENTER, NDIMS);

  primal::Sphere<double, NDIMS> S2(PRESCRIBED_CENTER);
  EXPECT_DOUBLE_EQ(S2.getRadius(), DEFAULT_RADIUS);
  check_array(S2.getCenter(), PRESCRIBED_CENTER, NDIMS);

  primal::Sphere<double, NDIMS> S4(PRESCRIBED_CENTER, PRESCRIBED_RADIUS);
  EXPECT_DOUBLE_EQ(S4.getRadius(), PRESCRIBED_RADIUS);
  check_array(S4.getCenter(), PRESCRIBED_CENTER, NDIMS);
}

//------------------------------------------------------------------------------
template <int NDIMS>
void check_signed_distance_and_orientation()
{
  double x[3] = {0.0, 0.0, 0.0};
  double signed_distance = 0.0;

  // STEP 0: test 2D sphere (i.e., circle) with radius 1.0 centered @ (0,0)
  primal::Sphere<double, NDIMS> sphere;
  const double radius = sphere.getRadius();
  const double* center = sphere.getCenter();
  EXPECT_TRUE(center != nullptr);

  // test sphere center
  signed_distance = sphere.computeSignedDistance(center);
  EXPECT_DOUBLE_EQ(signed_distance, (-1.0) * radius);
  EXPECT_EQ(sphere.getOrientation(center), primal::ON_NEGATIVE_SIDE);

  for(int j = 0; j < NDIMS; ++j)
  {
    // initialize test point, x
    memcpy(x, center, NDIMS * sizeof(double));

    // SHIFT RIGHT
    // shift right from center, but still within the sphere
    x[j] = center[j] + 0.5 * radius;
    signed_distance = sphere.computeSignedDistance(x);
    EXPECT_DOUBLE_EQ(signed_distance, (-0.5) * radius);
    EXPECT_EQ(sphere.getOrientation(x), primal::ON_NEGATIVE_SIDE);

    // shift right all the way on the sphere boundary
    x[j] = center[j] + radius;
    signed_distance = sphere.computeSignedDistance(x);
    EXPECT_DOUBLE_EQ(signed_distance, 0.0);
    EXPECT_EQ(sphere.getOrientation(x), primal::ON_BOUNDARY);

    // shift right outside sphere
    x[j] = center[j] + 2 * radius;
    signed_distance = sphere.computeSignedDistance(x);
    EXPECT_DOUBLE_EQ(signed_distance, radius);
    EXPECT_EQ(sphere.getOrientation(x), primal::ON_POSITIVE_SIDE);

    // SHIFT LEFT
    // shift left from center, but still within the sphere
    x[j] = center[j] - 0.5 * radius;
    signed_distance = sphere.computeSignedDistance(x);
    EXPECT_DOUBLE_EQ(signed_distance, (-0.5) * radius);
    EXPECT_EQ(sphere.getOrientation(x), primal::ON_NEGATIVE_SIDE);

    // shift left all the way on the sphere boundary
    x[j] = center[j] - radius;
    signed_distance = sphere.computeSignedDistance(x);
    EXPECT_DOUBLE_EQ(signed_distance, 0.0);
    EXPECT_EQ(sphere.getOrientation(x), primal::ON_BOUNDARY);

    // shift left outside sphere
    x[j] = center[j] - 2 * radius;
    signed_distance = sphere.computeSignedDistance(x);
    EXPECT_DOUBLE_EQ(signed_distance, radius);
    EXPECT_EQ(sphere.getOrientation(x), primal::ON_POSITIVE_SIDE);

  }  // END for all dimensions
}

//------------------------------------------------------------------------------
template <int NDIMS>
void check_sphere_intersection()
{
  double center[3] = {0.0, 0.0, 0.0};

  // STEP 0: test fully overlapping
  primal::Sphere<double, NDIMS> S0;
  EXPECT_TRUE(S0.intersectsWith(S0));

  // STEP 1: test inter-penetrating
  center[0] = 0.5;
  primal::Sphere<double, NDIMS> S1(center);
  EXPECT_TRUE(S0.intersectsWith(S1));

  // STEP 2: test abutting
  center[0] = 2.0;
  primal::Sphere<double, NDIMS> S2(center);
  EXPECT_TRUE(S0.intersectsWith(S2));

  // STEP 3: test non-intersecting
  center[0] = 4.0;
  primal::Sphere<double, NDIMS> S3(center);
  EXPECT_FALSE(S0.intersectsWith(S3));
}

//------------------------------------------------------------------------------
template <int NDIMS>
void check_copy_constructor()
{
  const double MAGIC_NUM = 42;
  double center[3] = {MAGIC_NUM, MAGIC_NUM, MAGIC_NUM};
  double radius = MAGIC_NUM;

  primal::Sphere<double, NDIMS> s1(center, radius);
  primal::Sphere<double, NDIMS> s2(s1);

  const double* c1 = s1.getCenter();
  const double* c2 = s2.getCenter();
  for(int i = 0; i < NDIMS; ++i)
  {
    EXPECT_DOUBLE_EQ(c1[i], c2[i]);
    EXPECT_DOUBLE_EQ(c2[i], MAGIC_NUM);
  }

  EXPECT_DOUBLE_EQ(s1.getRadius(), s2.getRadius());
  EXPECT_DOUBLE_EQ(s2.getRadius(), MAGIC_NUM);
}

//------------------------------------------------------------------------------
template <int NDIMS>
void check_assignment_operator()
{
  const double MAGIC_NUM = 42;
  double center[3] = {MAGIC_NUM, MAGIC_NUM, MAGIC_NUM};
  double radius = MAGIC_NUM;

  primal::Sphere<double, NDIMS> s1(center, radius);
  primal::Sphere<double, NDIMS> s2;

  /* test assignment */
  s2 = s1;

  const double* c1 = s1.getCenter();
  const double* c2 = s2.getCenter();
  for(int i = 0; i < NDIMS; ++i)
  {
    EXPECT_DOUBLE_EQ(c1[i], c2[i]);
    EXPECT_DOUBLE_EQ(c2[i], MAGIC_NUM);
  }

  EXPECT_DOUBLE_EQ(s1.getRadius(), s2.getRadius());
  EXPECT_DOUBLE_EQ(s2.getRadius(), MAGIC_NUM);
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST(primal_sphere, constructor)
{
  check_constructor<2>();
  check_constructor<3>();
}

//------------------------------------------------------------------------------
TEST(primal_sphere, copy_constructor)
{
  check_copy_constructor<2>();
  check_copy_constructor<3>();
}

//------------------------------------------------------------------------------
TEST(primal_sphere, assignment_operator)
{
  check_assignment_operator<2>();
  check_assignment_operator<3>();
}

//------------------------------------------------------------------------------
TEST(primal_sphere, signed_distance_and_orientation)
{
  check_signed_distance_and_orientation<2>();
  check_signed_distance_and_orientation<3>();
}

//------------------------------------------------------------------------------
TEST(primal_sphere, sphere_sphere_intersection)
{
  check_sphere_intersection<2>();
  check_sphere_intersection<3>();
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

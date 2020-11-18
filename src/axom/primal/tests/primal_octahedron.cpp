// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Octahedron.hpp"

#include "fmt/fmt.hpp"
#include "axom/slic/interface/slic.hpp"

#include <cmath>

using namespace axom;

/// Test fixture for testing primal::Octahedron
class OctahedronTest : public ::testing::Test
{
public:
  static const int DIM = 3;

  typedef double CoordType;
  typedef primal::Point<CoordType, DIM> QPoint;
  typedef primal::Octahedron<CoordType, DIM> QOct;

protected:
  virtual void SetUp()
  {
    EPS = 1e-12;

    // Define coordinates for first octahedron
    qData0[0] = QPoint::make_point(1, 0, 0);
    qData0[1] = QPoint::make_point(1, 1, 0);
    qData0[2] = QPoint::make_point(0, 1, 0);
    qData0[3] = QPoint::make_point(0, 1, 1);
    qData0[4] = QPoint::make_point(0, 0, 1);
    qData0[5] = QPoint::make_point(1, 0, 1);

    // Define coordinates for second octahedron
    qData1[0] = QPoint::make_point(0, 0, -1);
    qData1[1] = QPoint::make_point(1, 0, 0);
    qData1[2] = QPoint::make_point(0, 1, 0);
    qData1[3] = QPoint::make_point(0, 0, 1);
    qData1[4] = QPoint::make_point(-1, 0, 0);
    qData1[5] = QPoint::make_point(0, -1, 0);
  }

  QPoint qData0[6];
  QPoint qData1[6];
  double EPS;
};

//------------------------------------------------------------------------------
TEST_F(OctahedronTest, defaultConstructor)
{
  typedef OctahedronTest::QPoint QPoint;
  typedef OctahedronTest::QOct QOct;

  const QOct oct;

  // Test ostream operator
  SLIC_INFO("Empty octahedron coordinates: " << oct);

  // Check indirection operator
  EXPECT_EQ(QPoint::zero(), oct[0]);
  EXPECT_EQ(QPoint::zero(), oct[1]);
  EXPECT_EQ(QPoint::zero(), oct[2]);
  EXPECT_EQ(QPoint::zero(), oct[3]);
  EXPECT_EQ(QPoint::zero(), oct[4]);
  EXPECT_EQ(QPoint::zero(), oct[5]);
}

TEST_F(OctahedronTest, constructFromPoints)
{
  typedef OctahedronTest::QPoint QPoint;
  typedef OctahedronTest::QOct QOct;

  // Access the test data
  const QPoint* pt = this->qData0;

  QOct oct(pt[0], pt[1], pt[2], pt[3], pt[4], pt[5]);

  // Test ostream operator
  SLIC_INFO("Octahedron coordinates: " << oct);

  // Check indirection operator
  EXPECT_EQ(pt[0], oct[0]);
  EXPECT_EQ(pt[1], oct[1]);
  EXPECT_EQ(pt[2], oct[2]);
  EXPECT_EQ(pt[3], oct[3]);
  EXPECT_EQ(pt[4], oct[4]);
  EXPECT_EQ(pt[5], oct[5]);
}

TEST_F(OctahedronTest, equals)
{
  typedef OctahedronTest::QPoint QPoint;
  typedef OctahedronTest::QOct QOct;

  // Access the test data
  const QPoint* pt0 = this->qData0;
  const QPoint* pt1 = this->qData1;

  QOct oct0(pt0[0], pt0[1], pt0[2], pt0[3], pt0[4], pt0[5]);
  // Permute the points (but keep the same shape)
  QOct oct0a(pt0[0], pt0[1], pt0[2], pt0[3], pt0[4], pt0[5]);
  // Different shape
  QOct oct1(pt1[0], pt1[1], pt1[2], pt1[3], pt1[4], pt1[5]);

  // Oct should be equal to itself.
  EXPECT_TRUE(oct0.equals(oct0));
  // Oct should be equal to its points permuted.
  EXPECT_TRUE(oct0.equals(oct0a));
  // Oct should differ from an oct with different points.
  EXPECT_FALSE(oct0.equals(oct1));
  // Another oct should be equal to itself.
  EXPECT_TRUE(oct1.equals(oct1));
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
#include "axom/slic/core/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,
  axom::slic::setLoggingMsgLevel(axom::slic::message::Info);

  int result = RUN_ALL_TESTS();
  return result;
}

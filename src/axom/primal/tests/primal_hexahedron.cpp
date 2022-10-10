// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/slic.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Hexahedron.hpp"

#include <cmath>

namespace primal = axom::primal;

/// Test fixture for testing primal::Hexahedron
class HexahedronTest : public ::testing::Test
{
public:
  static const int DIM = 3;

  typedef double CoordType;
  typedef primal::Point<CoordType, DIM> QPoint;
  typedef primal::Hexahedron<CoordType, DIM> QHex;

protected:
  virtual void SetUp()
  {
    EPS = 1e-12;

    // Define coordinates for first hexahedron
    qData0[0] = QPoint::make_point(0, 0, 0);
    qData0[1] = QPoint::make_point(1, 0, 0);
    qData0[2] = QPoint::make_point(1, 0, 1);
    qData0[3] = QPoint::make_point(0, 0, 1);
    qData0[4] = QPoint::make_point(0, 1, 0);
    qData0[5] = QPoint::make_point(1, 1, 0);
    qData0[6] = QPoint::make_point(1, 1, 1);
    qData0[7] = QPoint::make_point(0, 1, 1);

    // Define coordinates for second hexahedron
    qData1[0] = QPoint::make_point(-1, 0, 0);
    qData1[1] = QPoint::make_point(0, 0, 0);
    qData1[2] = QPoint::make_point(0, 0, 1);
    qData1[3] = QPoint::make_point(-1, 0, 1);
    qData1[4] = QPoint::make_point(-1, 1, 0);
    qData1[5] = QPoint::make_point(0, 1, 0);
    qData1[6] = QPoint::make_point(0, 1, 1);
    qData1[7] = QPoint::make_point(-1, 1, 1);

    // Square frustum
    // (base side length 5, top side length 2, height 1)
    qData2[0] = QPoint::make_point(0, 0, 0);
    qData2[1] = QPoint::make_point(5, 0, 0);
    qData2[2] = QPoint::make_point(5, 0, 5);
    qData2[3] = QPoint::make_point(0, 0, 5);
    qData2[4] = QPoint::make_point(1.5, 1, 1.5);
    qData2[5] = QPoint::make_point(3.5, 1, 1.5);
    qData2[6] = QPoint::make_point(3.5, 1, 3.5);
    qData2[7] = QPoint::make_point(1.5, 1, 3.5);
  }

  QPoint qData0[8];
  QPoint qData1[8];
  QPoint qData2[8];
  double EPS;
};

//------------------------------------------------------------------------------
TEST_F(HexahedronTest, defaultConstructor)
{
  typedef HexahedronTest::QPoint QPoint;
  typedef HexahedronTest::QHex QHex;

  const QHex hex;

  // Test ostream operator
  SLIC_INFO("Empty Hexahedron coordinates: " << hex);

  // Check indirection operator
  EXPECT_EQ(QPoint::zero(), hex[0]);
  EXPECT_EQ(QPoint::zero(), hex[1]);
  EXPECT_EQ(QPoint::zero(), hex[2]);
  EXPECT_EQ(QPoint::zero(), hex[3]);
  EXPECT_EQ(QPoint::zero(), hex[4]);
  EXPECT_EQ(QPoint::zero(), hex[5]);
  EXPECT_EQ(QPoint::zero(), hex[6]);
  EXPECT_EQ(QPoint::zero(), hex[7]);
}

TEST_F(HexahedronTest, constructFromPoints)
{
  typedef HexahedronTest::QPoint QPoint;
  typedef HexahedronTest::QHex QHex;

  // Access the test data
  const QPoint* pt = this->qData0;

  QHex hex(pt[0], pt[1], pt[2], pt[3], pt[4], pt[5], pt[6], pt[7]);

  // Test ostream operator
  SLIC_INFO("Hexahedron coordinates: " << hex);

  // Check indirection operator
  EXPECT_EQ(pt[0], hex[0]);
  EXPECT_EQ(pt[1], hex[1]);
  EXPECT_EQ(pt[2], hex[2]);
  EXPECT_EQ(pt[3], hex[3]);
  EXPECT_EQ(pt[4], hex[4]);
  EXPECT_EQ(pt[5], hex[5]);
  EXPECT_EQ(pt[6], hex[6]);
  EXPECT_EQ(pt[7], hex[7]);
}

TEST_F(HexahedronTest, volume)
{
  typedef HexahedronTest::QPoint QPoint;
  typedef HexahedronTest::QHex QHex;

  // Access the test data
  const QPoint* pt0 = this->qData0;
  const QPoint* pt1 = this->qData1;
  const QPoint* pt2 = this->qData2;

  // Initialize hexahedrons
  QHex hex0(pt0[0], pt0[1], pt0[2], pt0[3], pt0[4], pt0[5], pt0[6], pt0[7]);
  QHex hex1(pt1[0], pt1[1], pt1[2], pt1[3], pt1[4], pt1[5], pt1[6], pt1[7]);
  QHex hex2(pt2[0], pt2[1], pt2[2], pt2[3], pt2[4], pt2[5], pt2[6], pt2[7]);

  // Check volume
  EXPECT_DOUBLE_EQ(hex0.volume(), 1);
  EXPECT_DOUBLE_EQ(hex1.volume(), 1);
  EXPECT_DOUBLE_EQ(hex2.volume(), 13);
}

TEST_F(HexahedronTest, equals)
{
  typedef HexahedronTest::QPoint QPoint;
  typedef HexahedronTest::QHex QHex;

  // Access the test data
  const QPoint* pt0 = this->qData0;
  const QPoint* pt1 = this->qData1;

  QHex hex0(pt0[0], pt0[1], pt0[2], pt0[3], pt0[4], pt0[5], pt0[6], pt0[7]);
  // Permute the points (but keep the same shape)
  QHex hex0a(pt0[0], pt0[1], pt0[2], pt0[3], pt0[4], pt0[5], pt0[6], pt0[7]);
  // Different shape
  QHex hex1(pt1[0], pt1[1], pt1[2], pt1[3], pt1[4], pt1[5], pt1[6], pt1[7]);

  // Hex should be equal to itself.
  EXPECT_TRUE(hex0.equals(hex0));
  // Hex should be equal to its points permuted.
  EXPECT_TRUE(hex0.equals(hex0a));
  // Hex should differ from an hex with different points.
  EXPECT_FALSE(hex0.equals(hex1));
  // Another hex should be equal to itself.
  EXPECT_TRUE(hex1.equals(hex1));
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger(axom::slic::message::Info);

  int result = RUN_ALL_TESTS();
  return result;
}
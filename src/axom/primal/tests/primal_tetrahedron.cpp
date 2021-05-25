// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/slic.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Tetrahedron.hpp"

#include "fmt/fmt.hpp"

#include <cmath>

using namespace axom;

/// Test fixture for testing primal::Tetrahedron
class TetrahedronTest : public ::testing::Test
{
public:
  static const int DIM = 3;

  typedef double CoordType;
  typedef primal::Point<CoordType, DIM> QPoint;
  typedef primal::Tetrahedron<CoordType, DIM> QTet;

protected:
  virtual void SetUp()
  {
    EPS = 1e-12;

    // Define coordinates for first tetrahedron
    qData0[0] = QPoint::make_point(0, 0, 0);
    qData0[1] = QPoint::make_point(1, 0, 0);
    qData0[2] = QPoint::make_point(1, 1, 0);
    qData0[3] = QPoint::make_point(1, 1, 1);

    // Define coordinates for second tetrahedron
    qData1[0] = QPoint::make_point(1, 0, 0);
    qData1[1] = QPoint::make_point(0, 1, 0);
    qData1[2] = QPoint::make_point(0, 0, 1);
    qData1[3] = QPoint::make_point(0, 0, 0);
  }

  QPoint qData0[4];
  QPoint qData1[4];
  double EPS;
};

//------------------------------------------------------------------------------
TEST_F(TetrahedronTest, defaultConstructor)
{
  typedef TetrahedronTest::QPoint QPoint;
  typedef TetrahedronTest::QTet QTet;

  const QTet tet;

  // Test ostream operator
  SLIC_INFO("Empty tetrahedron coordinates: " << tet);

  // Check indirection operator
  EXPECT_EQ(QPoint::zero(), tet[0]);
  EXPECT_EQ(QPoint::zero(), tet[1]);
  EXPECT_EQ(QPoint::zero(), tet[2]);
  EXPECT_EQ(QPoint::zero(), tet[3]);
}

TEST_F(TetrahedronTest, constructFromPoints)
{
  typedef TetrahedronTest::QPoint QPoint;
  typedef TetrahedronTest::QTet QTet;

  // Access the test data
  const QPoint* pt = this->qData0;

  QTet tet(pt[0], pt[1], pt[2], pt[3]);

  // Test ostream operator
  SLIC_INFO("Tetrahedron coordinates: " << tet);

  // Check indirection operator
  EXPECT_EQ(pt[0], tet[0]);
  EXPECT_EQ(pt[1], tet[1]);
  EXPECT_EQ(pt[2], tet[2]);
  EXPECT_EQ(pt[3], tet[3]);
}

TEST_F(TetrahedronTest, volume)
{
  typedef TetrahedronTest::QPoint QPoint;
  typedef TetrahedronTest::QTet QTet;

  // Access the test data
  const QPoint* pt = this->qData0;
  QTet tet(pt[0], pt[1], pt[2], pt[3]);

  double expVolume = 1. / 6.;
  EXPECT_EQ(expVolume, tet.signedVolume());
  EXPECT_EQ(tet.signedVolume(), tet.volume());
}

TEST_F(TetrahedronTest, degenerate)
{
  typedef TetrahedronTest::QPoint QPoint;
  typedef TetrahedronTest::QTet QTet;

  // Access the test data
  const QPoint* pt = this->qData0;
  QTet tet(pt[0], pt[1], pt[2], pt[3]);

  EXPECT_FALSE(tet.degenerate());

  // Make the tet degenerate by identifying two vertices
  tet[0] = tet[1];
  EXPECT_TRUE(tet.degenerate());
}

TEST_F(TetrahedronTest, barycentric)
{
  typedef TetrahedronTest::CoordType CoordType;
  typedef TetrahedronTest::QPoint QPoint;
  typedef TetrahedronTest::QTet QTet;

  typedef primal::Point<CoordType, 4> RPoint;
  typedef std::vector<std::pair<QPoint, RPoint>> TestVec;

  const QPoint* pt = this->qData1;
  QTet tet(pt[0], pt[1], pt[2], pt[3]);

  TestVec testData;

  // Test the four vertices
  const CoordType coord_list0[] = {1., 0., 0., 0.};
  const CoordType coord_list1[] = {0., 1., 0., 0.};
  const CoordType coord_list2[] = {0., 0., 1., 0.};
  const CoordType coord_list3[] = {0., 0., 0., 1.};
  testData.push_back(std::make_pair(pt[0], RPoint(coord_list0, 4)));
  testData.push_back(std::make_pair(pt[1], RPoint(coord_list1, 4)));
  testData.push_back(std::make_pair(pt[2], RPoint(coord_list2, 4)));
  testData.push_back(std::make_pair(pt[3], RPoint(coord_list3, 4)));

  // Test some of the edge midpoints
  const CoordType coord_list4[] = {0.5, 0.5, 0., 0.};
  const CoordType coord_list5[] = {0., 0.5, 0.5, 0.};
  const CoordType coord_list6[] = {0., 0., 0.5, 0.5};
  testData.push_back(std::make_pair(QPoint(0.5 * (pt[0].array() + pt[1].array())),
                                    RPoint(coord_list4, 4)));
  testData.push_back(std::make_pair(QPoint(0.5 * (pt[1].array() + pt[2].array())),
                                    RPoint(coord_list5, 4)));
  testData.push_back(std::make_pair(QPoint(0.5 * (pt[2].array() + pt[3].array())),
                                    RPoint(coord_list6, 4)));

  // Test the centroid
  const CoordType coord_list7[] = {1. / 4., 1. / 4., 1. / 4., 1. / 4.};
  testData.push_back(std::make_pair(
    QPoint(1. / 4. *
           (pt[0].array() + pt[1].array() + pt[2].array() + pt[3].array())),
    RPoint(coord_list7, 4)));

  // Test a point outside the tetrahedron
  const CoordType coord_list8[] = {-0.4, 1.2, 0.2, 0.};
  testData.push_back(std::make_pair(
    QPoint(-0.4 * pt[0].array() + 1.2 * pt[1].array() + 0.2 * pt[2].array()),
    RPoint(coord_list8, 4)));

  // Now run the actual tests
  for(TestVec::const_iterator it = testData.begin(); it != testData.end(); ++it)
  {
    const QPoint& query = it->first;
    const RPoint& expBary = it->second;
    RPoint bary = tet.physToBarycentric(query);

    SLIC_DEBUG(fmt::format(
      "Computed barycentric coordinates for tetrahedron {} and point {} are {}",
      tet,
      query,
      bary));

    EXPECT_NEAR(bary[0], expBary[0], this->EPS);
    EXPECT_NEAR(bary[1], expBary[1], this->EPS);
    EXPECT_NEAR(bary[2], expBary[2], this->EPS);
    EXPECT_NEAR(bary[3], expBary[3], this->EPS);
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
using axom::slic::SimpleLogger;

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);

  SimpleLogger logger;  // create & initialize test logger,
  axom::slic::setLoggingMsgLevel(axom::slic::message::Info);

  int result = RUN_ALL_TESTS();
  return result;
}

/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */



#include "gtest/gtest.h"


#include <cmath>

#include "fmt/fmt.hpp"
#include "slic/slic.hpp"

#include "quest/Point.hpp"
#include "quest/Triangle.hpp"


//------------------------------------------------------------------------------
TEST( quest_triangle, triangle_area_2D)
{
    static const int DIM = 2;
    static const double EPS = 1e-12;
    typedef double CoordType;
    typedef quest::Point<CoordType, DIM> QPoint;
    typedef quest::Triangle<CoordType, DIM> QTri;

    QPoint pt[3] = {
                   QPoint::make_point(0,0),
                   QPoint::make_point(0,1),
                   QPoint::make_point(1,0),
            };

    QTri tri(pt[0],pt[1],pt[2]);
    EXPECT_NEAR(tri.area(), 0.5, EPS );

    tri = QTri(pt[1],pt[2],pt[0]);
    EXPECT_NEAR(tri.area(), 0.5, EPS );

    tri = QTri(pt[2],pt[1],pt[0]);
    EXPECT_NEAR(tri.area(), 0.5, EPS );

    tri = QTri(pt[0],pt[2],pt[1]);
    EXPECT_NEAR(tri.area(), 0.5, EPS );
}

TEST( quest_triangle, triangle_area_3D)
{
    static const int DIM = 3;
    static const double EPS = 1e-12;
    typedef double CoordType;
    typedef quest::Point<CoordType, DIM> QPoint;
    typedef quest::Triangle<CoordType, DIM> QTri;

    QPoint pt[4] = {
                   QPoint::make_point(0,0,0),
                   QPoint::make_point(1,0,0),
                   QPoint::make_point(0,1,0),
                   QPoint::make_point(0,0,1),
            };

    QTri tri(pt[0],pt[1],pt[2]);
    EXPECT_NEAR(tri.area(), 0.5, EPS );

    tri = QTri(pt[0],pt[2],pt[3]);
    EXPECT_NEAR(tri.area(), 0.5, EPS );

    tri = QTri(pt[0],pt[1],pt[3]);
    EXPECT_NEAR(tri.area(), 0.5, EPS );

    tri = QTri(pt[1],pt[2],pt[3]);
    EXPECT_NEAR(tri.area(), std::sqrt(3)/2., EPS );
}


//------------------------------------------------------------------------------
TEST( quest_triangle, triangle_barycentric)
{
  static const int DIM = 3;
  static const double EPS = 1e-12;
  typedef double CoordType;
  typedef quest::Point<CoordType, DIM> QPoint;
  typedef quest::Triangle<CoordType, DIM> QTri;

  QPoint pt[3] = {
                 QPoint::make_point(1,0,0),
                 QPoint::make_point(0,1,0),
                 QPoint::make_point(0,0,1),
          };

  QTri tri(pt[0],pt[1],pt[2]);

  typedef std::vector<std::pair<QPoint,QPoint> > TestVec;
  TestVec testData;

  // Test the three vertices
  testData.push_back( std::make_pair( pt[0], QPoint::make_point(1.,0.,0.)));
  testData.push_back( std::make_pair( pt[1], QPoint::make_point(0.,1.,0.)));
  testData.push_back( std::make_pair( pt[2], QPoint::make_point(0.,0.,1.)));

  // Test the three edge midpoints
  testData.push_back( std::make_pair(
          QPoint( 0.5 * (pt[0].array() + pt[1].array())),
          QPoint::make_point(0.5,0.5,0.)));
  testData.push_back( std::make_pair(
          QPoint( 0.5 * (pt[0].array() + pt[2].array())),
          QPoint::make_point(0.5,0.,0.5)));
  testData.push_back( std::make_pair(
          QPoint( 0.5 * (pt[1].array() + pt[2].array())),
          QPoint::make_point(0.,0.5,0.5)));

  // Test the triangle midpoint
  testData.push_back( std::make_pair(
          QPoint( 1./3. * (pt[0].array() + pt[1].array() + pt[2].array())),
          QPoint::make_point(1./3.,1./3.,1./3.)));


  // Now run the actual tests
  for(TestVec::const_iterator it= testData.begin(); it != testData.end(); ++it)
  {
    const QPoint& query = it->first;
    const QPoint& expBary = it->second;
    QPoint bary = tri.computeBarycenterCoords(query);

    SLIC_DEBUG(fmt::format(
            "Computed barycentric coordinates for triangle {} and point {} are {}",
            tri, query, bary));
    EXPECT_NEAR(bary[0], expBary[0], EPS );
    EXPECT_NEAR(bary[1], expBary[1], EPS );
    EXPECT_NEAR(bary[2], expBary[2], EPS );
  }

}


//----------------------------------------------------------------------
//----------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using asctoolkit::slic::UnitTestLogger;

int main(int argc, char * argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}



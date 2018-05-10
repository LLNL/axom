/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#include "gtest/gtest.h"

#include "axom/config.hpp"

#include "primal/Point.hpp"
#include "primal/Tetrahedron.hpp"

#include "fmt/format.h"
#include "slic/slic.hpp"

#include <cmath>

using namespace axom;


//------------------------------------------------------------------------------
TEST( primal_tetrahedron, basic)
{
  const int DIM = 3;

  typedef double CoordType;
  typedef primal::Point< CoordType, DIM > QPoint;
  typedef primal::Tetrahedron< CoordType, DIM > QTet;

  const QPoint pt[4] = {
    QPoint::make_point( 0,0,0),
    QPoint::make_point( 1,0,0),
    QPoint::make_point( 1,1,0),
    QPoint::make_point( 1,1,1)
  };

  // Test tetrahedron construction
  QTet tet(pt[0],pt[1],pt[2],pt[3]);

  // Test ostream operator
  SLIC_INFO("Tetrahedron coordinates: " << tet);

  // Check indirection operator
  EXPECT_EQ(pt[0],tet[0]);
  EXPECT_EQ(pt[1],tet[1]);
  EXPECT_EQ(pt[2],tet[2]);
  EXPECT_EQ(pt[3],tet[3]);
}


TEST( primal_tetrahedron, tetrahedron_barycentric)
{
  const int DIM = 3;
  const double EPS = 1e-12;

  typedef double CoordType;
  typedef primal::Point< CoordType, DIM > QPoint;
  typedef primal::Point< CoordType, 4 > RPoint;
  typedef primal::Tetrahedron< CoordType, DIM > QTet;

  const QPoint pt[4] = {
    QPoint::make_point( 1,0,0),
    QPoint::make_point( 0,1,0),
    QPoint::make_point( 0,0,1),
    QPoint::make_point( 0,0,0)
  };

  QTet tet(pt[0],pt[1],pt[2],pt[3]);

  typedef std::vector< std::pair< QPoint,RPoint > > TestVec;
  TestVec testData;

  // Test the three vertices
  const CoordType coord_list0[] = {1.,0.,0.,0.};
  testData.push_back(  std::make_pair( pt[0], RPoint( coord_list0, 4 ) ));
  const CoordType coord_list1[] = {0.,1.,0.,0.};
  testData.push_back(  std::make_pair( pt[1], RPoint( coord_list1, 4 ) ));
  const CoordType coord_list2[] = {0.,0.,1.,0.};
  testData.push_back(  std::make_pair( pt[2], RPoint( coord_list2, 4 ) ));
  const CoordType coord_list3[] = {0.,0.,0.,1.};
  testData.push_back(  std::make_pair( pt[3], RPoint( coord_list3, 4 ) ));

  // Test the edge midpoints
  const CoordType coord_list4[] = {0.5,0.5,0.,0.};
  testData.push_back(  std::make_pair(
                         QPoint( 0.5 * (pt[0].array() + pt[1].array())),
                         RPoint( coord_list4, 4 )));
  const CoordType coord_list5[] = {0.,0.5,0.5,0.};
  testData.push_back( std::make_pair(
                        QPoint( 0.5 * (pt[1].array() + pt[2].array())),
                        RPoint( coord_list5, 4 )));
  const CoordType coord_list6[] = {0.,0.,0.5,0.5};
  testData.push_back( std::make_pair(
                        QPoint( 0.5 * (pt[2].array() + pt[3].array())),
                        RPoint( coord_list6, 4 )));

  // Test the midpoint
  const CoordType coord_list7[] = {1./4.,1./4.,1./4.,1./4.};
  testData.push_back( std::make_pair(
                        QPoint( 1./4. *
                                (pt[0].array() + pt[1].array() +
                                 pt[2].array() + pt[3].array())),
                        RPoint( coord_list7, 4 )));

  // Test a point outside the tetrahedron
  const CoordType coord_list8[] = {-0.4,1.2,0.2,0.};
  testData.push_back(std::make_pair(
                       QPoint(-0.4*pt[0].array() + 1.2*pt[1].array() + 0.2*
                              pt[2].array()),
                       RPoint( coord_list8, 4 )));


  // Now run the actual tests
  for (TestVec::const_iterator it= testData.begin() ;
       it != testData.end() ;
       ++it)
  {
    const QPoint& query = it->first;
    const RPoint& expBary = it->second;
    RPoint bary = tet.physToBarycentric(query);

    SLIC_DEBUG(fmt::format(
                 "Computed barycentric coordinates for tetrahedron {} and point {} are {}",
                 tet, query, bary));

    EXPECT_NEAR(bary[0],  expBary[0], EPS );
    EXPECT_NEAR(bary[1],  expBary[1], EPS );
    EXPECT_NEAR(bary[2],  expBary[2], EPS );
    EXPECT_NEAR(bary[3],  expBary[3], EPS );
  }

}


//----------------------------------------------------------------------
//----------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,
  axom::slic::setLoggingMsgLevel( axom::slic::message::Info);

  int result = RUN_ALL_TESTS();
  return result;
}

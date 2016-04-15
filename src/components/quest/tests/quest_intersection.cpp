/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

/*!
 *******************************************************************************
 * \file
 *
 * \date Jan 5, 2016
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */

#include "gtest/gtest.h"

#include "quest/BoundingBox.hpp"
#include "quest/Intersection.hpp"
#include "quest/Point.hpp"
#include "quest/Ray.hpp"
#include "quest/Segment.hpp"
#include "quest/Triangle.hpp"
#include "quest/Vector.hpp"


TEST( quest_intersection, ray_segment_intersection )
{
  typedef quest::Point< double,2 >   PointType;
  typedef quest::Segment< double,2 > SegmentType;
  typedef quest::Vector< double,2 >  VectorType;
  typedef quest::Ray< double,2 >     RayType;

  // STEP 0: construct segment
  PointType A(0.0);
  PointType B(1.0,1);
  SegmentType S( A, B );

  // STEP 1: construct ray
  PointType origin = PointType::make_point( 0.5,-0.5 );
  VectorType direction;
  direction[0] = 0.0;
  direction[1] = 0.5;
  RayType R( origin,direction.unitVector() );

  // STEP 2: compute intersection
  PointType ip;
  bool intersects = quest::intersect( R, S, ip );
  EXPECT_TRUE( intersects );
  EXPECT_DOUBLE_EQ(0.5,ip[0]);
  EXPECT_DOUBLE_EQ(0.0,ip[1]);

  // STEP 3: construct non-intersecting ray
  origin[1] = 0.5; // shift R up
  RayType R2( origin, direction.unitVector() );
  bool intersects2 = quest::intersect( R2, S, ip );
  EXPECT_FALSE( intersects2 );
}


TEST( quest_intersection, triangle_aabb_intersection )
{
  static int const DIM = 3;
  typedef quest::Point< double,DIM >   PointType;
  typedef quest::Triangle< double,DIM > TriangleType;
  typedef quest::BoundingBox< double,DIM > BoundingBoxType;

  double xArr[3] = { 1., 0., 0.};
  double yArr[3] = { 0., 1., 0.};
  double zArr[3] = { 0., 0., 1.};

  PointType ptX(xArr);
  PointType ptY(yArr);
  PointType ptZ(zArr);

  TriangleType unitTri( ptX, ptY, ptZ );
  BoundingBoxType unitBB( PointType::zero(), PointType::ones());

  EXPECT_TRUE( quest::intersect(unitTri, unitBB));

  // Let's first move the bounding box around
  BoundingBoxType v0_BB( ptX );
  v0_BB.expand(.1);
  SLIC_INFO("Testing v0 bounding box: " << v0_BB << " against unit triangle");
  EXPECT_TRUE( v0_BB.contains(ptX));
  EXPECT_TRUE( quest::intersect(unitTri, v0_BB) );

  BoundingBoxType v1_BB( ptY );
  v1_BB.expand(.1);
  SLIC_INFO("Testing v1 bounding box: " << v1_BB << " against unit triangle");
  EXPECT_TRUE( v1_BB.contains(ptY));
  EXPECT_TRUE( quest::intersect(unitTri, v1_BB) );

  BoundingBoxType v2_BB( ptZ );
  v2_BB.expand(.1);
  SLIC_INFO("Testing v2 bounding box: " << v2_BB << " against unit triangle");
  EXPECT_TRUE( v2_BB.contains(ptZ));
  EXPECT_TRUE( quest::intersect(unitTri, v2_BB) );


  BoundingBoxType mid_BB( PointType::zero());
  mid_BB.addPoint( PointType(0.9));
  SLIC_INFO("Testing bounding box: " << mid_BB << " against unit triangle.  Note -- BB should intersect interior of triangle");
  EXPECT_TRUE( quest::intersect(unitTri, mid_BB) );

  BoundingBoxType high_BB( PointType::ones());
  high_BB.addPoint( PointType(0.5));
  SLIC_INFO("Testing bounding box: " << high_BB << " against unit triangle.  Note -- BB should not intersect interior of triangle");
  EXPECT_FALSE( quest::intersect(unitTri, high_BB) );

  BoundingBoxType out_BB( PointType::ones());
  out_BB.addPoint( PointType(2));
  SLIC_INFO("Testing bounding box: " << out_BB << " against unit triangle.  Note -- BB should not intersect triangle");
  EXPECT_FALSE( quest::intersect(unitTri, out_BB) );


  BoundingBoxType negBB(PointType(-5), PointType(-10));
  SLIC_INFO("Testing bounding box: " << negBB << " against unit triangle.  Note -- BB should not intersect triangle");
  EXPECT_FALSE( quest::intersect(unitTri, negBB) );


  // Test new triangle whose edge crosses the BB
  double t2_0[3] = { 10., 0., 0.};
  double t2_1[3] = { -10.,0., 0.};
  double t2_2[3] = { 0., 100., 0};

  TriangleType xyTri( t2_0, t2_1, t2_2);
  BoundingBoxType bbOrigin(PointType::zero() );
  bbOrigin.expand(1.);
  SLIC_INFO("Testing bounding box: " << bbOrigin << " against triangle " << xyTri << ".  Note -- BB should not intersect triangle");
  EXPECT_TRUE( quest::intersect(xyTri, bbOrigin) );


  BoundingBoxType bbOrigin2(PointType::zero() );
  bbOrigin.addPoint( PointType(-1.));
  bbOrigin.addPoint( PointType::make_point(-1.,1.,1.));
  SLIC_INFO("Testing bounding box: " << bbOrigin2<< " against triangle " << xyTri << ".  Note -- BB should not intersect triangle");
  EXPECT_TRUE( quest::intersect(xyTri, bbOrigin2) );

  BoundingBoxType bbAbove(PointType::ones() );
  bbAbove.addPoint( PointType(2.));
  SLIC_INFO("Testing bounding box: " << bbAbove << " against triangle " << xyTri << ".  Note -- BB should not intersect triangle");
  EXPECT_FALSE( quest::intersect(xyTri, bbAbove) );

  BoundingBoxType bbBelow;
  bbBelow.addPoint( PointType(-1.));
  bbBelow.addPoint( PointType(-2.));
  SLIC_INFO("Testing bounding box: " << bbBelow << " against triangle " << xyTri << ".  Note -- BB should not intersect triangle");
  EXPECT_FALSE( quest::intersect(xyTri, bbBelow) );

  BoundingBoxType bbPoint_OnTri;
  bbPoint_OnTri.addPoint( PointType::make_point(0.,1.,0.));
  SLIC_INFO("Testing point bounding box: " << bbPoint_OnTri << " against triangle " << xyTri << ".  Note -- BB is a point on triangle");
  EXPECT_TRUE( quest::intersect(xyTri, bbPoint_OnTri) );

  BoundingBoxType bbPoint_OutsideTri;
  bbPoint_OutsideTri.addPoint( PointType::make_point(1.,1.,1.));
  SLIC_INFO("Testing point bounding box: " << bbPoint_OutsideTri << " against triangle " << xyTri << ".  Note -- BB is a point outside triangle");
  EXPECT_FALSE( quest::intersect(xyTri, bbPoint_OutsideTri) );

  BoundingBoxType bbInvalid;
  SLIC_INFO("Testing point bounding box: " << bbInvalid << " against triangle " << xyTri << ".  Note -- BB is invalid (empty)");
  EXPECT_FALSE( quest::intersect(xyTri, bbInvalid) );
}



TEST( quest_intersection, triangle_aabb_intersection_fromData )
{
    static int const DIM = 3;
    typedef quest::Point< double,DIM >   PointType;
    typedef quest::Triangle< double,DIM > TriangleType;
    typedef quest::BoundingBox< double,DIM > BoundingBoxType;


    PointType v0 = PointType::make_point(-31.015,63.7756,55.0043);
    PointType v1 = PointType::make_point(-29.0086,59.2982,58.0078);
    PointType v2 = PointType::make_point(-29.2009,70.1039,61.3229);

    TriangleType tri(v0,v1,v2);

    BoundingBoxType box0(PointType::make_point(-39.2793,46.3735,53.3791), PointType::make_point(-26.1692,60.1549,57.0148));
    BoundingBoxType box1(PointType::make_point(-39.2793,60.1549,53.3791), PointType::make_point(-26.1692,73.9362,57.0148));
    BoundingBoxType box2(PointType::make_point(-39.2793,46.3735,57.0148), PointType::make_point(-26.1692,60.1549,60.6506));
    BoundingBoxType box3(PointType::make_point(-39.2793,60.1549,57.0148), PointType::make_point(-26.1692,73.9362,60.6506));
    BoundingBoxType box4(PointType::make_point(-39.2793,46.3735,60.6506), PointType::make_point(-26.1692,60.1549,64.2863));
    BoundingBoxType box5(PointType::make_point(-39.2793,60.1549,60.6506), PointType::make_point(-26.1692,73.9362,64.2863));

    SLIC_INFO("Testing point bounding box: " << box0 << " against triangle " << tri );
    EXPECT_FALSE( quest::intersect(tri, box0 ));

    SLIC_INFO("Testing point bounding box: " << box1 << " against triangle " << tri );
    EXPECT_TRUE( quest::intersect(tri, box1 ));

    //
    asctoolkit::slic::setLoggingMsgLevel( asctoolkit::slic::message::Debug);

    SLIC_INFO("Testing point bounding box: " << box2 << " against triangle " << tri );
    EXPECT_TRUE( quest::intersect(tri, box2 ));

    asctoolkit::slic::setLoggingMsgLevel( asctoolkit::slic::message::Warning);

    SLIC_INFO("Testing point bounding box: " << box3 << " against triangle " << tri );
    EXPECT_TRUE( quest::intersect(tri, box3 ));

    SLIC_INFO("Testing point bounding box: " << box4 << " against triangle " << tri );
    EXPECT_FALSE( quest::intersect(tri, box4 ));

    SLIC_INFO("Testing point bounding box: " << box5 << " against triangle " << tri );
    EXPECT_TRUE( quest::intersect(tri, box5 ));


}

TEST( quest_intersection, triangle_aabb_intersection_fromData2 )
{
    static int const DIM = 3;
    typedef quest::Point< double,DIM >   PointType;
    typedef quest::Triangle< double,DIM > TriangleType;
    typedef quest::BoundingBox< double,DIM > BoundingBoxType;

    // Triangle 569
    TriangleType tri(PointType::make_point(0,5,0), PointType::make_point(-0.665356,4.93844,-0.411212), PointType::make_point(-0.665356,4.93844,0.411212));

    // {pt: (8,15,8); level: 4}
    BoundingBoxType box0(PointType::make_point(0,4.375,0), PointType::make_point(0.625,5,0.625));
    // {pt: (6,15,7); level: 4}
    BoundingBoxType box1(PointType::make_point(-1.25,4.375,-0.625), PointType::make_point(-0.625,5,0));
    // {pt: (6,15,8); level: 4}
    BoundingBoxType box2(PointType::make_point(-1.25,4.375,0), PointType::make_point(-0.625,5,0.625));

    // Block index {pt: (16,31,16); level: 5}
    BoundingBoxType box3(PointType::make_point(0,4.6875,0), PointType::make_point(0.3125,5,0.3125));

    // Block index {pt: (8,15,8); level: 4}
    BoundingBoxType box4(PointType::make_point(0,4.375,0), PointType::make_point(0.625,5,0.625));

    asctoolkit::slic::setLoggingMsgLevel( asctoolkit::slic::message::Info);


    SLIC_INFO("Testing point bounding box: " << box0 << " against triangle " << tri );
    EXPECT_TRUE( quest::intersect(tri, box0));

    SLIC_INFO("Testing point bounding box: " << box1 << " against triangle " << tri );
    EXPECT_TRUE( quest::intersect(tri, box1));

    SLIC_INFO("Testing point bounding box: " << box2 << " against triangle " << tri );
    EXPECT_TRUE( quest::intersect(tri, box2));

    SLIC_INFO("Testing point bounding box: " << box3 << " against triangle " << tri );
    EXPECT_TRUE( quest::intersect(tri, box3));

    SLIC_INFO("Testing point bounding box: " << box4 << " against triangle " << tri );
    EXPECT_TRUE( quest::intersect(tri, box4));

    asctoolkit::slic::setLoggingMsgLevel( asctoolkit::slic::message::Warning);
}


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using asctoolkit::slic::UnitTestLogger;

int main(int argc, char * argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,

  asctoolkit::slic::setLoggingMsgLevel( asctoolkit::slic::message::Warning);

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}

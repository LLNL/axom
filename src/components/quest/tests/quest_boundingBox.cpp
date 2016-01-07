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

#include "quest/Point.hpp"
#include "quest/BoundingBox.hpp"

//------------------------------------------------------------------------------
TEST( quest_boundingBox, bb_default_constructor)
{
  static const int DIM = 2;
  typedef double CoordType;
  typedef quest::Point<CoordType, DIM> QPoint;
  typedef quest::BoundingBox<CoordType, DIM> QBBox;

  QBBox bbox;
  EXPECT_FALSE( bbox.isValid() ) << "Default constructed bounding box is invalid";
  EXPECT_FALSE( bbox.contains( QPoint() ))
      << "Default constructed bounding box should not contain any points";
  EXPECT_FALSE( bbox.contains( QPoint(1000) ))
      << "Default constructed bounding box should not contain any points";
}

TEST( quest_boundingBox, bb_ctor_from_singlePt)
{
  static const int DIM = 3;
  typedef double CoordType;
  typedef quest::Point<CoordType, DIM> QPoint;
  typedef quest::BoundingBox<CoordType, DIM> QBBox;

  QPoint pt1;
  QPoint pt2 (2);

  QBBox bbox1(pt1);

  EXPECT_TRUE(bbox1.isValid());
  EXPECT_TRUE(bbox1.contains(pt1));
  EXPECT_FALSE(bbox1.contains(pt2));
  EXPECT_EQ(bbox1.getMin(), bbox1.getMax())
      <<"BBox only has a single point, so bb.getMin()==bb.getMax()";

  QBBox bbox2(pt2);

  EXPECT_TRUE(bbox2.isValid());
  EXPECT_TRUE(bbox2.contains(pt2));
  EXPECT_FALSE(bbox2.contains(pt1));
  EXPECT_EQ(bbox2.getMin(), bbox2.getMax())
      <<"BBox only has a single point, so bb.getMin()==bb.getMax()";

}


TEST( quest_boundingBox, bb_ctor_from_twoPoints)
{
  static const int DIM = 3;
  typedef double CoordType;
  typedef quest::Point<CoordType, DIM> QPoint;
  typedef quest::BoundingBox<CoordType, DIM> QBBox;

  QPoint pt1(1);
  QPoint pt2(3);
  QPoint midPt = QPoint::midpoint(pt1,pt2);
  QPoint outPt(5);

  QBBox bbox1(pt1, pt2);

  EXPECT_TRUE(bbox1.isValid());
  EXPECT_TRUE(bbox1.contains(pt1));
  EXPECT_TRUE(bbox1.contains(pt2));
  EXPECT_TRUE(bbox1.contains(midPt));
  EXPECT_FALSE(bbox1.contains(outPt));

  //
  std::cout<<"\n** Testing from pairs of points that are "
          <<"max and min of the bounding box"<<std::endl;
  QBBox bbox2(pt2, pt1);

  EXPECT_TRUE(bbox2.isValid());
  EXPECT_TRUE(bbox2.contains(pt1));
  EXPECT_TRUE(bbox2.contains(pt2));
  EXPECT_TRUE(bbox2.contains(midPt));
  EXPECT_FALSE(bbox2.contains(outPt));



  //
  std::cout<<"\n** Testing from pairs of points that are "
          <<"not the smallest and largest of the bounding box"<<std::endl;
  const int val = 10;
  QPoint pt101 = QPoint::make_point(val,0,val);
  QPoint pt010 = QPoint::make_point(0,val,0);;
  QPoint midPt2 = QPoint::midpoint(pt101,pt010);

  QBBox bbox3(pt101, pt010);

  EXPECT_TRUE(bbox3.isValid());
  EXPECT_TRUE(bbox3.contains(pt101));
  EXPECT_TRUE(bbox3.contains(pt010));
  EXPECT_TRUE(bbox3.contains(midPt2));
  EXPECT_TRUE(bbox3.contains( QPoint() ));
  EXPECT_TRUE(bbox3.contains( QPoint(val) ));
  EXPECT_TRUE(bbox3.contains( QPoint::make_point(0, val,val) ));

}

TEST( quest_boundingBox, bb_addPoint)
{
  static const int DIM = 3;
  typedef double CoordType;
  typedef quest::Point<CoordType, DIM> QPoint;
  typedef quest::BoundingBox<CoordType, DIM> QBBox;

  QPoint pt1(1);
  QPoint pt2(3);
  QPoint pt10(10);
  QPoint pt20(20);
  QPoint pt30(30);

  QBBox bbox1(pt1);
  bbox1.addPoint(pt2);

  EXPECT_TRUE(bbox1.isValid());
  EXPECT_TRUE(bbox1.contains(pt1));
  EXPECT_TRUE(bbox1.contains(pt2));
  EXPECT_TRUE(bbox1.contains( QPoint::midpoint(pt1,pt2)));

  // Testing that if we add a missing point, it is within bounds
  EXPECT_FALSE(bbox1.contains(pt10));
  bbox1.addPoint(pt10);
  EXPECT_TRUE(bbox1.contains(pt10));

  // Testing that if we add a point, then points outside the bounds remain outside
  EXPECT_FALSE(bbox1.contains(pt30));
  bbox1.addPoint(pt20);                 // note: adding 20, but testing 30
  EXPECT_FALSE(bbox1.contains(pt30));

}

TEST( quest_boundingBox, bb_test_clear)
{
  static const int DIM = 3;
  typedef double CoordType;
  typedef quest::Point<CoordType, DIM> QPoint;
  typedef quest::BoundingBox<CoordType, DIM> QBBox;

  QPoint pt1(1);
  QPoint pt2(3);

  QBBox bbox(pt1,pt2);

  EXPECT_TRUE(bbox.isValid());
  EXPECT_TRUE(bbox.contains(pt1));
  EXPECT_TRUE(bbox.contains(pt2));
  EXPECT_TRUE(bbox.contains( QPoint::midpoint(pt1,pt2)));

  bbox.clear();
  EXPECT_FALSE(bbox.isValid())
          << "After clear() the bounding box should be invalid";
  EXPECT_FALSE(bbox.contains(pt1))
          << "After clear() the bounding box should not contain any points";
  EXPECT_FALSE(bbox.contains(pt2))
          << "After clear() the bounding box should not contain any points";
  EXPECT_FALSE(bbox.contains( QPoint::midpoint(pt1,pt2)))
          << "After clear() the bounding box should not contain any points";

}


TEST( quest_boundingBox, bb_copy_and_assignment)
{
  static const int DIM = 3;
  typedef double CoordType;
  typedef quest::Point<CoordType, DIM> QPoint;
  typedef quest::BoundingBox<CoordType, DIM> QBBox;

  QPoint pt1(1);
  QPoint pt2(3);
  QPoint midPt = QPoint::midpoint(pt1,pt2);

  QBBox bbox1(pt1, pt2);
  EXPECT_TRUE(bbox1.isValid());

  QBBox bbox2(bbox1);               // copy .ctor
  EXPECT_TRUE(bbox2.isValid());

  EXPECT_EQ(bbox1.getMin(), bbox2.getMin());
  EXPECT_EQ(bbox1.getMax(), bbox2.getMax());
  EXPECT_TRUE( bbox2.contains(midPt));


  QBBox bbox3(pt2);                 // some initialization that we don't care about
  EXPECT_TRUE(bbox3.isValid());
  EXPECT_FALSE(bbox3.contains(pt1));
  bbox3 = bbox1;                    // assignment operation

  EXPECT_EQ(bbox1.getMin(), bbox3.getMin());
  EXPECT_EQ(bbox1.getMax(), bbox3.getMax());
  EXPECT_TRUE( bbox3.contains(midPt));



  QBBox bbox4;
  QBBox bbox5;
  EXPECT_EQ(bbox4.getMin(), bbox5.getMin());
  EXPECT_EQ(bbox4.getMax(), bbox5.getMax());
  EXPECT_NE(bbox4.getMin(), bbox1.getMin()) << "Empty bb should not be equal to an initialized one.";
  EXPECT_NE(bbox4.getMax(), bbox1.getMax()) << "Empty bb should not be equal to an initialized one.";

}

TEST( quest_boundingBox, bb_test_equality)
{
    static const int DIM = 3;
    typedef double CoordType;
    typedef quest::Point<CoordType, DIM> QPoint;
    typedef quest::BoundingBox<CoordType, DIM> QBBox;

    QPoint pt1(1);
    QPoint pt2(3);
    QPoint midPt = QPoint::midpoint(pt1,pt2);

    QBBox bbox1(pt1, pt2);
    EXPECT_TRUE(bbox1.isValid());
    EXPECT_TRUE(bbox1.contains(pt1));
    EXPECT_TRUE(bbox1.contains(pt2));

    QBBox bbox2;
    EXPECT_FALSE(bbox2.isValid()) << "Default constructed bbox is invalid";
    EXPECT_FALSE(bbox1 == bbox2) << "Default constructed bbox should differ from valid bbox (operator==)";
    EXPECT_TRUE(bbox1 != bbox2) << "Default constructed bbox should differ from valid bbox (operator!=)";

    bbox2.addPoint(midPt);
    EXPECT_TRUE(bbox2.isValid()) << "bbox should be valid after adding a point";
    EXPECT_FALSE(bbox1 == bbox2);
    EXPECT_TRUE(bbox1 != bbox2);

    bbox2.addPoint(pt2);
    EXPECT_TRUE(bbox2.isValid()) << "bbox should be valid after adding a point";
    EXPECT_FALSE(bbox1 == bbox2);
    EXPECT_TRUE(bbox1 != bbox2);

    bbox2.addPoint(pt1);
    EXPECT_TRUE(bbox1 == bbox2) << "after adding both endpoints, bounds should now be equal (operator==)";
    EXPECT_FALSE(bbox1 != bbox2) << "after adding both endpoints, bounds should now be equal (operator!=)";
}

TEST( quest_boundingBox, bb_add_box)
{
    static const int DIM = 3;
    typedef double CoordType;
    typedef quest::Point<CoordType, DIM> QPoint;
    typedef quest::BoundingBox<CoordType, DIM> QBBox;

    //
    std::cout <<"Testing addBox() for two simple bounding boxes"<<std::endl;
    QBBox bbox1 ( QPoint(1), QPoint(3) );       // first box
    EXPECT_TRUE(bbox1.isValid());
    EXPECT_TRUE(bbox1.contains(QPoint(2)));
    EXPECT_FALSE(bbox1.contains(QPoint(4)));

    QBBox bbox2 ( QPoint(5), QPoint(7) );       // second box
    EXPECT_TRUE(bbox2.isValid());
    EXPECT_TRUE(bbox2.contains(QPoint(6)));
    EXPECT_FALSE(bbox2.contains(QPoint(4)));

    bbox1.addBox(bbox2);                        // adding the boxes
    EXPECT_TRUE(bbox1.isValid());
    EXPECT_TRUE(bbox1.contains(QPoint(6)))
            << "After addBox() we should now contain points in the other box";
    EXPECT_TRUE(bbox1.contains(QPoint(4)))
            << "After addBox() we should now contain points in between the two bounding boxes";
    EXPECT_FALSE(bbox1.contains(QPoint(10)))
            << "Points outside both ranges are still outside";


    //
    std::cout <<"Testing addBox() for differently arranged boxes"<<std::endl;
    QBBox bbox3 ( QPoint::make_point(1,3,3)
                , QPoint::make_point(3,1,1) );       // first box in 2nd test
    QBBox bbox4 ( QPoint::make_point(7,5,7)
                , QPoint::make_point(5,7,7) );       // second box in 2nd test
    EXPECT_TRUE(bbox3.isValid());
    EXPECT_TRUE(bbox3.contains(QPoint::make_point(2,3,1)));
    EXPECT_FALSE(bbox3.contains(QPoint(4)));

    bbox3.addBox( bbox4);
    EXPECT_TRUE(bbox3.isValid());
    EXPECT_TRUE(bbox3.contains(QPoint::make_point(2,3,1)));
    EXPECT_TRUE(bbox3.contains(QPoint::make_point(6,6,7)));
    EXPECT_TRUE(bbox3.contains(QPoint(4)));


    //
    std::cout <<"Testing addBox() to initially empty box"<<std::endl;
    QBBox bbox5;
    EXPECT_FALSE(bbox5.isValid());
    EXPECT_NE(bbox5, bbox1);
    bbox5.addBox(bbox1);
    EXPECT_EQ( bbox5, bbox1);

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


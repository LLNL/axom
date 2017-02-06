/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

#include <limits>

#include "gtest/gtest.h"

#include "quest/Point.hpp"
#include "quest/BoundingBox.hpp"
#include "quest/VirtualGrid.hpp"

//-----------------------------------------------------------------------------
TEST( quest_virtual_grid, default_constructor)
{   static const int DIM = 3;
    typedef double CoordType;
    typedef quest::Point<CoordType, DIM> QPoint;
    
    quest::VirtualGrid<QPoint,DIM> valid;

    EXPECT_TRUE(valid.getNumBins() == (100 * 100 * 100));//bins = 100^NDIMS
    EXPECT_TRUE(valid.binEmpty(0));

}

TEST( quest_virtual_grid, argument_constructor)
{
    static const int DIM = 3;
    typedef double CoordType;
    typedef quest::Point<CoordType, DIM> QPoint;
    
    QPoint pt1 =QPoint::make_point(0,0,0);
    double step[DIM] = {10,10,10};
    int resolution[DIM] = {10, 10, 10};

    quest::VirtualGrid<QPoint,DIM> valid(pt1,step,resolution);
    EXPECT_TRUE(valid.getNumBins() == (10 * 10 * 10));
    EXPECT_TRUE(valid.binEmpty(0));

}
TEST( quest_virtual_grid, indexing)
{
     static const int DIM = 3;
    typedef double CoordType;
    typedef quest::Point<CoordType, DIM> QPoint;
    
    quest::VirtualGrid<QPoint,DIM> valid;

    
    QPoint pt1 = QPoint::make_point(1.5,0,0);
    EXPECT_TRUE(valid.getIndex(pt1) == 1);
    
    QPoint pt2 = QPoint::make_point(0,1.5,0);
    EXPECT_TRUE(valid.getIndex(pt2) == 100);

}

TEST(quest_virtual_grid, add_stuff){
    static const int DIM = 3;
    typedef double CoordType;
    typedef quest::Point<CoordType, DIM> QPoint;
    typedef quest::BoundingBox<CoordType, DIM> QBBox;
    quest::VirtualGrid<QPoint,DIM> valid;

    QPoint pt1 = QPoint::make_point(1.1,1.1,1.1);
    QBBox bbox1(pt1);

    valid.insert(bbox1,pt1);
    int index = valid.getIndex(pt1);
    EXPECT_TRUE(valid.getBinContents(index).size() == 1);


}

TEST(quest_virtual_grid, delete_stuff){
    static const int DIM = 3;
    typedef double CoordType;
    typedef quest::Point<CoordType, DIM> QPoint;
    typedef quest::BoundingBox<CoordType, DIM> QBBox;
    quest::VirtualGrid<QPoint,DIM> valid;

    QPoint pt1 = QPoint::make_point(1.1,1.1,1.1);
    QBBox bbox1(pt1);

    valid.insert(bbox1,pt1);
    int index = valid.getIndex(pt1);
    valid.clear(index);
    EXPECT_TRUE(valid.binEmpty(index));   
}


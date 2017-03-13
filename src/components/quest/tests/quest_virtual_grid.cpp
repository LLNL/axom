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

#include "primal/Point.hpp"
#include "primal/BoundingBox.hpp"

#include "quest/VirtualGrid.hpp"

//-----------------------------------------------------------------------------
TEST( quest_virtual_grid, point_constructor)
{
    static const int DIM = 3;
    typedef double CoordType;
    typedef axom::primal::Point<CoordType, DIM> QPoint;
    
    QPoint origin = QPoint::make_point(0,0,0);
    const int spacing = 10;
    double step[DIM] = {spacing, spacing, spacing};
    const int resolution = 10;
    int res[DIM] = {resolution, resolution, resolution};

    quest::VirtualGrid<QPoint,DIM> valid(origin, step, res);
    EXPECT_TRUE(valid.getNumBins() == (spacing * spacing * spacing));
    EXPECT_TRUE(valid.binEmpty(0));

}

TEST( quest_virtual_grid, double_array_constructor)
{
    static const int DIM = 3;
    typedef double CoordType;
    typedef axom::primal::Point<CoordType, DIM> QPoint;
    
    double origin[DIM] = {0, 0, 0};
    const int spacing = 10;
    double step[DIM] = {spacing, spacing, spacing};
    const int resolution = 10;
    int res[DIM] = {resolution, resolution, resolution};

    quest::VirtualGrid<QPoint,DIM> valid(origin, step, res);
    EXPECT_TRUE(valid.getNumBins() == (spacing * spacing * spacing));
    EXPECT_TRUE(valid.binEmpty(0));

}

TEST( quest_virtual_grid, indexing)
{
    static const int DIM = 3;
    typedef double CoordType;
    typedef axom::primal::Point<CoordType, DIM> QPoint;
    
    double origin[DIM] = {0, 0, 0};
    const int spacing = 1;
    double step[DIM] = {spacing, spacing, spacing};
    const int resolution = 100;
    int res[DIM] = {resolution, resolution, resolution};
    quest::VirtualGrid<QPoint,DIM> valid(origin, step, res);

    // valid has 100 bins in each dimension, and each bin has a
    // width of 1.0.  The bins are laid out in row-major order.
    // So, each increment in the x-dimension increases the bin index by 1,
    // an increment in the y-dimension increases the bin by 100, and an
    // increment in the z-dimension increases the bin by 100*100 = 10000.
    
    QPoint pt1 = QPoint::make_point(1.5,0,0);
    int expectedBin = 1;  // The 0th z-slab, the 0th y-row, the 1st x-bin
    EXPECT_EQ(valid.getBinIndex(pt1), expectedBin);
    
    QPoint pt2 = QPoint::make_point(0,1.5,0);
    expectedBin = 100;  // The 0th z-slab, the 1st y-row, the 0th x-bin
    EXPECT_EQ(valid.getBinIndex(pt2), expectedBin);

    QPoint pt3 = QPoint::make_point(99.5, 0, 99.5);
    // The (0-based) 99th z-slab, the 0th y-row, the 99th x-bin:
    expectedBin = 99*100*100 + 0 + 99;
    EXPECT_EQ(valid.getBinIndex(pt3), expectedBin);

    QPoint pt4 = QPoint::make_point(16.1, 99.6, 89.2);
    // The 89th z-slab, the 99th y-row, the 16th x-bin:
    expectedBin = 89*100*100 + 99*100 + 16;
    EXPECT_EQ(valid.getBinIndex(pt4), expectedBin);

    // Now go outside the grid, and get some invalid ones.

    QPoint pt5 = QPoint::make_point(12.5, 100.1, 0);
    // Above 100 is over the fence.
    expectedBin = quest::VirtualGrid<QPoint,DIM>::INVALID_BIN_INDEX;  
    EXPECT_EQ(valid.getBinIndex(pt5), expectedBin);

    QPoint pt6 = QPoint::make_point(-0.5, 12, 54.3);
    // Below 0 is over (under?) the fence.
    expectedBin = quest::VirtualGrid<QPoint,DIM>::INVALID_BIN_INDEX;  
    EXPECT_EQ(valid.getBinIndex(pt6), expectedBin);
}

TEST(quest_virtual_grid, add_stuff){
    static const int DIM = 3;
    typedef double CoordType;
    typedef axom::primal::Point<CoordType, DIM> QPoint;
    typedef axom::primal::BoundingBox<CoordType, DIM> QBBox;

    double origin[DIM] = {0, 0, 0};
    const int spacing = 1;
    double step[DIM] = {spacing, spacing, spacing};
    const int resolution = 100;
    int res[DIM] = {resolution, resolution, resolution};
    quest::VirtualGrid<QPoint,DIM> valid(origin, step, res);

    QPoint pt1 = QPoint::make_point(2.5, 2.5, 2.5);
    QBBox bbox1(pt1);
    valid.insert(bbox1,pt1);
    int index = valid.getBinIndex(pt1);
    EXPECT_TRUE(valid.getBinContents(index).size() == 1);

    QPoint pt2 = QPoint::make_point(2.1,2.1,2.1);
    QPoint pt3 = QPoint::make_point(4.2,2.9,2.1);
    QBBox bbox2(pt2, pt3);
    valid.insert(bbox2, pt2);
    for (int j = 2; j < 5; ++j) {
      for (int i = 1; i < 6; ++i) {
        QPoint tpoint1 = QPoint::make_point(i + 0.5, j + 0.5, 1.8);
        int tidx = valid.getBinIndex(tpoint1);
        EXPECT_EQ(valid.getBinContents(tidx).size(), 0);

        QPoint tpoint2 = QPoint::make_point(i + 0.5, j + 0.5, 3.2);
        tidx = valid.getBinIndex(tpoint2);
        EXPECT_EQ(valid.getBinContents(tidx).size(), 0);
      }
    }
    for (int i = 1; i < 6; ++i) {
      QPoint tpoint1 = QPoint::make_point(i + 0.5, 1.8, 2.5);
      int tidx = valid.getBinIndex(tpoint1);
      EXPECT_EQ(valid.getBinContents(tidx).size(), 0);

      QPoint tpoint2 = QPoint::make_point(i + 0.5, 3.2, 2.5);
      tidx = valid.getBinIndex(tpoint2);
      EXPECT_EQ(valid.getBinContents(tidx).size(), 0);
    }
    QPoint pt4 = QPoint::make_point(2.3, 2.3, 2.6);
    index = valid.getBinIndex(pt4);
    EXPECT_EQ(valid.getBinContents(index).size(), 2);
    QPoint pt5 = QPoint::make_point(3.4, 2.3, 2.6);
    index = valid.getBinIndex(pt5);
    EXPECT_EQ(valid.getBinContents(index).size(), 1);
    QPoint pt6 = QPoint::make_point(4.7, 2.3, 2.6);
    index = valid.getBinIndex(pt6);
    EXPECT_EQ(valid.getBinContents(index).size(), 1);

}

TEST(quest_virtual_grid, delete_stuff){
    static const int DIM = 3;
    typedef double CoordType;
    typedef axom::primal::Point<CoordType, DIM> QPoint;
    typedef axom::primal::BoundingBox<CoordType, DIM> QBBox;

    double origin[DIM] = {0, 0, 0};
    const int spacing = 1;
    double step[DIM] = {spacing, spacing, spacing};
    const int resolution = 100;
    int res[DIM] = {resolution, resolution, resolution};
    quest::VirtualGrid<QPoint,DIM> valid(origin, step, res);

    QPoint pt1 = QPoint::make_point(1.1,1.1,1.1);
    QBBox bbox1(pt1);

    valid.insert(bbox1,pt1);
    int index = valid.getBinIndex(pt1);
    valid.clear(index);
    EXPECT_TRUE(valid.binEmpty(index));   
}


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
    EXPECT_TRUE(valid.getBinIndex(pt1) == expectedBin);
    
    QPoint pt2 = QPoint::make_point(0,1.5,0);
    expectedBin = 100;  // The 0th z-slab, the 1st y-row, the 0th x-bin
    EXPECT_TRUE(valid.getBinIndex(pt2) == expectedBin);
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

    QPoint pt1 = QPoint::make_point(1.1,1.1,1.1);
    QBBox bbox1(pt1);

    valid.insert(bbox1,pt1);
    int index = valid.getBinIndex(pt1);
    EXPECT_TRUE(valid.getBinContents(index).size() == 1);
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


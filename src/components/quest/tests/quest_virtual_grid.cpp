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
TEST( quest_virtual_grid, bbox_constructor)
{
    static const int DIM = 3;
    typedef double CoordType;
    typedef axom::primal::Point<CoordType, DIM> QPoint;

    double p_max[DIM] = {10, 10, 10};
    double p_min[DIM] = {0, 0, 0};
    const int resolution = 4;
    int res[DIM] = {resolution, resolution, resolution};

    quest::VirtualGrid<QPoint,DIM> valid(p_min, p_max, res);
    EXPECT_EQ(valid.getNumBins(),  resolution * resolution * resolution);
    EXPECT_TRUE(valid.isBinEmpty(0));

}

TEST( quest_virtual_grid, indexing)
{
    static const int DIM = 3;
    typedef double CoordType;
    typedef axom::primal::Point<CoordType, DIM> QPoint;
    
    double p_max[DIM] = {100, 100, 100};
    double p_min[DIM] = {0, 0, 0};
    const int resolution = 100;
    int res[DIM] = {resolution, resolution, resolution};

    quest::VirtualGrid<QPoint,DIM> valid(p_min, p_max, res);

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

// For insertions and deletions, we use these next few helper functions
// to track the count of objects in each bin.  

// Verify the count in each bin against the map maintained "by hand".
template< typename T, int NDIMS >
void checkBinCounts(quest::VirtualGrid<T, NDIMS> & v,
                   std::map<int, int> & bincounts)
{
  int bcount = v.getNumBins();
  for (int i = 0; i < bcount; ++i) {
    bool binAgrees =
      (bincounts.count(i) < 1 && v.isBinEmpty(i)) ||
      (bincounts[i] == ((int)v.getBinContents(i).size()));
    EXPECT_TRUE(binAgrees) << "Difference at bin " << i << ": v has " <<
      v.getBinContents(i).size() << " and bincounts has " <<
      (((int)bincounts.count(i)) < 1 ? 0 : bincounts[i]);
  }
}

// Increment the count for a bin
void incr(std::map<int, int> & m, int idx)
{
  int dat = 0;
  if (m.count(idx) > 0) {
    dat = m[idx];
  }
  m[idx] = dat + 1;
}

// Zero out the count for a bin
void zero(std::map<int, int> & m, int idx)
{
  m.erase(idx);
}

TEST(quest_virtual_grid, add_stuff_3D){
    static const int DIM = 3;
    typedef double CoordType;
    typedef axom::primal::Point<CoordType, DIM> QPoint;
    typedef axom::primal::BoundingBox<CoordType, DIM> QBBox;

    double origin[DIM] = {0, 0, 0};
    const int mpt = 6;
    double maxpoint[DIM] = {mpt, mpt, mpt};
    const int resolution = 6;
    int res[DIM] = {resolution, resolution, resolution};
    quest::VirtualGrid<QPoint,DIM> valid(origin, maxpoint, res);

    std::map<int, int> check;

    QPoint pt1 = QPoint::make_point(2.5, 2.5, 2.5);
    QBBox bbox1(pt1);
    valid.insert(bbox1,pt1);
    int index = valid.getBinIndex(pt1);
    incr(check, index);
    {
      SCOPED_TRACE("One bin");
      checkBinCounts(valid, check);
    }

    // Add a point with a bounding box overlapping the first one
    QPoint pt2 = QPoint::make_point(2.1,2.1,2.1);
    QPoint pt3 = QPoint::make_point(4.2,2.9,2.1);
    QBBox bbox2(pt2, pt3);
    valid.insert(bbox2, pt2);
    incr(check, valid.getBinIndex(pt2));
    incr(check, valid.getBinIndex(QPoint::make_point(3.5, 2.5, 2.1)));
    incr(check, valid.getBinIndex(pt3));
    {
      SCOPED_TRACE("Added an overlapping bin");
      checkBinCounts(valid, check);
    }

    // Try inserting something partially outside the index
    QPoint pt4 = QPoint::make_point(2.7, 5.2, -3.1);
    QPoint pt5 = QPoint::make_point(3.1, 9.0, 2.8);
    QBBox bbox3(pt4, pt5);
    valid.insert(bbox3, pt4);
    for (int k = 0; k < 3; ++k) {
      for (int i = 2; i < 4; ++i) {
        incr(check, valid.getBinIndex(QPoint::make_point(i + 0.5, 5.5, k + 0.5)));
      }
    }
    {
      SCOPED_TRACE("Bbox partially overlaps grid");
      checkBinCounts(valid, check);
    }

    // Try inserting something completely outside the index
    QPoint pt6 = QPoint::make_point(-2, 7.2, -3.1);
    QPoint pt7 = QPoint::make_point(-0.1, 9.0, 2.8);
    QBBox bbox4(pt6, pt7);
    valid.insert(bbox4, pt6);
    {
      SCOPED_TRACE("Bbox completely misses grid");
      checkBinCounts(valid, check);
    }

    // Try inserting something completely overlapping the box
    QPoint pt8 = QPoint::make_point(-2, -0.2, -3.1);
    QPoint pt9 = QPoint::make_point(7, 6.2, 6.1);
    QBBox bbox5(pt8, pt9);
    valid.insert(bbox5, pt8);
    for (int i = 0; i < valid.getNumBins(); ++i) {
      incr(check, i);
    }
    {
      SCOPED_TRACE("Bbox completely overlaps grid");
      checkBinCounts(valid, check);
    }
}

TEST(quest_virtual_grid, delete_stuff_3D){
    static const int DIM = 3;
    typedef double CoordType;
    typedef axom::primal::Point<CoordType, DIM> QPoint;
    typedef axom::primal::BoundingBox<CoordType, DIM> QBBox;

    double origin[DIM] = {0, 0, 0};
    const int mpt = 6;
    double maxpoint[DIM] = {mpt, mpt, mpt};
    const int resolution = 6;
    int res[DIM] = {resolution, resolution, resolution};
    quest::VirtualGrid<QPoint,DIM> valid(origin, maxpoint, res);

    std::map<int, int> check;

    // Insert something small (one bin), then clear its bin
    QPoint pt1 = QPoint::make_point(1.1,1.1,1.1);
    QBBox bbox1(pt1);
    valid.insert(bbox1,pt1);
    int index = valid.getBinIndex(pt1);
    incr(check, index);
    valid.clear(index);
    zero(check, index);
    {
      SCOPED_TRACE("Inserted one, then deleted it");
      checkBinCounts(valid, check);
    }

    // Insert some objects, clear some bins
    valid.insert(bbox1,pt1);
    incr(check, valid.getBinIndex(pt1));
    QPoint pt2 = QPoint::make_point(0.1, 0.1, -1.1);
    QPoint pt3 = QPoint::make_point(1.9, 2.3, 7.1);
    QBBox bbox2(pt2, pt3);
    valid.insert(bbox2,pt2);
    for (int k = 0; k < 6; ++k) {
      for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 2; ++i) {
          incr(check, valid.getBinIndex(QPoint::make_point(i + 0.5, j + 0.5, k + 0.5)));
        }
      }
    }
    {
      SCOPED_TRACE("Insert two, clear several bins");
      QPoint test = QPoint::make_point(0.5, 2.5, 5.5);
      int tidx = valid.getBinIndex(test);
      valid.clear(tidx);
      zero(check, tidx);
      test = QPoint::make_point(0.5, 1.5, 0.5);
      tidx = valid.getBinIndex(test);
      valid.clear(tidx);
      zero(check, tidx);
      test = QPoint::make_point(0.5, 2.5, 5.5);
      tidx = valid.getBinIndex(test);
      valid.clear(tidx);
      zero(check, tidx);

      checkBinCounts(valid, check);
    }
    {
      SCOPED_TRACE("Try clearing the invalid bin");
      valid.clear(quest::VirtualGrid<QPoint,DIM>::INVALID_BIN_INDEX);
      checkBinCounts(valid, check);
    }
}

TEST(quest_virtual_grid, add_stuff_2D){
    static const int DIM = 2;
    typedef double CoordType;
    typedef axom::primal::Point<CoordType, DIM> QPoint;
    typedef axom::primal::BoundingBox<CoordType, DIM> QBBox;

    double origin[DIM] = {0, 0};
    const int mpt = 6;
    double maxpoint[DIM] = {mpt, mpt};
    const int resolution = 6;
    int res[DIM] = {resolution, resolution};
    quest::VirtualGrid<QPoint,DIM> valid(origin, maxpoint, res);

    std::map<int, int> check;

    QPoint pt1 = QPoint::make_point(2.5, 2.5);
    QBBox bbox1(pt1);
    valid.insert(bbox1,pt1);
    int index = valid.getBinIndex(pt1);
    incr(check, index);
    {
      SCOPED_TRACE("One bin");
      checkBinCounts(valid, check);
    }

    // Add a point with a bounding box overlapping the first one
    QPoint pt2 = QPoint::make_point(2.1,2.1);
    QPoint pt3 = QPoint::make_point(4.2,2.9);
    QBBox bbox2(pt2, pt3);
    valid.insert(bbox2, pt2);
    incr(check, valid.getBinIndex(pt2));
    incr(check, valid.getBinIndex(QPoint::make_point(3.5, 2.5)));
    incr(check, valid.getBinIndex(pt3));
    {
      SCOPED_TRACE("Added an overlapping bin");
      checkBinCounts(valid, check);
    }

    // Try inserting something partially outside the index
    QPoint pt4 = QPoint::make_point(2.7, 5.2);
    QPoint pt5 = QPoint::make_point(3.1, 9.0);
    QBBox bbox3(pt4, pt5);
    valid.insert(bbox3, pt4);
    for (int i = 2; i < 4; ++i) {
      incr(check, valid.getBinIndex(QPoint::make_point(i + 0.5, 5.5)));
    }
    {
      SCOPED_TRACE("Bbox partially overlaps grid");
      checkBinCounts(valid, check);
    }

    // Try inserting something completely outside the index
    QPoint pt6 = QPoint::make_point(-2, 7.2);
    QPoint pt7 = QPoint::make_point(-0.1, 9.0);
    QBBox bbox4(pt6, pt7);
    valid.insert(bbox4, pt6);
    {
      SCOPED_TRACE("Bbox completely misses grid");
      checkBinCounts(valid, check);
    }

    // Try inserting something completely overlapping the box
    QPoint pt8 = QPoint::make_point(-2, -0.2);
    QPoint pt9 = QPoint::make_point(7, 6.2);
    QBBox bbox5(pt8, pt9);
    valid.insert(bbox5, pt8);
    for (int i = 0; i < valid.getNumBins(); ++i) {
      incr(check, i);
    }
    {
      SCOPED_TRACE("Bbox completely overlaps grid");
      checkBinCounts(valid, check);
    }
}

TEST(quest_virtual_grid, delete_stuff_2D){
    static const int DIM = 2;
    typedef double CoordType;
    typedef axom::primal::Point<CoordType, DIM> QPoint;
    typedef axom::primal::BoundingBox<CoordType, DIM> QBBox;

    double origin[DIM] = {0, 0};
    const int mpt = 6;
    double maxpoint[DIM] = {mpt, mpt};
    const int resolution = 6;
    int res[DIM] = {resolution, resolution};
    quest::VirtualGrid<QPoint,DIM> valid(origin, maxpoint, res);

    std::map<int, int> check;

    // Insert something small (one bin), then clear its bin
    QPoint pt1 = QPoint::make_point(1.1, 1.1);
    QBBox bbox1(pt1);
    valid.insert(bbox1,pt1);
    int index = valid.getBinIndex(pt1);
    incr(check, index);
    valid.clear(index);
    zero(check, index);
    {
      SCOPED_TRACE("Inserted one, then deleted it");
      checkBinCounts(valid, check);
    }

    // Insert some objects, clear some bins
    valid.insert(bbox1,pt1);
    incr(check, valid.getBinIndex(pt1));
    QPoint pt2 = QPoint::make_point(0.1, 0.1);
    QPoint pt3 = QPoint::make_point(1.9, 2.3);
    QBBox bbox2(pt2, pt3);
    valid.insert(bbox2,pt2);
    for (int j = 0; j < 3; ++j) {
      for (int i = 0; i < 2; ++i) {
        incr(check, valid.getBinIndex(QPoint::make_point(i + 0.5, j + 0.5)));
      }
    }
    {
      SCOPED_TRACE("Insert two, clear several bins");
      QPoint test = QPoint::make_point(0.5, 2.5);
      int tidx = valid.getBinIndex(test);
      valid.clear(tidx);
      zero(check, tidx);
      test = QPoint::make_point(0.5, 1.5);
      tidx = valid.getBinIndex(test);
      valid.clear(tidx);
      zero(check, tidx);
      test = QPoint::make_point(0.5, 2.5);
      tidx = valid.getBinIndex(test);
      valid.clear(tidx);
      zero(check, tidx);

      checkBinCounts(valid, check);
    }
    {
      SCOPED_TRACE("Try clearing the invalid bin");
      valid.clear(quest::VirtualGrid<QPoint,DIM>::INVALID_BIN_INDEX);
      checkBinCounts(valid, check);
    }
}


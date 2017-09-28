/*
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


#include "gtest/gtest.h"

#include "primal/Point.hpp"
#include "primal/BoundingBox.hpp"

#include "quest/ImplicitGrid.hpp"

#include "slic/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

#include <cstdlib>    // for srand
#include <vector>
#include <algorithm>  // for std::find


namespace implicit_grid_1D
{
  const int DIM = 1;

  typedef axom::quest::ImplicitGrid<DIM> GridT;

  typedef GridT::GridCell GridCell;
  typedef GridT::SpacePoint SpacePt;
  typedef GridT::SpatialBoundingBox BBox;
}

namespace implicit_grid_2D
{
  const int DIM = 2;

  typedef axom::quest::ImplicitGrid<DIM> GridT;

  typedef GridT::GridCell  GridCell;
  typedef GridT::SpacePoint SpacePt;
  typedef GridT::SpatialBoundingBox BBox;
}

namespace implicit_grid_3D
{
  const int DIM = 3;

  typedef axom::quest::ImplicitGrid<DIM> GridT;

  typedef GridT::GridCell  GridCell;
  typedef GridT::SpacePoint SpacePt;
  typedef GridT::SpatialBoundingBox BBox;
}


TEST( quest_implicit_grid, implicit_grid_ctor_1D)
{
    SLIC_INFO("Test the constructor for a 1D quest::ImplicitGrid object");

    using namespace implicit_grid_1D;

    GridCell res(10);
    BBox bbox(SpacePt::zero(), SpacePt::ones());
    int numElts = 100;

    // Default constructor followed by initialize()
    GridT grid1;
    EXPECT_FALSE( grid1.isInitialized() );

    grid1.initialize(bbox, res, numElts);
    EXPECT_TRUE( grid1.isInitialized() );

    // Initializing constructor
    GridT grid2( bbox, res, numElts);
    EXPECT_TRUE( grid2.isInitialized() );
}

TEST( quest_implicit_grid, implicit_grid_ctor_2D)
{
    SLIC_INFO("Test the constructor for a 2D quest::ImplicitGrid object");

    using namespace implicit_grid_2D;

    GridCell res(10);
    BBox bbox(SpacePt::zero(), SpacePt::ones());
    int numElts = 100;

    // Default constructor followed by initialize()
    GridT grid1;
    EXPECT_FALSE( grid1.isInitialized() );

    grid1.initialize(bbox, res, numElts);
    EXPECT_TRUE( grid1.isInitialized() );

    // Initializing constructor
    GridT grid2( bbox, res, numElts);
    EXPECT_TRUE( grid2.isInitialized() );
}

TEST( quest_implicit_grid, implicit_grid_ctor_3D)
{
    SLIC_INFO("Test the constructor for a 3D quest::ImplicitGrid object");

    using namespace implicit_grid_3D;

    GridCell res(10);
    BBox bbox(SpacePt::zero(), SpacePt::ones());
    int numElts = 100;

    // Default constructor followed by initialize()
    GridT grid1;
    EXPECT_FALSE( grid1.isInitialized() );

    grid1.initialize(bbox, res, numElts);
    EXPECT_TRUE( grid1.isInitialized() );

    // Initializing constructor
    GridT grid2( bbox, res, numElts);
    EXPECT_TRUE( grid2.isInitialized() );
}

TEST( quest_implicit_grid, implicit_grid_insert_contains)
{
    SLIC_INFO("Testing insertion on a quest::ImplicitGrid object");

    using namespace implicit_grid_3D;

    typedef axom::primal::BoundingBox<int, DIM> RangeBox;

    // Note: A 10 x 10 x 10 implicit grid in the unit cube.
    //       Grid cells have a spacing of .1 along each dimension
    GridCell res(10);
    BBox bbox( SpacePt(0.), SpacePt(1.) );
    const int maxElts = 10;

    GridT grid( bbox, res, maxElts);

    /// Test insertion of several objects under the assumption
    /// that we are in a unit cube with a 10^3 grid resolution

    SpacePt pts[] = {
        // Object 1 -- completely inside grid
        SpacePt::make_point(.15, .25, .05),
        SpacePt::make_point(.45, .25, .35),

        // Object 2 -- extends past upper grid boundaries
        SpacePt::make_point(.85, .85, .85),
        SpacePt::make_point(1.25, 1.25, 1.25),

        // Object 3 -- single point at center of grid cell
        SpacePt::make_point(.25, .35, .15),

        // Object 4 -- single point on grid lines
        SpacePt::make_point(.7, .7, .7),

        // Object 5 -- completely outside grid
        SpacePt::make_point(1.85, 1.85, 1.85),
        SpacePt::make_point(2.25, 2.25, 2.25),

        // Object 6 -- overlaps entire grid
        SpacePt::make_point(-1.85, -1.85, -1.85),
        SpacePt::make_point(2.25, 2.25, 2.25)
    };

    // Explicitly define expected ranges of cells
    // where we expect objects to be indexed
    GridCell rangePts[] = {
        GridCell::make_point(1,2,0), // obj1
        GridCell::make_point(4,2,3),

        GridCell::make_point(8,8,8), // obj2
        GridCell::make_point(9,9,9),

        GridCell::make_point(2,3,1), // obj3

        GridCell::make_point(6,6,6), // obj4
        GridCell::make_point(7,7,7)
    };

    const int NUM_DEFINED_OBJS = 7;

    BBox objBox[] = {
        BBox(),     // unused
        BBox( pts[0], pts[1]),
        BBox( pts[2], pts[3]),
        BBox( pts[4]),
        BBox( pts[5]),
        BBox( pts[6], pts[7]),
        BBox( pts[8], pts[9])
    };

    RangeBox ranges[] = {
        RangeBox(),  // unused
        RangeBox( rangePts[0], rangePts[1] ),
        RangeBox( rangePts[2], rangePts[3] ),
        RangeBox( rangePts[4] ),
        RangeBox( rangePts[5], rangePts[6] ),
        RangeBox(),  // unused
        RangeBox(),  // unused
    };

    // Insert objects into implicit grid.
    // Note: We do not insert an object with id 0
    for(int i=1; i< NUM_DEFINED_OBJS; ++i)
    {
      grid.insert( objBox[i], i);
    }

    // Test that points are contained in the expected grid cells
    for(int i=0; i< res[0]; ++i)
    {
      for(int j=0; j< res[1]; ++j)
      {
        for(int k=0; k< res[2]; ++k)
        {
          GridCell gridCell = GridCell::make_point(i,j,k);

          // Item 0 -- was not inserted, so is never in grid
          {
            int idx = 0;
            ASSERT_FALSE( grid.contains(gridCell, idx));
          }

          // Objects 1,2,3 and 4 are expected to be in grid
          // when corresponding ranges contain the query point
          {
            for(int idx=1; idx< 5; ++idx)
            {
              bool exp = ranges[idx].contains(gridCell);
              ASSERT_EQ( exp, grid.contains(gridCell, idx));
            }
          }

          // Object 5 is always outside grid
          {
            int idx = 5;
            EXPECT_FALSE( grid.contains(gridCell, idx));
          }

          // Object 6 is always inside grid
          {
            int idx = 6;
            EXPECT_TRUE( grid.contains(gridCell, idx));
          }

          // Test an index that is out of range
          {
            int idx = maxElts + 1;
            ASSERT_FALSE( grid.contains(gridCell, idx));
          }
        }
      }
    }
}

TEST( quest_implicit_grid, implicit_grid_query)
{
    SLIC_INFO("Query implicit grid using array version");

    using namespace implicit_grid_3D;

    typedef GridT::IndexType IndexType;
    typedef GridT::BitsetType CandidateBitset;
    typedef std::vector<IndexType> CandidateVector;

    // Note: A 10 x 10 x 10 implicit grid in the unit cube.
    //       Grid cells have a spacing of .1 along each dimension
    GridCell res(10);
    BBox bbox(SpacePt(0.), SpacePt(1.));
    int maxElts = 10;

    GridT grid( bbox, res, maxElts);

    BBox objBox1(
        SpacePt::make_point(.15, .25, .05),
        SpacePt::make_point(.45, .25, .35));
    grid.insert(objBox1, 1);

    {
      SpacePt queryPts[] = {
          // First three query points inside only obj1
          objBox1.getMin(),
          objBox1.getMax(),
          objBox1.getCentroid(),

          // Next two points are not in obj1, but in same grid cells
          SpacePt::make_point(0.11, 0.21, 0.01),
          SpacePt::make_point(0.49, 0.29, 0.39),

          // Next four points are not in obj1 or any of same cells
          SpacePt(0.55),
          SpacePt::make_point(.99, 0.25, 0.25),  // outside coord 0
          SpacePt::make_point(.35, 0.99, 0.25),  // outside coord 1
          SpacePt::make_point(.35, 0.25, 0.99),  // outside coord 2
      };

      // Test some points that are expected to match
      for(int i=0; i<5; ++i)
      {
        const std::size_t expSize = 1;
        const IndexType expIdx = 1;

        const SpacePt& queryPt = queryPts[i];

        // Test getCandidates() which returns a bitset
        CandidateBitset candidateBits = grid.getCandidates(queryPt);
        EXPECT_EQ( expSize, candidateBits.count());
        EXPECT_TRUE( candidateBits[expIdx] );

        // Test getCandidatesAsArray() which returns a vector
        CandidateVector candidateVec = grid.getCandidatesAsArray(queryPt);
        EXPECT_EQ( expSize, candidateVec.size() );
        EXPECT_EQ( expIdx, candidateVec[0] );
      }

      // Test some points that are expected to not match
      for(int i=6; i<9; ++i)
      {
        const std::size_t expSize = 0;

        const SpacePt& queryPt = queryPts[i];

        // Test getCandidates() which returns a bitset
        CandidateBitset candidateBits = grid.getCandidates(queryPt);
        EXPECT_EQ( expSize, candidateBits.count());

        // Test getCandidatesAsArray() which returns a vector
        CandidateVector candidateVec = grid.getCandidatesAsArray(queryPt);
        EXPECT_EQ( expSize, candidateVec.size() );
      }
    }

    BBox objBox2(
        SpacePt::make_point(.75, .85, .85),
        SpacePt::make_point(.85, .85, .85));
    grid.insert(objBox2, 2);

    BBox objBox3(
        SpacePt::make_point(.85, .85, .75),
        SpacePt::make_point(.85, .85, .85));
    grid.insert(objBox3, 3);

    SpacePt queryPts[] = {
        SpacePt::make_point(.85, .85, .75), // in obj3 only
        SpacePt::make_point(.85, .85, .85)  // in obj2 and obj3
    };

    // test a point that should contain only obj2
    {
      const std::size_t expSize = 1;
      SpacePt queryPt = SpacePt::make_point(.75, .85, .85);

      CandidateBitset candidateBits = grid.getCandidates(queryPt);
      EXPECT_EQ( expSize, candidateBits.count());
      EXPECT_FALSE( candidateBits[1] );
      EXPECT_TRUE( candidateBits[2] );
      EXPECT_FALSE( candidateBits[3] );

      CandidateVector candidateVec = grid.getCandidatesAsArray(queryPt);
      EXPECT_EQ( expSize, candidateVec.size() );
      EXPECT_EQ( IndexType(2), candidateVec[0] );
    }

    // test a point that should contain only obj3
    {
      const std::size_t expSize = 1;
      SpacePt queryPt = SpacePt::make_point(.85, .85, .75);

      CandidateBitset candidateBits = grid.getCandidates(queryPt);
      EXPECT_EQ( expSize, candidateBits.count());
      EXPECT_FALSE( candidateBits[1] );
      EXPECT_FALSE( candidateBits[2] );
      EXPECT_TRUE( candidateBits[3] );

      CandidateVector candidateVec = grid.getCandidatesAsArray(queryPt);
      EXPECT_EQ( expSize, candidateVec.size() );
      EXPECT_EQ( IndexType(3), candidateVec[0] );
    }

    // test a point that should contain obj2 and obj3, but not obj1
    {
      const std::size_t expSize = 2;
      SpacePt queryPt = SpacePt::make_point(.85, .85, .85);

      CandidateBitset candidateBits = grid.getCandidates(queryPt);
      EXPECT_EQ( expSize, candidateBits.count());
      EXPECT_FALSE( candidateBits[1] );
      EXPECT_TRUE( candidateBits[2] );
      EXPECT_TRUE( candidateBits[3] );

      CandidateVector candidateVec = grid.getCandidatesAsArray(queryPt);
      EXPECT_EQ( expSize, candidateVec.size() );

      bool has1 = candidateVec.end() !=
          std::find(candidateVec.begin(),candidateVec.end(),IndexType(1));
      EXPECT_FALSE( has1 );

      bool has2 = candidateVec.end() !=
          std::find(candidateVec.begin(),candidateVec.end(),IndexType(2));
      EXPECT_TRUE( has2 );

      bool has3 = candidateVec.end() !=
          std::find(candidateVec.begin(),candidateVec.end(),IndexType(3));
      EXPECT_TRUE( has3 );
    }
}


//----------------------------------------------------------------------

int main(int argc, char * argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,
  axom::slic::setLoggingMsgLevel( axom::slic::message::Info );

  std::srand( 105 );

  result = RUN_ALL_TESTS();

  return result;
}

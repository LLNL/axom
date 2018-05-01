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


/*!
 * \file
 * \brief Unit tests for quest's ImplicitGrid class
 *
 * Uses gtest's TYPED_TESTS to test ImplicitGrid in 1,2 and 3 dimensions.
 */

#include "gtest/gtest.h"

#include "primal/Point.hpp"
#include "primal/BoundingBox.hpp"

#include "quest/ImplicitGrid.hpp"

#include "slic/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

#include <vector>
#include <algorithm>  // for std::find


/*!
 * Templated test fixture for ImplicitGrid tests
 *
 * \tparam T Is expected to wrap an integer, in its Value field
 */
template<typename T>
class ImplicitGridTest : public ::testing::Test
{
public:
  typedef T GridT;
  typedef typename GridT::GridCell GridCell;
  typedef typename GridT::SpacePoint SpacePt;
  typedef typename GridT::SpatialBoundingBox BBox;

  static const int DIM = SpacePt::DIMENSION;
};

/*! Typelist for TypedTests on ImplicitGrid */
typedef ::testing::Types <
    axom::quest::ImplicitGrid<1>,
    axom::quest::ImplicitGrid<2>,
    axom::quest::ImplicitGrid<3> > MyTypes;

TYPED_TEST_CASE( ImplicitGridTest, MyTypes );


TYPED_TEST( ImplicitGridTest, implicit_grid_ctor)
{
  const int DIM = TestFixture::DIM;
  typedef typename TestFixture::GridCell GridCell;
  typedef typename TestFixture::BBox BBox;
  typedef typename TestFixture::GridT GridT;
  typedef typename TestFixture::SpacePt SpacePt;

  SLIC_INFO("Test ImplicitGrid constructor in " << DIM << "D");

  GridCell res(10);
  BBox bbox(SpacePt::zero(), SpacePt::ones());
  int numElts = 100;

  // Tests default constructor followed by initialize()
  GridT grid1;
  EXPECT_FALSE( grid1.isInitialized() );

  grid1.initialize(bbox, &res, numElts);
  EXPECT_TRUE( grid1.isInitialized() );


  // Tests initializing constructor
  GridT grid2( bbox, &res, numElts);
  EXPECT_TRUE( grid2.isInitialized() );


  // Tests initializing from primitive types
  GridT grid3( bbox.getMin().data(), bbox.getMax().data(), res.data(), numElts);
  EXPECT_TRUE( grid3.isInitialized() );

}


TYPED_TEST( ImplicitGridTest, implicit_grid_resolution)
{
  const int DIM = TestFixture::DIM;
  typedef typename TestFixture::GridCell GridCell;
  typedef typename TestFixture::BBox BBox;
  typedef typename TestFixture::GridT GridT;
  typedef typename TestFixture::SpacePt SpacePt;

  SLIC_INFO("Test ImplicitGrid resolution in " << DIM << "D");

  BBox bbox(SpacePt::zero(), SpacePt::ones());

  // Set the number of elements so that the DIM^th root is an integer
  const int dimRes = 8;
  int numElts = 1;
  for(int i=0 ; i< DIM ; ++i)
    numElts *= dimRes;

  // Tests explicitly set grid resolution
  {
    GridCell res= GridCell::make_point(12, 8, 16);
    GridT grid( bbox, &res, numElts);
    EXPECT_EQ( res, grid.gridResolution() );
    EXPECT_EQ( numElts, grid.numIndexElements() );
  }

  // Tests implicitly set grid resolution
  {
    GridT grid( bbox, AXOM_NULLPTR, numElts);
    GridCell expRes(dimRes);
    EXPECT_EQ( expRes, grid.gridResolution() );
    EXPECT_EQ( numElts, grid.numIndexElements() );
  }

  // Test that grid resolution is at least one in each dim
  {
    // Using NULL pointer for resolution
    int zeroMeshElts = 0;
    GridT grid( bbox, AXOM_NULLPTR, zeroMeshElts );
    GridCell expRes = GridCell::ones();
    EXPECT_EQ( expRes, grid.gridResolution() );
    EXPECT_EQ( 0, grid.numIndexElements() );
  }
  {
    // Using a resolution of zero in each dim
    GridCell zeroRes = GridCell::zero();
    GridT grid( bbox, &zeroRes, numElts );
    GridCell expRes = GridCell::ones();
    EXPECT_EQ( expRes, grid.gridResolution() );
    EXPECT_EQ( numElts, grid.numIndexElements() );
  }
}


TYPED_TEST( ImplicitGridTest, implicit_grid_insert_contains)
{
  const int DIM = TestFixture::DIM;
  typedef typename TestFixture::GridCell GridCell;
  typedef typename TestFixture::BBox BBox;
  typedef typename TestFixture::GridT GridT;
  typedef typename TestFixture::SpacePt SpacePt;
  typedef axom::primal::BoundingBox<int, DIM> RangeBox;

  SLIC_INFO("Testing ImplicitGrid insert() and contains() in " << DIM << "D");

  // Note: A 10 x 10 x 10 implicit grid in the unit cube.
  //       Grid cells have a spacing of .1 along each dimension

  GridCell res(10);
  BBox bbox( SpacePt(0.), SpacePt(1.) );
  const int maxElts = 10;
  const int numDefinedObjs = 7;

  GridT grid( bbox, &res, maxElts);

  /// Test insertion of several objects under assumption that we're
  /// in a unit cube with 10 grid cells per dimension resolution
  BBox objBox[maxElts];
  RangeBox ranges[maxElts];

  {
    // setup test data
    // Note: make_point() ignores coordinates higher than its dimension
    SpacePt pts[] = {
      // Object 1 -- completely inside grid
      SpacePt::make_point(.15, .25, .05),
      SpacePt::make_point(.45, .25, .35),

      // Object 2 -- extends past upper grid boundaries
      SpacePt::make_point(.85, .85, .85),
      SpacePt::make_point(1.25, 1.25, 1.25),

      // Object 3 -- single point at center of a grid cell
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
      GridCell::make_point(1,2,0),   // obj1
      GridCell::make_point(4,2,3),

      GridCell::make_point(8,8,8),   // obj2
      GridCell::make_point(9,9,9),

      GridCell::make_point(2,3,1),   // obj3

      GridCell::make_point(6,6,6),   // obj4
      GridCell::make_point(7,7,7)
    };

    objBox[0] = BBox();       // unused
    objBox[1] = BBox( pts[0], pts[1]);
    objBox[2] = BBox( pts[2], pts[3]);
    objBox[3] = BBox( pts[4]);
    objBox[4] = BBox( pts[5]);
    objBox[5] = BBox( pts[6], pts[7]);
    objBox[6] = BBox( pts[8], pts[9]);

    ranges[0] = RangeBox();   // unused
    ranges[1] = RangeBox( rangePts[0], rangePts[1] );
    ranges[2] = RangeBox( rangePts[2], rangePts[3] );
    ranges[3] = RangeBox( rangePts[4] );
    ranges[4] = RangeBox( rangePts[5], rangePts[6] );
    ranges[5] = RangeBox();   // unused
    ranges[6] = RangeBox();   // unused
  }

  // Insert objects into implicit grid.
  // Note: We do not insert an object with id 0
  for(int i=1 ; i< numDefinedObjs ; ++i)
  {
    grid.insert( objBox[i], i);
  }

  // Test that points are contained in the expected grid cells
  const int i_max = DIM >= 1 ? grid.gridResolution()[0] : 1;
  const int j_max = DIM >= 2 ? grid.gridResolution()[1] : 1;
  const int k_max = DIM >= 3 ? grid.gridResolution()[2] : 1;
  for(int i=0 ; i< i_max ; ++i)
  {
    for(int j=0 ; j< j_max ; ++j)
    {
      for(int k=0 ; k< k_max ; ++k)
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
          for(int idx=1 ; idx< 5 ; ++idx)
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

TYPED_TEST( ImplicitGridTest, implicit_grid_get_candidates)
{
  const int DIM = TestFixture::DIM;
  typedef typename TestFixture::GridCell GridCell;
  typedef typename TestFixture::BBox BBox;
  typedef typename TestFixture::GridT GridT;
  typedef typename TestFixture::SpacePt SpacePt;

  SLIC_INFO("Test ImplicitGrid getCandidates() in " << DIM << "D");

  typedef typename GridT::IndexType IndexType;
  typedef typename GridT::BitsetType CandidateBitset;
  typedef std::vector<IndexType> CandidateVector;

  // Note: A 10 x 10 x 10 implicit grid in the unit cube.
  //       Grid cells have a spacing of .1 along each dimension
  GridCell res(10);
  BBox bbox(SpacePt(0.), SpacePt(1.));
  const int maxElts = 10;

  GridT grid( bbox, &res, maxElts);

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
      SpacePt::make_point(.99, 0.25, 0.25),    // outside coord 0
      SpacePt::make_point(.35, 0.99, 0.25),    // outside coord 1
      SpacePt::make_point(.35, 0.25, 0.99),    // outside coord 2
    };

    // Test some points that are expected to match
    for(int i=0 ; i<5 ; ++i)
    {
      const std::size_t expSize = 1;
      const IndexType expIdx = 1;

      const SpacePt& queryPt = queryPts[i];

      // Test getCandidates() which returns a bitset
      CandidateBitset candidateBits = grid.getCandidates(queryPt);
      EXPECT_EQ( expSize, candidateBits.count());
      EXPECT_TRUE( candidateBits.test(expIdx) );

      // Test getCandidatesAsArray() which returns a vector
      CandidateVector candidateVec = grid.getCandidatesAsArray(queryPt);
      EXPECT_EQ( expSize, candidateVec.size() );
      EXPECT_EQ( expIdx, candidateVec[0] );
    }

    // Test some points that are expected to not match
    for(int i=5 ; i< 6+DIM ; ++i)
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
    SpacePt(.85));
  grid.insert(objBox2, 2);

  BBox objBox3(
    SpacePt::make_point(.85, .85, .75),
    SpacePt(.95));
  grid.insert(objBox3, 3);

  // test a point that should contain only obj2
  {
    const std::size_t expSize = 1;
    SpacePt queryPt = SpacePt::make_point(.75, .85, .85);

    CandidateBitset candidateBits = grid.getCandidates(queryPt);
    EXPECT_EQ( expSize, candidateBits.count());
    EXPECT_FALSE( candidateBits.test(1) );
    EXPECT_TRUE( candidateBits.test(2) );
    EXPECT_FALSE( candidateBits.test(3) );

    CandidateVector candidateVec = grid.getCandidatesAsArray(queryPt);
    EXPECT_EQ( expSize, candidateVec.size() );
    EXPECT_EQ( IndexType(2), candidateVec[0] );
  }

  // test a point that should contain only obj3
  {
    const std::size_t expSize = 1;
    SpacePt queryPt = SpacePt(.91);

    CandidateBitset candidateBits = grid.getCandidates(queryPt);
    EXPECT_EQ( expSize, candidateBits.count());
    EXPECT_FALSE( candidateBits.test(1) );
    EXPECT_FALSE( candidateBits.test(2) );
    EXPECT_TRUE( candidateBits.test(3) );

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
    EXPECT_FALSE( candidateBits.test(1) );
    EXPECT_TRUE( candidateBits.test(2) );
    EXPECT_TRUE( candidateBits.test(3) );

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

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,
  axom::slic::setLoggingMsgLevel( axom::slic::message::Info );

  result = RUN_ALL_TESTS();

  return result;
}

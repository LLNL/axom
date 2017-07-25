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

#include <cstdlib>
#include <limits>


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

    GridCell res(10, DIM);
    BBox bbox(SpacePt::zero(), SpacePt::ones());
    int numElts = 100;

    GridT grid1;
    EXPECT_FALSE( grid1.isInitialized() );

    grid1.initialize(bbox, res, numElts);
    EXPECT_TRUE( grid1.isInitialized() );


    GridT grid2( bbox, res, numElts);
    EXPECT_TRUE( grid2.isInitialized() );
}

TEST( quest_implicit_grid, implicit_grid_ctor_2D)
{
    SLIC_INFO("Test the constructor for a 2D quest::ImplicitGrid object");

    using namespace implicit_grid_2D;

    GridCell res(10, DIM);
    BBox bbox(SpacePt::zero(), SpacePt::ones());
    int numElts = 100;

    GridT grid1;
    EXPECT_FALSE( grid1.isInitialized() );

    grid1.initialize(bbox, res, numElts);
    EXPECT_TRUE( grid1.isInitialized() );


    GridT grid2( bbox, res, numElts);
    EXPECT_TRUE( grid2.isInitialized() );
}

TEST( quest_implicit_grid, implicit_grid_ctor_3D)
{
    SLIC_INFO("Test the constructor for a 3D quest::ImplicitGrid object");

    using namespace implicit_grid_3D;

    GridCell res(10, DIM);
    BBox bbox(SpacePt::zero(), SpacePt::ones());
    int numElts = 100;

    GridT grid1;
    EXPECT_FALSE( grid1.isInitialized() );

    grid1.initialize(bbox, res, numElts);
    EXPECT_TRUE( grid1.isInitialized() );


    GridT grid2( bbox, res, numElts);
    EXPECT_TRUE( grid2.isInitialized() );
}

TEST( quest_implicit_grid, implicit_grid_insert)
{
    SLIC_INFO("Test the constructor for a quest::ImplicitGrid object");

    using namespace implicit_grid_3D;

    typedef axom::primal::BoundingBox<int, DIM> RangeBox;

    // Note: A 10 x 10 x 10 implicit grid in the unit cube.
    //       Grid cells have a spacing of .1 along each dimension
    GridCell res = GridCell::make_point(10,10,10);
    BBox bbox(SpacePt::zero(), SpacePt::ones());
    int numElts = 10;

    GridT grid( bbox, res, numElts);

    /// Test insertion of several objects under the assumption that we are in a unit cube
    /// with a 10^3 grid resolution

    // Element with arbitrary bounding box within the grid
    BBox objBox1( SpacePt::make_point(.15, .25, .05), SpacePt::make_point(.45, .25, .35));
    grid.insert(objBox1, 1);
    RangeBox range1(GridCell::make_point(1,2,0), GridCell::make_point(4,2,3));

    // Some of the element's bounding box extends past the end of the grid
    BBox objBox2( SpacePt::make_point(.85, .85, .85), SpacePt::make_point(1.25, 1.25, 1.25));
    grid.insert(objBox2, 2);
    RangeBox range2(GridCell::make_point(8,8,8), GridCell::make_point(9,9,9));

    // Bounding box of element is completely outside grid
    BBox objBox3( SpacePt::make_point(1.85, 1.85, 1.85), SpacePt::make_point(2.25, 2.25, 2.25));
    grid.insert(objBox3, 3);

    // Bounding box of element overlaps with all grid elements
    BBox objBox4( SpacePt::make_point(-1.85, -1.85, -1.85), SpacePt::make_point(2.25, 2.25, 2.25));
    grid.insert(objBox4, 4);

    // Single point within a cell
    BBox objBox5( SpacePt::make_point(.25, .35, .15), SpacePt::make_point(.25, .35, .15));
    grid.insert(objBox5, 5);
    RangeBox range5(GridCell::make_point(2,3,1));

    // Single point on a grid cell boundary gets added to all incident grid cells
    BBox objBox6( SpacePt::make_point(.7, .7, .7));
    grid.insert(objBox6, 6);
    RangeBox range6(GridCell::make_point(6,6,6), GridCell::make_point(7,7,7));

    // Test that points are contained in the expected grid cells
    for(int i=0; i< res[0]; ++i)
    {
      for(int j=0; j< res[1]; ++j)
      {
        for(int k=0; k< res[2]; ++k)
        {
          GridCell gridCell = GridCell::make_point(i,j,k);

          ASSERT_FALSE( grid.contains(gridCell, 0));  // Not inserted
          ASSERT_FALSE( grid.contains(gridCell, 10)); // Not a valid element ID

          ASSERT_EQ( range1.contains(gridCell), grid.contains(gridCell, 1));
          ASSERT_EQ( range2.contains(gridCell), grid.contains(gridCell, 2));
          ASSERT_FALSE( grid.contains(gridCell, 3));
          ASSERT_TRUE ( grid.contains(gridCell, 4));
          ASSERT_EQ( range5.contains(gridCell), grid.contains(gridCell, 5));
          ASSERT_EQ( range6.contains(gridCell), grid.contains(gridCell, 6));
        }
      }
    }
}

TEST( quest_implicit_grid, implicit_grid_query)
{
    SLIC_INFO("Test the constructor for a quest::ImplicitGrid object");

    using namespace implicit_grid_3D;

    typedef axom::primal::BoundingBox<int, DIM> RangeBox;

    // Note: A 10 x 10 x 10 implicit grid in the unit cube.
    //       Grid cells have a spacing of .1 along each dimension
    GridCell res = GridCell::make_point(10,10,10);
    BBox bbox(SpacePt::zero(), SpacePt::ones());
    int numElts = 10;

    GridT grid( bbox, res, numElts);

    BBox objBox1( SpacePt::make_point(.15, .25, .05), SpacePt::make_point(.45, .25, .35));
    grid.insert(objBox1, 1);
    RangeBox range1(GridCell::make_point(1,2,0), GridCell::make_point(4,2,3));

    std::vector<int> candidatesMin = grid.getCandidatesAsArray( objBox1.getMin() );
    std::vector<int> candidatesMax = grid.getCandidatesAsArray( objBox1.getMax() );
    std::vector<int> candidatesMid = grid.getCandidatesAsArray( objBox1.centroid() );

    ASSERT_TRUE( std::find(candidatesMin.begin(), candidatesMin.end(), 1) != candidatesMin.end() );
    ASSERT_TRUE( std::find(candidatesMax.begin(), candidatesMax.end(), 1) != candidatesMax.end() );
    ASSERT_TRUE( std::find(candidatesMid.begin(), candidatesMid.end(), 1) != candidatesMid.end() );
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

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

#include "primal/Point.hpp"
#include "primal/BoundingBox.hpp"
#include "primal/RectangularLattice.hpp"

#include "slic/slic.hpp"
#include "slic/UnitTestLogger.hpp"

#include <cstdlib>
#include <limits>

// Define some helpful typedefs for 1D rectangular lattices
namespace lattice_1D {

const int DIM = 1;

typedef axom::primal::RectangularLattice< DIM > LatticeT;

typedef LatticeT::GridCell GridCell;
typedef LatticeT::SpacePoint SpacePt;
typedef LatticeT::SpaceVector SpaceVector;
typedef LatticeT::SpatialBoundingBox BBox;
typedef axom::primal::NumericArray< int,DIM > IntArray;
}

// Define some helpful typedefs for 2D rectangular lattices
namespace lattice_2D {

const int DIM = 2;

typedef axom::primal::RectangularLattice< DIM > LatticeT;

typedef LatticeT::GridCell GridCell;
typedef LatticeT::SpacePoint SpacePt;
typedef LatticeT::SpaceVector SpaceVector;
typedef LatticeT::SpatialBoundingBox BBox;
typedef axom::primal::NumericArray< int,DIM > IntArray;
}

// Define some helpful typedefs for 3D rectangular lattices
namespace lattice_3D {

const int DIM = 3;

typedef axom::primal::RectangularLattice< DIM > LatticeT;

typedef LatticeT::GridCell GridCell;
typedef LatticeT::SpacePoint SpacePt;
typedef LatticeT::SpaceVector SpaceVector;
typedef LatticeT::SpatialBoundingBox BBox;
typedef axom::primal::NumericArray< int,DIM > IntArray;
}

TEST( primal_rectangle_lattice, grid_ctor)
{
  SLIC_INFO("Testing lattice constructors in 1D, 2D and 3D");
  //1D
  {
    using namespace lattice_1D;
    SpacePt origin(1.1);
    SpaceVector spacing( SpacePt(.1) );

    LatticeT defaultGrid;
    ASSERT_EQ(SpacePt::zero(), defaultGrid.origin());
    ASSERT_EQ(SpaceVector(SpacePt(1)), defaultGrid.spacing());

    LatticeT grid(origin, spacing);
    ASSERT_EQ(origin, grid.origin());
    ASSERT_EQ(spacing, grid.spacing());
  }

  //2D
  {
    using namespace lattice_2D;
    SpacePt origin(1.1);
    SpaceVector spacing( SpacePt(.1) );

    LatticeT defaultGrid;
    ASSERT_EQ(SpacePt::zero(), defaultGrid.origin());
    ASSERT_EQ(SpaceVector(SpacePt(1)), defaultGrid.spacing());

    LatticeT grid(origin, spacing);
    ASSERT_EQ(origin, grid.origin());
    ASSERT_EQ(spacing, grid.spacing());
  }

  //3D
  {
    using namespace lattice_3D;
    SpacePt origin(1.1);
    SpaceVector spacing( SpacePt(.1) );

    LatticeT defaultGrid;
    ASSERT_EQ(SpacePt::zero(), defaultGrid.origin());
    ASSERT_EQ(SpaceVector(SpacePt(1)), defaultGrid.spacing());

    LatticeT grid(origin, spacing);
    ASSERT_EQ(origin, grid.origin());
    ASSERT_EQ(spacing, grid.spacing());
  }
}

TEST( primal_rectangle_lattice, from_bounding_box)
{
  SLIC_INFO("Testing lattice creation from 1D, 2D and 3D bounding boxes");

  //1D
  {
    using namespace lattice_1D;

    BBox bbox( SpacePt(1.25), SpacePt(2.5));
    IntArray res(5);
    LatticeT lattice = axom::primal::rectangular_lattice_from_bounding_box(bbox, res);

    EXPECT_DOUBLE_EQ(1.25, lattice.origin()[0]);
    EXPECT_DOUBLE_EQ( .25,  lattice.spacing()[0]);
  }

  //2D
  {
    using namespace lattice_2D;

    BBox bbox( SpacePt(1.25), SpacePt(2.5));
    IntArray res(5);
    LatticeT lattice = axom::primal::rectangular_lattice_from_bounding_box(bbox, res);

    EXPECT_DOUBLE_EQ(1.25, lattice.origin()[0]);
    EXPECT_DOUBLE_EQ( .25,  lattice.spacing()[0]);
  }

  //3D
  {
    using namespace lattice_3D;

    SpacePt bbMin = SpacePt::make_point(1.25, 2.5, 5.);
    SpacePt bbMax = SpacePt::make_point(2.5, 5., 10.);
    BBox bbox(bbMin, bbMax);

    int resData[3] = { 5, 50, 500 };
    IntArray res(resData);
    LatticeT lattice = axom::primal::rectangular_lattice_from_bounding_box(bbox, res);

    EXPECT_DOUBLE_EQ(1.25, lattice.origin()[0]);
    EXPECT_DOUBLE_EQ(2.5,  lattice.origin()[1]);
    EXPECT_DOUBLE_EQ(5.,   lattice.origin()[2]);

    EXPECT_DOUBLE_EQ( .25,  lattice.spacing()[0]);
    EXPECT_DOUBLE_EQ( .05,  lattice.spacing()[1]);
    EXPECT_DOUBLE_EQ( .01,  lattice.spacing()[2]);
  }
}

TEST( primal_rectangle_lattice, convert_point_cell_1D)
{
  SLIC_INFO("Testing point conversion of 1D lattice");

  using namespace lattice_1D;

  const double SPACING = 0.1;
  const double HALF_SPACING = SPACING / 2.;
  const double EPS = 1e-10;

  SpacePt origin = SpacePt::zero();
  SpaceVector spacing = SpaceVector( SpacePt( SPACING ) );

  LatticeT lattice(origin, spacing);

  ASSERT_EQ(origin, lattice.origin());
  ASSERT_EQ(spacing, lattice.spacing());

  // Test that we can map points to cells and cells to points
  for (int i=-10; i<=10; ++i) {
    // Point near the lower bounds of a cell
    // Note: We have to add an epsilon to guarantee we find the right cell
    SpacePt lowerPoint(SPACING * i + EPS);

    GridCell lowerCell = lattice.getGridCell( lowerPoint );
    EXPECT_EQ(i, lowerCell[0]);

    SpacePt ptFromLoweCell = lattice.getSpacePoint(lowerCell);
    EXPECT_DOUBLE_EQ(SPACING * i, ptFromLoweCell[0]);

    // Point near the middle of a cell
    SpacePt midPoint(SPACING * i + HALF_SPACING);

    GridCell midCell = lattice.getGridCell( midPoint );
    EXPECT_EQ(i, midCell[0]);

    SpacePt ptFromMidCell = lattice.getSpacePoint(midCell);
    EXPECT_DOUBLE_EQ(SPACING * i, ptFromMidCell[0]);

    // Point near the upper bounds of a cell
    SpacePt upperPoint(SPACING * (i+1) - EPS);
    GridCell upperCell = lattice.getGridCell( upperPoint );
    EXPECT_EQ(i, upperCell[0]);

    SpacePt ptFromUpperCell = lattice.getSpacePoint(upperCell);
    EXPECT_DOUBLE_EQ(SPACING * i, ptFromUpperCell[0]);

  }

}

TEST( primal_rectangle_lattice, convert_point_cell_2D)
{
  SLIC_INFO("Testing point conversion of 2D lattice");

  using namespace lattice_2D;

  SpacePt origin = SpacePt::make_point(1.1, 2.2);
  SpaceVector spacing = SpaceVector::make_vector(.1, .2);
  LatticeT lattice(origin, spacing);

  // Test a few hand-selected points in 2D

  {
    SpacePt pt = lattice.getSpacePoint( GridCell::zero() );
    EXPECT_DOUBLE_EQ(1.1, pt[0]);
    EXPECT_DOUBLE_EQ(2.2, pt[1]);
  }

  {
    SpacePt pt = lattice.getSpacePoint( GridCell::make_point(1,0));
    EXPECT_DOUBLE_EQ(1.2, pt[0]);
    EXPECT_DOUBLE_EQ(2.2, pt[1]);
  }
  {
    SpacePt pt = lattice.getSpacePoint( GridCell::make_point(0,1));
    EXPECT_DOUBLE_EQ(1.1, pt[0]);
    EXPECT_DOUBLE_EQ(2.4, pt[1]);
  }
  {
    SpacePt pt = lattice.getSpacePoint( GridCell::make_point(1,1));
    EXPECT_DOUBLE_EQ(1.2, pt[0]);
    EXPECT_DOUBLE_EQ(2.4, pt[1]);
  }

  {
    SpacePt pt = lattice.getSpacePoint( GridCell::make_point(-1,0));
    EXPECT_DOUBLE_EQ(1.0, pt[0]);
    EXPECT_DOUBLE_EQ(2.2, pt[1]);
  }
  {
    SpacePt pt = lattice.getSpacePoint( GridCell::make_point(0,-1));
    EXPECT_DOUBLE_EQ(1.1, pt[0]);
    EXPECT_DOUBLE_EQ(2.0, pt[1]);
  }
  {
    SpacePt pt = lattice.getSpacePoint( GridCell::make_point(-1,-1));
    EXPECT_DOUBLE_EQ(1.0, pt[0]);
    EXPECT_DOUBLE_EQ(2.0, pt[1]);
  }

}

TEST( primal_rectangle_lattice, negative_spacing_2D)
{
  SLIC_INFO("Testing point conversion of 3D lattice");

  using namespace lattice_2D;

  SpacePt origin = SpacePt::make_point(1.1, 2.2);
  SpaceVector spacing = SpaceVector::make_vector(-.1, -.2);
  LatticeT lattice(origin, spacing);

  // Test a few hand-selected points in 2D
  // Note: Spacing is negative

  {
    SpacePt pt = lattice.getSpacePoint( GridCell::zero() );
    EXPECT_DOUBLE_EQ(1.1, pt[0]);
    EXPECT_DOUBLE_EQ(2.2, pt[1]);
  }

  {
    SpacePt pt = lattice.getSpacePoint( GridCell::make_point(1,1));
    EXPECT_DOUBLE_EQ(1.0, pt[0]);
    EXPECT_DOUBLE_EQ(2.0, pt[1]);
  }

  {
    SpacePt pt = lattice.getSpacePoint( GridCell::make_point(-1,-1));
    EXPECT_DOUBLE_EQ(1.2, pt[0]);
    EXPECT_DOUBLE_EQ(2.4, pt[1]);
  }

}

TEST( primal_rectangle_lattice, cell_bounding_box_3D)
{
  SLIC_INFO("Testing cell bounding box in 3D");

  using namespace lattice_3D;

  SpacePt origin(1);
  SpaceVector spacing( SpacePt(2));
  LatticeT lattice(origin, spacing);

  // Get cell bounding box of cell at origin
  {
    SpacePt bbMin = origin;
    SpacePt bbMax = SpacePt(3);
    BBox expCellBBox(bbMin, bbMax);

    BBox cellBBox = lattice.getCellBounds( GridCell::zero() );

    for (int i=0; i < DIM; ++i) {
      EXPECT_DOUBLE_EQ(expCellBBox.getMin()[i], cellBBox.getMin()[i]);
      EXPECT_DOUBLE_EQ(expCellBBox.getMax()[i], cellBBox.getMax()[i]);
    }
  }

  // Get cell bounding box of cell at lattice point (5,4,3)
  {
    GridCell cell = GridCell::make_point(5,4,3);

    SpacePt bbMin = SpacePt::make_point(11., 9., 7.);
    SpacePt bbMax = SpacePt::make_point(13., 11., 9.);
    BBox expCellBBox(bbMin, bbMax);

    BBox cellBBox = lattice.getCellBounds( cell );

    for (int i=0; i < DIM; ++i) {
      EXPECT_DOUBLE_EQ(expCellBBox.getMin()[i], cellBBox.getMin()[i]);
      EXPECT_DOUBLE_EQ(expCellBBox.getMax()[i], cellBBox.getMax()[i]);
    }
  }
}

//----------------------------------------------------------------------

int main(int argc, char * argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  using axom::slic::UnitTestLogger;
  UnitTestLogger logger;  // create & initialize test logger,
  axom::slic::setLoggingMsgLevel( axom::slic::message::Info );

  // finalized when exiting main scope

  std::srand( 105 );

  result = RUN_ALL_TESTS();

  return result;
}

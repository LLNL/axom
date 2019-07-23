// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MIR_UTILITIES_TEST_H_
#define MIR_UTILITIES_TEST_H_

#include "gtest/gtest.h"

#include "axom/slic.hpp"
#include "axom/mir.hpp"

namespace mir = axom::mir;

//----------------------------------------------------------------------

TEST(mir_shape_tests, shape_dimesionality)
{
  EXPECT_EQ( mir::utilities::isShapeThreeDimensional(mir::Shape::Triangle), false );
  EXPECT_EQ( mir::utilities::isShapeThreeDimensional(mir::Shape::Quad), false );
  EXPECT_EQ( mir::utilities::isShapeThreeDimensional(mir::Shape::Tetrahedron), true );
  EXPECT_EQ( mir::utilities::isShapeThreeDimensional(mir::Shape::Pyramid), true );
  EXPECT_EQ( mir::utilities::isShapeThreeDimensional(mir::Shape::Triangular_Prism), true );
  EXPECT_EQ( mir::utilities::isShapeThreeDimensional(mir::Shape::Hexahedron), true );
}

//----------------------------------------------------------------------

TEST(mir_shape_tests, check_shape_central_vertex)
{
  EXPECT_EQ( mir::utilities::isCenterVertex(mir::Shape::Triangle, 6), false);

  EXPECT_EQ( mir::utilities::isCenterVertex(mir::Shape::Quad, 8), false);

  EXPECT_EQ( mir::utilities::isCenterVertex(mir::Shape::Tetrahedron, 9), false);
  EXPECT_EQ( mir::utilities::isCenterVertex(mir::Shape::Tetrahedron, 10), true);
  EXPECT_EQ( mir::utilities::isCenterVertex(mir::Shape::Tetrahedron, 11), false);

  EXPECT_EQ( mir::utilities::isCenterVertex(mir::Shape::Pyramid, 12), false);
  EXPECT_EQ( mir::utilities::isCenterVertex(mir::Shape::Pyramid, 13), true);
  EXPECT_EQ( mir::utilities::isCenterVertex(mir::Shape::Pyramid, 14), false);
  
  EXPECT_EQ( mir::utilities::isCenterVertex(mir::Shape::Triangular_Prism, 14), false);
  EXPECT_EQ( mir::utilities::isCenterVertex(mir::Shape::Triangular_Prism, 15), true);
  EXPECT_EQ( mir::utilities::isCenterVertex(mir::Shape::Triangular_Prism, 16), false);
  
  EXPECT_EQ( mir::utilities::isCenterVertex(mir::Shape::Hexahedron, 19), false);
  EXPECT_EQ( mir::utilities::isCenterVertex(mir::Shape::Hexahedron, 20), true);
  EXPECT_EQ( mir::utilities::isCenterVertex(mir::Shape::Hexahedron, 21), false);
}

//----------------------------------------------------------------------

TEST(mir_shape_tests, determine_central_vertex)
{
  EXPECT_EQ(  mir::utilities::getCenterVertex(mir::Shape::Triangle), -1 );
  EXPECT_EQ(  mir::utilities::getCenterVertex(mir::Shape::Quad), -1 );

  EXPECT_EQ(  mir::utilities::getCenterVertex(mir::Shape::Tetrahedron), 10 );
  EXPECT_EQ(  mir::utilities::getCenterVertex(mir::Shape::Pyramid), 13 );
  EXPECT_EQ(  mir::utilities::getCenterVertex(mir::Shape::Triangular_Prism), 15 );
  EXPECT_EQ(  mir::utilities::getCenterVertex(mir::Shape::Hexahedron), 20 );
}

//----------------------------------------------------------------------

TEST(mir_compute_averages, float_value)
{
  std::vector<axom::float64> values = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0};

  axom::float64 average = mir::utilities::computeAverageFloat( values );

  EXPECT_DOUBLE_EQ(  average, 2.5  );
}

//----------------------------------------------------------------------

TEST(mir_compute_averages, point_value)
{
  std::vector<mir::Point2> points = {
    mir::Point2::make_point(0.0, 0.0, 0.0),
    mir::Point2::make_point(1.0, 0.0, 0.0),
    mir::Point2::make_point(0.0, 1.0, 0.0),
    mir::Point2::make_point(0.0, 0.0, 1.0),
    mir::Point2::make_point(1.0, 1.0, 0.0),
    mir::Point2::make_point(1.0, 0.0, 1.0),
    mir::Point2::make_point(0.0, 1.0, 1.0),
    mir::Point2::make_point(1.0, 1.0, 1.0)
  };

  mir::Point2 centroid = mir::utilities::computeAveragePoint( points );

  EXPECT_DOUBLE_EQ(  centroid[0], 0.5  );
  EXPECT_DOUBLE_EQ(  centroid[1], 0.5  );
  EXPECT_DOUBLE_EQ(  centroid[2], 0.5  );
}

//----------------------------------------------------------------------

TEST(mir_compute_shape_volumes, triangle)
{
  std::vector<mir::Point2> points = {
    mir::Point2::make_point(0.0, 0.0, 0.0),
    mir::Point2::make_point(1.0, 0.0, 0.0),
    mir::Point2::make_point(0.0, 1.0, 0.0)
  };
  
  axom::float64 area = mir::utilities::computeShapeVolume( mir::Shape::Triangle, points.data() );

  EXPECT_DOUBLE_EQ(  area, 0.5  );
}

//----------------------------------------------------------------------

TEST(mir_compute_shape_volumes, quad)
{
  std::vector<mir::Point2> points = {
    mir::Point2::make_point(0.0, 0.0, 0.0),
    mir::Point2::make_point(1.0, 0.0, 0.0),
    mir::Point2::make_point(0.0, 1.0, 0.0),
    mir::Point2::make_point(1.0, 1.0, 0.0),
  };
  
  axom::float64 area = mir::utilities::computeShapeVolume( mir::Shape::Quad, points.data() );

  EXPECT_DOUBLE_EQ(  area, 1.0  );
}

//----------------------------------------------------------------------

TEST(mir_compute_shape_volumes, tetrahedron)
{
  std::vector<mir::Point2> points = {
    mir::Point2::make_point(0.0, 0.0, 0.0),
    mir::Point2::make_point(0.0, 0.0, 1.0),
    mir::Point2::make_point(1.0, 0.0, 0.0),
    mir::Point2::make_point(1.0, 1.0, 0.0),
  };
  
  axom::float64 volume = mir::utilities::computeShapeVolume( mir::Shape::Tetrahedron, points.data() );

  EXPECT_DOUBLE_EQ(  volume, 1.0/ (6.0)  );
}

//----------------------------------------------------------------------

TEST(mir_compute_shape_volumes, pyramid)
{
  std::vector<mir::Point2> points = {
    mir::Point2::make_point(0.0, 0.0, 1.0),
    mir::Point2::make_point(1.0, 0.0, 1.0),
    mir::Point2::make_point(1.0, 0.0, 0.0),
    mir::Point2::make_point(0.0, 0.0, 0.0),
    mir::Point2::make_point(0.5, 0.5, 0.5)
  };
  
  axom::float64 volume = mir::utilities::computeShapeVolume( mir::Shape::Pyramid, points.data() );

  EXPECT_DOUBLE_EQ(  volume, 1.0/ (6.0)  );
}

//----------------------------------------------------------------------

TEST(mir_compute_shape_volumes, hexahedron)
{
  std::vector<mir::Point2> points = {
    mir::Point2::make_point(0.0, 0.0, 1.0),
    mir::Point2::make_point(1.0, 0.0, 1.0),
    mir::Point2::make_point(1.0, 0.0, 0.0),
    mir::Point2::make_point(0.0, 0.0, 0.0),

    mir::Point2::make_point(0.0, 1.0, 1.0),
    mir::Point2::make_point(1.0, 1.0, 1.0),
    mir::Point2::make_point(1.0, 1.0, 0.0),
    mir::Point2::make_point(0.0, 1.0, 0.0)
  };
  
  axom::float64 volume = mir::utilities::computeShapeVolume( mir::Shape::Hexahedron, points.data() );

  EXPECT_DOUBLE_EQ(  volume, 1.0  );
}

//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::UnitTestLogger logger;  // create & initialize test logger,

  result = RUN_ALL_TESTS();
  return result;
}


#endif //  MIR_UTILITIES_TEST_H_

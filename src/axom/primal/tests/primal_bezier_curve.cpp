// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/* /file bezier_test.cpp
 * /brief This file tests the BezierCurve.hpp and eval_bezier.hpp files
 */

#include "gtest/gtest.h"

#include "axom/primal/geometry/BezierCurve.hpp"

// _using_start
using namespace axom;
using namespace primal;
// _using_end
//----------------------------------------------------------------------------------
TEST( primal_beziercurve, bezier_constructor )
{
  static const int DIM = 3;
  typedef double CoordType;
  typedef primal::BezierCurve< CoordType, DIM > BezierCurveType;
  SLIC_INFO("primal: testing default bezier constructor ");
  BezierCurveType bCurve;

  EXPECT_EQ(bCurve.getControlPoints(),
            (std::vector< Point< CoordType, DIM > >()));
}

//----------------------------------------------------------------------------------
TEST( primal_beziercurve, bezier_order_constructor )
{
  static const int DIM = 3;
  typedef double CoordType;
  typedef primal::BezierCurve< CoordType, DIM > BezierCurveType;
  SLIC_INFO("primal: testing bezier order constructor ");

  BezierCurveType bCurve(1);

  EXPECT_EQ(bCurve.getOrder(), 1);
}

//----------------------------------------------------------------------------------
TEST( primal_beziercurve, bezier_add_controlpoints )
{
  static const int DIM =3;
  typedef double CoordType;
  typedef primal::Point< CoordType, DIM > PointType;
  typedef primal::BezierCurve< CoordType, DIM > BezierCurveType;
  SLIC_INFO("primal: testing adding control points to empty beziercurve");

  BezierCurveType bCurve;

  CoordType coords[6] = {.6,1.2,1.0,0.0,1.6,1.8};
  bCurve.addControlpoint(
    PointType::make_point( coords[0], coords[1], coords[2] ) );
  bCurve.addControlpoint(
    PointType::make_point( coords[3], coords[4], coords[5] ) );

  EXPECT_EQ(bCurve.getOrder(), 1);
  for (int i=0 ; i<DIM ; i++)
  {
    for (int j=0 ; j<2 ; j++)
    {
      EXPECT_TRUE(utilities::isNearlyEqual(bCurve[j][i],coords[j*DIM+i],
                                           1.0e-14));
    }
  }
}

//----------------------------------------------------------------------------------
TEST( primal_beziercurve, bezier_point_array_constructor )
{
  static const int DIM =3;
  typedef double CoordType;
  typedef primal::Point< CoordType, DIM > PointType;
  typedef primal::BezierCurve< CoordType, DIM > BezierCurveType;
  SLIC_INFO("\nprimal: testing point array constructor");

  CoordType coords[6] = {.6,1.2,1.0,0.0,1.6,1.8};
  PointType pts[2] = {
    PointType::make_point(coords[0], coords[1], coords[2]),
    PointType::make_point(coords[3], coords[4], coords[5])
  };

  BezierCurveType bCurve(pts,1);
  EXPECT_EQ(bCurve.getOrder(), 1);
  for (int i=0 ; i<DIM ; i++)
  {
    for (int j=0 ; j<2 ; j++)
    {
      EXPECT_TRUE(utilities::isNearlyEqual(bCurve[j][i],coords[j*DIM+i],
                                           1.0e-14));
    }
  }
}

//----------------------------------------------------------------------------------
TEST( primal_beziercurve, bezier_coordinate_array_constructor )
{
  static const int DIM =3;
  typedef double CoordType;
  typedef primal::BezierCurve< CoordType, DIM > BezierCurveType;
  SLIC_INFO("\nprimal: testing coordinate array constructor");

  CoordType coords[6] = {.6,1.2,1.0,0.0,1.6,1.8};
  CoordType coordsToPass[6] = {coords[0], coords[3], coords[1],
                               coords[4], coords[2], coords[5]};

  BezierCurveType bCurve(coordsToPass,1);
  EXPECT_EQ(bCurve.getOrder(), 1);
  for (int i=0 ; i<DIM ; i++)
  {
    for (int j=0 ; j<2 ; j++)
    {
      EXPECT_TRUE(utilities::isNearlyEqual(bCurve[j][i],coords[j*DIM+i],
                                           1.0e-14));
    }
  }
}

//----------------------------------------------------------------------------------
TEST( primal_beziercurve, bezier_evaluation )
{
  static const int DIM =3;
  typedef double CoordType;
  typedef primal::Point< CoordType, DIM > PointType;
  typedef primal::BezierCurve< CoordType, DIM > BezierCurveType;
  SLIC_INFO("\nprimal: testing bezier evaluation");

  const int nbr_points = 4;
  PointType data[nbr_points];
  data[0] = PointType::make_point(0.6, 1.2, 1.0);
  data[1] = PointType::make_point(1.3, 1.6, 1.8);
  data[2] = PointType::make_point(2.9, 2.4, 2.3);
  data[3] = PointType::make_point(3.2, 3.5, 3.0);
  BezierCurveType b2Curve(data, nbr_points-1);

  PointType midtval = PointType::make_point(2.05,2.0875,2.0375);
  for ( int i=0 ; i<DIM ; i++)
  {
    EXPECT_TRUE(utilities::isNearlyEqual(b2Curve[0][i],
                                         b2Curve.eval_bezier(0.0)[i]));
    EXPECT_TRUE(utilities::isNearlyEqual(b2Curve[nbr_points-1][i],
                                         b2Curve.eval_bezier(1.0)[i]));
    EXPECT_TRUE(utilities::isNearlyEqual((midtval[i]),
                                         b2Curve.eval_bezier(0.5)[i]));
  }
}

//----------------------------------------------------------------------------------
TEST( primal_beziercurve, bezier_split )
{
  static const int DIM =3;
  typedef double CoordType;
  typedef primal::Point< CoordType, DIM > PointType;
  typedef primal::BezierCurve< CoordType, DIM > BezierCurveType;
  SLIC_INFO("\nprimal: testing bezier splitting");

  const int nbr_points = 4;
  PointType data[nbr_points];
  data[0] = PointType::make_point(0.6, 1.2, 1.0);
  data[1] = PointType::make_point(1.3, 1.6, 1.8);
  data[2] = PointType::make_point(2.9, 2.4, 2.3);
  data[3] = PointType::make_point(3.2, 3.5, 3.0);
  BezierCurveType b2Curve(data, nbr_points-1);

  BezierCurveType b3Curve(3); // Checks split_bezier with order constructor
  BezierCurveType b4Curve; // Checks split_bezier with default constructor
  b2Curve.split_bezier(.5,b3Curve,b4Curve);

  CoordType b3Coords[12] = {.6, .95, 1.525, 2.05,
                            1.2, 1.4, 1.7, 2.0875,
                            1.0, 1.4, 1.725, 2.0375};
  CoordType b4Coords[12] = {2.05, 2.575, 3.05, 3.2,
                            2.0875, 2.475, 2.95, 3.5,
                            2.0375, 2.35, 2.65, 3.0};
  BezierCurveType b3True(b3Coords,3);
  BezierCurveType b4True(b4Coords,3);
  for ( int i=0 ; i<DIM ; i++)
  {
    for ( int j=0 ; j<(nbr_points-1) ; j++)
    {
      EXPECT_DOUBLE_EQ(b3Curve[j][i], b3True[j][i]);
      EXPECT_DOUBLE_EQ(b4Curve[j][i], b4True[j][i]);
    }
  }

}

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::UnitTestLogger logger;     // create & initialize test logger,

  result = RUN_ALL_TESTS();

  return result;
}

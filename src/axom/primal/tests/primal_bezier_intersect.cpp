// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/* /file bezier_test.cpp
 * /brief This file tests the BezierCurve.hpp and eval_bezier.hpp files
 */

#include "gtest/gtest.h"

#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/operators/intersect_bezier.hpp"

// _using_start
using namespace axom;
using namespace primal;
// _using_end
//----------------------------------------------------------------------------------
TEST( primal_bezier_inter, bezier_inters )
{
  static const int DIM =2;
  typedef double CoordType;
  typedef primal::Point< CoordType, DIM > PointType;
  typedef primal::BezierCurve< CoordType, DIM > BezierCurveType;
  SLIC_INFO("\nprimal: testing bezier intersection");

  const int nbr_points = 4;
  PointType data[nbr_points];
  data[0] = PointType::make_point(0.0, 0.0);
  data[1] = PointType::make_point(1.0, 0.0);
  data[2] = PointType::make_point(2.0, 0.0);
  data[3] = PointType::make_point(3.0, 0.0);
  BezierCurveType b2Curve(data, nbr_points-1);

  data[0] = PointType::make_point(0.0, 0.5);
  data[1] = PointType::make_point(1.0,-1.0);
  data[2] = PointType::make_point(2.0, 1.0);
  data[3] = PointType::make_point(3.0,-0.5);
  BezierCurveType b5Curve(data, nbr_points-1);

  CoordType actualintersects[3];
  actualintersects[0] = 0.17267316464601146;
  actualintersects[1] = 0.5;
  actualintersects[2] = 0.827326835353989;

  std::vector<CoordType> t1;
  std::vector<CoordType> s1;
  intersect_bezier(b2Curve, b5Curve, s1, t1);
  for (int i=0 ; i<static_cast<int>(s1.size()) ; i++)
  {
    EXPECT_DOUBLE_EQ(s1[i], actualintersects[i]);
    EXPECT_DOUBLE_EQ(t1[i], actualintersects[i]);
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

// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"

#include "gtest/gtest.h"

namespace primal = axom::primal;

TEST(primal_quadrilateral, area_2D)
{
  constexpr int DIM = 2;
  constexpr double EPS = 1e-12;
  using CoordinateType = double;
  using PointType = primal::Point<CoordinateType, DIM>;
  using QuadrilateralType = primal::Quadrilateral<CoordinateType, DIM>;

  // This is a concave quadrilateral
  QuadrilateralType quad{PointType{-1.0, 0.1},
                         PointType{ 0.0, 1.0},
                         PointType{ 2.0, 0.0},
                         PointType{-0.1, 0.5}};

  const CoordinateType signedArea = quad.signedArea();
  const CoordinateType area = quad.area();

  EXPECT_NEAR(-0.755, signedArea, EPS);
  EXPECT_NEAR(0.755, area, EPS);
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger(axom::slic::message::Info);

  int result = RUN_ALL_TESTS();
  return result;
}

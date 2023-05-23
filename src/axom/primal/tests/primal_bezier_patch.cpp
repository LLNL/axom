// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*! 
 * \file primal_bezier_patch.cpp
 * \brief This file tests primal's Bezier patch functionality
 */

#include "gtest/gtest.h"

#include "axom/slic.hpp"

#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/geometry/BezierPatch.hpp"
#include "axom/primal/operators/squared_distance.hpp"

namespace primal = axom::primal;

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, initial_test)
{
  const int DIM = 3;
  using CoordType = double;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType, DIM>;
  using CoordsVec = BezierPatchType::CoordsVec;
  using CoordsMat = BezierPatchType::CoordsMat;
  using PointType = primal::Point<CoordType, DIM>;

  /* clang-format off */
  PointType PatchNodes[4 * 3] = {PointType {0, 0, 0}, PointType {0, 4,  0}, PointType {0, 8, -3},
                                 PointType {2, 0, 6}, PointType {2, 4,  0}, PointType {2, 8,  0}, 
                                 PointType {4, 0, 0}, PointType {4, 4,  0}, PointType {4, 8,  3},
                                 PointType {6, 0, 0}, PointType {6, 4, -3}, PointType {6, 8,  0}};
  
  double PatchWeights[4 * 3] = {1.0, 2.0, 3.0,
                                2.0, 3.0, 4.0,
                                3.0, 4.0, 5.0,
                                4.0, 5.0, 6.0};
  /* clang-format on */

  BezierPatchType bPatch(PatchNodes, PatchWeights, 3, 2);
  BezierPatchType bPatch1(3, 2), bPatch2(3, 2), bPatch3(3, 2), bPatch4(3, 2);

  double t = 0.25, s = 0.75;

  bPatch.split(t, s, bPatch1, bPatch2, bPatch3, bPatch4);

  std::cout << bPatch.evaluate(t, s) << std::endl;
  std::cout << bPatch1 << std::endl;
  std::cout << bPatch2 << std::endl;
  std::cout << bPatch3 << std::endl;
  std::cout << bPatch4 << std::endl;
}

//------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}

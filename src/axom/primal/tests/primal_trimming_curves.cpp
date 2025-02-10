// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*! 
 * \file primal_nurbs_patch.cpp
 * \brief This file tests primal's NURBS patch functionality
 */

#include "gtest/gtest.h"

#include "axom/slic.hpp"

#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/geometry/BezierPatch.hpp"
#include "axom/primal/geometry/NURBSPatch.hpp"
#include "axom/primal/operators/squared_distance.hpp"

#include "axom/core/numerics/matvecops.hpp"

namespace primal = axom::primal;

//------------------------------------------------------------------------------
TEST(primal_nurbspatch, basic_operations)
{
  SLIC_INFO("Testing access of trimming curves");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using ParameterPointType = primal::Point<CoordType, 2>;
  using NURBSPatchType = primal::NURBSPatch<CoordType, DIM>;

  const int degree_u = 1;
  const int degree_v = 1;

  const int npts_u = 3;
  const int npts_v = 3;

  // Construct from 1D axom::Array
  axom::Array<PointType> controlPointsArray({PointType {0.0, 0.0, 1.0},
                                             PointType {0.0, 1.0, 0.0},
                                             PointType {0.0, 2.0, 0.0},
                                             PointType {1.0, 0.0, 0.0},
                                             PointType {1.0, 1.0, -1.0},
                                             PointType {1.0, 2.0, 0.0},
                                             PointType {2.0, 0.0, 0.0},
                                             PointType {2.0, 1.0, 0.0},
                                             PointType {2.0, 2.0, 1.0}});
  axom::Array<CoordType> weightsArray(
    {1.0, 2.0, 3.0, 2.0, 3.0, 4.0, 3.0, 4.0, 5.0});

  NURBSPatchType nPatch(controlPointsArray, npts_u, npts_v, degree_u, degree_v);

  EXPECT_FALSE(nPatch.isTrimmed());
  EXPECT_EQ(nPatch.getNumTrimmingCurves(), 0);

  // Add simple trimming curves
  nPatch.makeSimpleTrimmed();
  EXPECT_TRUE(nPatch.isTrimmed());
  EXPECT_EQ(nPatch.getNumTrimmingCurves(), 4);

  // Remove all trimming curves
  nPatch.makeUntrimmed();
  EXPECT_FALSE(nPatch.isTrimmed());

  // Add a simple trimming curve from a NURBS curve
  axom::Array<ParameterPointType> controlPointsArray {
    ParameterPointType {0.25, 0.25},
    ParameterPointType {0.75, 0.25},
    ParameterPointType {0.75, 0.75},
    ParameterPointType {0.25, 0.75},
    ParameterPointType {0.25, 0.25}};
}

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}

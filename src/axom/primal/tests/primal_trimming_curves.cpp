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

#include "axom/primal/geometry/NURBSPatch.hpp"

namespace primal = axom::primal;

//------------------------------------------------------------------------------
TEST(primal_nurbspatch_trimming, basic_operations)
{
  SLIC_INFO("Testing access of trimming curves");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSPatchType = primal::NURBSPatch<CoordType, DIM>;

  using ParameterPointType = primal::Point<CoordType, 2>;
  using TrimmingCurveType = primal::NURBSCurve<CoordType, 2>;

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
  axom::Array<ParameterPointType> trimmingCurveControlPoints {
    ParameterPointType {0.25, 0.25},
    ParameterPointType {0.75, 0.25},
    ParameterPointType {0.75, 0.75},
    ParameterPointType {0.25, 0.75},
    ParameterPointType {0.25, 0.25}};

  TrimmingCurveType trimmingCurve(trimmingCurveControlPoints, 2);

  nPatch.addTrimmingCurve(trimmingCurve);
  EXPECT_TRUE(nPatch.isTrimmed());
  EXPECT_EQ(nPatch.getNumTrimmingCurves(), 1);

  // Copy the patch
  NURBSPatchType nPatchCopy(nPatch);

  // Check that the copy is the same
  EXPECT_TRUE(nPatchCopy.isTrimmed());
  EXPECT_EQ(nPatchCopy.getNumTrimmingCurves(), 1);
  EXPECT_TRUE(nPatchCopy == nPatch);

  // Delete all trimming curves from the copy
  nPatchCopy.makeUntrimmed();
  EXPECT_FALSE(nPatchCopy.isTrimmed());

  // Check that the two patches are not equal
  EXPECT_FALSE(nPatchCopy == nPatch);
}

//------------------------------------------------------------------------------
TEST(primal_nurbspatch_trimming, visibility_queries)
{
  SLIC_INFO("Testing queries for trimming curves");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSPatchType = primal::NURBSPatch<CoordType, DIM>;

  using ParameterPointType = primal::Point<CoordType, 2>;
  using TrimmingCurveType = primal::NURBSCurve<CoordType, 2>;

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

  NURBSPatchType nPatch(controlPointsArray, npts_u, npts_v, degree_u, degree_v);

  // Add a simple trimming curve from a NURBS curve
  axom::Array<ParameterPointType> trimmingCurveControlPoints {
    ParameterPointType {0.25, 0.25},
    ParameterPointType {0.75, 0.25},
    ParameterPointType {0.75, 0.75},
    ParameterPointType {0.25, 0.75},
    ParameterPointType {0.25, 0.25}};

  TrimmingCurveType trimmingCurve(trimmingCurveControlPoints, 2);
  nPatch.addTrimmingCurve(trimmingCurve);

  // Check that a point on the patch is visible
  EXPECT_TRUE(nPatch.isVisible(ParameterPointType {0.5, 0.5}));
  EXPECT_TRUE(nPatch.isVisible(ParameterPointType {0.4, 0.5}));
  EXPECT_TRUE(nPatch.isVisible(ParameterPointType {0.4, 0.6}));

  // Check that a point outside the patch is invisible
  EXPECT_FALSE(nPatch.isVisible(ParameterPointType {0.0, 0.0}));
  EXPECT_FALSE(nPatch.isVisible(ParameterPointType {0.0, 0.1}));
  EXPECT_FALSE(nPatch.isVisible(ParameterPointType {0.9, 0.0}));
}

//------------------------------------------------------------------------------
TEST(primal_nurbspatch_trimming, normal_evaluation)
{
  SLIC_INFO("Testing queries for trimming curves");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using NURBSPatchType = primal::NURBSPatch<CoordType, DIM>;

  using ParameterPointType = primal::Point<CoordType, 2>;
  using TrimmingCurveType = primal::NURBSCurve<CoordType, 2>;

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

  NURBSPatchType nPatch(controlPointsArray, npts_u, npts_v, degree_u, degree_v);

  // Add a simple trimming curve from a NURBS curve
  axom::Array<ParameterPointType> trimmingCurveControlPoints {
    ParameterPointType {0.25, 0.25},
    ParameterPointType {0.75, 0.25},
    ParameterPointType {0.75, 0.75},
    ParameterPointType {0.25, 0.75},
    ParameterPointType {0.25, 0.25}};

  TrimmingCurveType trimmingCurve(trimmingCurveControlPoints, 2);
  nPatch.addTrimmingCurve(trimmingCurve);

  // Copy the patch, and flip the normals of the copy
  NURBSPatchType nPatchCopy(nPatch);
//   nPatchCopy.flipNormals();

//   // Check that the normals are flipped
//   constexpr int npts = 11;
//   double u_pts[npts], v_pts[npts];
//   axom::numerics::linspace(0.0, 1.0, u_pts, npts);
//   axom::numerics::linspace(0.0, 1.0, v_pts, npts);

//   for(auto u : u_pts)
//   {
//     for(auto v : v_pts)
//     {
//       auto normal = nPatch.normal(u, v);
//       auto normalCopy = nPatchCopy.normal(u, v);

//       for(int i = 0; i < DIM; ++i)
//       {
//         EXPECT_DOUBLE_EQ(-normal[i], normalCopy[i]);
//       }
//     }
//   }
}

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}

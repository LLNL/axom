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
TEST(primal_bezierpatch, constructor)
{
  using CoordType = double;
  using BezierPatchType = primal::BezierPatch<CoordType>;
  using CoordsMat = BezierPatchType::CoordsMat;

  {
    SLIC_INFO("Testing default BezierPatch constructor ");
    BezierPatchType bPatch;
    EXPECT_FALSE(bPatch.isRational());

    int expOrder_u = -1, expOrder_v = -1;
    EXPECT_EQ(expOrder_u, bPatch.getOrder_u());
    EXPECT_EQ(expOrder_v, bPatch.getOrder_v());

    EXPECT_EQ(expOrder_u + 1, bPatch.getControlPoints().shape()[0]);
    EXPECT_EQ(expOrder_v + 1, bPatch.getControlPoints().shape()[1]);

    EXPECT_EQ(CoordsMat(), bPatch.getControlPoints());
  }

  {
    SLIC_INFO("Testing BezierPatch order constructor ");
    BezierPatchType bPatch(1, 1);
    EXPECT_FALSE(bPatch.isRational());

    int expOrder_u = 1, expOrder_v = 1;
    EXPECT_EQ(expOrder_u, bPatch.getOrder_u());
    EXPECT_EQ(expOrder_v, bPatch.getOrder_v());

    EXPECT_EQ(expOrder_u + 1, bPatch.getControlPoints().shape()[0]);
    EXPECT_EQ(expOrder_v + 1, bPatch.getControlPoints().shape()[1]);
  }
}

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, set_order)
{
  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType>;

  SLIC_INFO("Test adding control points to an empty Bezier patch");

  BezierPatchType bPatch;
  EXPECT_EQ(bPatch.getOrder_u(), -1);
  EXPECT_EQ(bPatch.getOrder_v(), -1);

  const int order_u = 1;
  const int order_v = 1;

  bPatch.setOrder(order_u, order_v);
  EXPECT_EQ(bPatch.getOrder_u(), order_u);
  EXPECT_EQ(bPatch.getOrder_v(), order_v);

  PointType controlPoints[4] = {PointType {0.0, 0.0, 1.0},
                                PointType {0.0, 1.0, 0.0},
                                PointType {1.0, 0.0, 0.0},
                                PointType {1.0, 1.0, -1.0}};

  bPatch(0, 0) = controlPoints[0];
  bPatch(0, 1) = controlPoints[1];
  bPatch(1, 0) = controlPoints[2];
  bPatch(1, 1) = controlPoints[3];

  for(int p = 0; p <= bPatch.getOrder_u(); ++p)
    for(int q = 0; q <= bPatch.getOrder_v(); ++q)
    {
      auto& pt = bPatch(p, q);
      for(int i = 0; i < DIM; ++i)
        EXPECT_DOUBLE_EQ(controlPoints[p * (order_u + 1) + q][i], pt[i]);
    }
}

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, array_constructors)
{
  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType>;

  SLIC_INFO("Testing point array constructor");

  const int order_u = 1;
  const int order_v = 1;

  PointType controlPoints[4] = {PointType {0.0, 0.0, 1.0},
                                PointType {0.0, 1.0, 0.0},
                                PointType {1.0, 0.0, 0.0},
                                PointType {1.0, 1.0, -1.0}};

  CoordType weights[4] = {1.0, 2.0, 1.0, 0.5};

  BezierPatchType nonrational_patch(controlPoints, order_u, order_v);
  EXPECT_EQ(nonrational_patch.getOrder_u(), order_u);
  EXPECT_EQ(nonrational_patch.getOrder_v(), order_v);
  EXPECT_FALSE(nonrational_patch.isRational());

  for(int p = 0; p <= nonrational_patch.getOrder_u(); ++p)
    for(int q = 0; q <= nonrational_patch.getOrder_v(); ++q)
    {
      auto& pt = nonrational_patch(p, q);
      for(int i = 0; i < DIM; ++i)
        EXPECT_DOUBLE_EQ(controlPoints[p * (order_u + 1) + q][i], pt[i]);
    }

  BezierPatchType nonrational_patch_again(controlPoints, nullptr, order_u, order_v);
  EXPECT_EQ(nonrational_patch_again.getOrder_u(), order_u);
  EXPECT_EQ(nonrational_patch_again.getOrder_v(), order_v);
  EXPECT_FALSE(nonrational_patch_again.isRational());

  for(int p = 0; p <= nonrational_patch_again.getOrder_u(); ++p)
    for(int q = 0; q <= nonrational_patch_again.getOrder_v(); ++q)
    {
      auto& pt = nonrational_patch_again(p, q);
      for(int i = 0; i < DIM; ++i)
        EXPECT_DOUBLE_EQ(controlPoints[p * (order_u + 1) + q][i], pt[i]);
    }

  BezierPatchType rational_patch(controlPoints, weights, order_u, order_v);
  EXPECT_EQ(rational_patch.getOrder_u(), order_u);
  EXPECT_EQ(rational_patch.getOrder_v(), order_v);
  EXPECT_TRUE(rational_patch.isRational());

  for(int p = 0; p <= rational_patch.getOrder_u(); ++p)
    for(int q = 0; q <= rational_patch.getOrder_v(); ++q)
    {
      auto& pt = rational_patch(p, q);
      for(int i = 0; i < DIM; ++i)
        EXPECT_DOUBLE_EQ(controlPoints[p * (order_u + 1) + q][i], pt[i]);
      EXPECT_DOUBLE_EQ(weights[p * (order_u + 1) + q],
                       rational_patch.getWeight(p, q));
    }
}

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, axom_array_constructors)
{
  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType>;

  SLIC_INFO("Testing point array constructor");

  const int order_u = 1;
  const int order_v = 1;

  // Construct with 1D axom arrays
  axom::Array<PointType> controlPoints({PointType {0.0, 0.0, 1.0},
                                        PointType {0.0, 1.0, 0.0},
                                        PointType {1.0, 0.0, 0.0},
                                        PointType {1.0, 1.0, -1.0}});

  axom::Array<CoordType> weights({1.0, 2.0, 1.0, 0.5});

  BezierPatchType nonrational_patch(controlPoints, order_u, order_v);
  EXPECT_EQ(nonrational_patch.getOrder_u(), order_u);
  EXPECT_EQ(nonrational_patch.getOrder_v(), order_v);
  EXPECT_FALSE(nonrational_patch.isRational());

  for(int p = 0; p <= nonrational_patch.getOrder_u(); ++p)
    for(int q = 0; q <= nonrational_patch.getOrder_v(); ++q)
    {
      auto& pt = nonrational_patch(p, q);
      for(int i = 0; i < DIM; ++i)
        EXPECT_DOUBLE_EQ(controlPoints[p * (order_u + 1) + q][i], pt[i]);
    }

  BezierPatchType rational_patch(controlPoints, weights, order_u, order_v);
  EXPECT_EQ(rational_patch.getOrder_u(), order_u);
  EXPECT_EQ(rational_patch.getOrder_v(), order_v);
  EXPECT_TRUE(rational_patch.isRational());

  for(int p = 0; p <= rational_patch.getOrder_u(); ++p)
    for(int q = 0; q <= rational_patch.getOrder_v(); ++q)
    {
      auto& pt = rational_patch(p, q);
      for(int i = 0; i < DIM; ++i)
        EXPECT_DOUBLE_EQ(controlPoints[p * (order_u + 1) + q][i], pt[i]);
      EXPECT_DOUBLE_EQ(weights[p * (order_u + 1) + q],
                       rational_patch.getWeight(p, q));
    }

  // Construct with 2D axom arrays
  axom::Array<PointType, 2> controlPoints_2D(2, 2);
  controlPoints_2D(0, 0) = PointType {0.0, 0.0, 1.0};
  controlPoints_2D(0, 1) = PointType {0.0, 1.0, 0.0};
  controlPoints_2D(1, 0) = PointType {1.0, 0.0, 0.0};
  controlPoints_2D(1, 1) = PointType {1.0, 1.0, -1.0};

  axom::Array<CoordType, 2> weights_2D(2, 2);
  weights_2D(0, 0) = 1.0;
  weights_2D(0, 1) = 2.0;
  weights_2D(1, 0) = 1.0;
  weights_2D(1, 1) = 0.5;

  BezierPatchType nonrational_patch_2D(controlPoints_2D, order_u, order_v);
  EXPECT_EQ(nonrational_patch_2D.getOrder_u(), order_u);
  EXPECT_EQ(nonrational_patch_2D.getOrder_v(), order_v);
  EXPECT_FALSE(nonrational_patch_2D.isRational());

  for(int p = 0; p <= nonrational_patch_2D.getOrder_u(); ++p)
    for(int q = 0; q <= nonrational_patch_2D.getOrder_v(); ++q)
      for(int i = 0; i < DIM; ++i)
        EXPECT_DOUBLE_EQ(controlPoints_2D(p, q)[i],
                         nonrational_patch_2D(p, q)[i]);

  BezierPatchType rational_patch_2D(controlPoints_2D, weights_2D, order_u, order_v);
  EXPECT_EQ(rational_patch_2D.getOrder_u(), order_u);
  EXPECT_EQ(rational_patch_2D.getOrder_v(), order_v);
  EXPECT_TRUE(rational_patch_2D.isRational());

  for(int p = 0; p <= rational_patch_2D.getOrder_u(); ++p)
    for(int q = 0; q <= rational_patch_2D.getOrder_v(); ++q)
    {
      for(int i = 0; i < DIM; ++i)
        EXPECT_DOUBLE_EQ(controlPoints_2D(p, q)[i], rational_patch_2D(p, q)[i]);
      EXPECT_DOUBLE_EQ(weights_2D(p, q), rational_patch_2D.getWeight(p, q));
    }
}

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, make_rational)
{
  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType>;

  const int order_u = 1;
  const int order_v = 1;

  SLIC_INFO("Test rational property");

  PointType controlPoints[4] = {PointType {0.0, 0.0, 1.0},
                                PointType {0.0, 1.0, 0.0},
                                PointType {1.0, 0.0, 0.0},
                                PointType {1.0, 1.0, -1.0}};

  BezierPatchType rPatch(controlPoints, order_u, order_v);
  EXPECT_FALSE(rPatch.isRational());
  EXPECT_EQ(rPatch.getWeights().size(), 0);

  rPatch.makeRational();
  EXPECT_TRUE(rPatch.isRational());
  EXPECT_EQ(rPatch.getWeights().shape(), rPatch.getControlPoints().shape());

  // makeRational should set all weights to 1
  for(int p = 0; p <= order_u; ++p)
    for(int q = 0; q <= order_v; ++q)
      EXPECT_DOUBLE_EQ(rPatch.getWeight(p, q), 1.0);

  // With all weights 1, the surface should be the same as if unweighted
  BezierPatchType bPatch(controlPoints, order_u, order_v);
  for(double u = 0; u <= 1.0; u += 0.1)
    for(double v = 0; v <= 1.0; v += 0.1)
    {
      auto pt1 = rPatch.evaluate(u, v);
      auto pt2 = bPatch.evaluate(u, v);
      for(int i = 0; i < DIM; ++i) EXPECT_NEAR(pt1[i], pt2[i], 1e-10);
    }

  rPatch.makeNonrational();
  EXPECT_FALSE(rPatch.isRational());
  EXPECT_EQ(rPatch.getWeights().size(), 0);
}

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, evaluate)
{
  SLIC_INFO("Testing bezier Patch evaluation");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType>;

  const int order_u = 3;
  const int order_v = 2;

  // clang-format off
  PointType controlPoints[(order_u + 1) * (order_v + 1)] = {
    PointType {0, 0, 0}, PointType {0, 4,  0}, PointType {0, 8, -3},
    PointType {2, 0, 6}, PointType {2, 4,  0}, PointType {2, 8,  0},
    PointType {4, 0, 0}, PointType {4, 4,  0}, PointType {4, 8,  3},
    PointType {6, 0, 0}, PointType {6, 4, -3}, PointType {6, 8,  0}};
  // clang-format on

  BezierPatchType bPatch(controlPoints, order_u, order_v);

  // Evaluate the patch at each of the four corners, where the patch interpolates
  for(int i = 0; i < DIM; ++i)
  {
    EXPECT_DOUBLE_EQ(bPatch.evaluate(0, 0)[i], controlPoints[0][i]);
    EXPECT_DOUBLE_EQ(bPatch.evaluate(0, 1)[i], controlPoints[2][i]);
    EXPECT_DOUBLE_EQ(bPatch.evaluate(1, 0)[i], controlPoints[9][i]);
    EXPECT_DOUBLE_EQ(bPatch.evaluate(1, 1)[i], controlPoints[11][i]);
  }

  // Evaluate the patch at some interior points
  PointType interior_1 {3.0, 4.0, 0.5625};
  PointType interior_2 {1.5, 6.0, -171.0 / 512.0};
  PointType interior_3 {4.5, 2.0, 39.0 / 512.0};

  for(int i = 0; i < DIM; ++i)
  {
    EXPECT_DOUBLE_EQ(bPatch.evaluate(0.5, 0.5)[i], interior_1[i]);
    EXPECT_DOUBLE_EQ(bPatch.evaluate(0.25, 0.75)[i], interior_2[i]);
    EXPECT_DOUBLE_EQ(bPatch.evaluate(0.75, 0.25)[i], interior_3[i]);
  }

  // Evaluate calculation is asymmetric across axes. Swap them, and repeat the tests.
  bPatch.swapAxes();
  EXPECT_EQ(bPatch.getOrder_u(), order_v);
  EXPECT_EQ(bPatch.getOrder_v(), order_u);

  for(int i = 0; i < DIM; ++i)
  {
    EXPECT_DOUBLE_EQ(bPatch.evaluate(0, 0)[i], controlPoints[0][i]);
    EXPECT_DOUBLE_EQ(bPatch.evaluate(0, 1)[i], controlPoints[9][i]);
    EXPECT_DOUBLE_EQ(bPatch.evaluate(1, 0)[i], controlPoints[2][i]);
    EXPECT_DOUBLE_EQ(bPatch.evaluate(1, 1)[i], controlPoints[11][i]);
  }

  for(int i = 0; i < DIM; ++i)
  {
    EXPECT_DOUBLE_EQ(bPatch.evaluate(0.5, 0.5)[i], interior_1[i]);
    EXPECT_DOUBLE_EQ(bPatch.evaluate(0.75, 0.25)[i], interior_2[i]);
    EXPECT_DOUBLE_EQ(bPatch.evaluate(0.25, 0.75)[i], interior_3[i]);
  }
}

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, evaluate_tall)
{
  SLIC_INFO("Testing bezier Patch evaluation with axes swapped");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType>;

  const int order_u = 2;
  const int order_v = 3;

  // clang-format off
  PointType controlPoints[(order_u + 1) * (order_v + 1)] = {
    PointType {0, 0,  0}, PointType {2, 0, 6}, PointType {4, 0, 0}, PointType {6, 0, 0},
    PointType {0, 4,  0}, PointType {2, 4, 0}, PointType {4, 4, 0}, PointType {6, 4, -3},
    PointType {0, 8, -3}, PointType {2, 8, 0}, PointType {4, 8, 3}, PointType {6, 8,  0}};
  // clang-format on

  BezierPatchType bPatch(controlPoints, order_u, order_v);

  // Evaluate the patch at each of the four corners, where the patch interpolates
  for(int i = 0; i < DIM; ++i)
  {
    EXPECT_DOUBLE_EQ(bPatch.evaluate(0, 0)[i], controlPoints[0][i]);
    EXPECT_DOUBLE_EQ(bPatch.evaluate(0, 1)[i], controlPoints[3][i]);
    EXPECT_DOUBLE_EQ(bPatch.evaluate(1, 0)[i], controlPoints[8][i]);
    EXPECT_DOUBLE_EQ(bPatch.evaluate(1, 1)[i], controlPoints[11][i]);
  }

  // Evaluate the patch at some interior points
  PointType interior_1 {3.0, 4.0, 0.5625};          // (0.5, 0.5)
  PointType interior_2 {1.5, 6.0, -171.0 / 512.0};  // (0.75, 0.25)
  PointType interior_3 {4.5, 2.0, 39.0 / 512.0};    // (0.25, 0.75)

  for(int i = 0; i < DIM; ++i)
  {
    EXPECT_DOUBLE_EQ(bPatch.evaluate(0.5, 0.5)[i], interior_1[i]);
    EXPECT_DOUBLE_EQ(bPatch.evaluate(0.75, 0.25)[i], interior_2[i]);
    EXPECT_DOUBLE_EQ(bPatch.evaluate(0.25, 0.75)[i], interior_3[i]);
  }
}

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, evaluation_degenerate)
{
  SLIC_INFO("Testing Bezier patch evaluation with one degenerate axis");
  // Should reduce to a Bezier curve along the nonempty dimension

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType>;

  const int order = 3;
  PointType data[order + 1] = {PointType {0.6, 1.2, 1.0},
                               PointType {1.3, 1.6, 1.8},
                               PointType {2.9, 2.4, 2.3},
                               PointType {3.2, 3.5, 3.0}};

  BezierCurveType bCurve(data, order);
  BezierPatchType bPatch(data, order, 0);

  for(double t = 0; t <= 1; t += 0.01)
    for(int i = 0; i < DIM; ++i)
    {
      EXPECT_NEAR(bCurve.evaluate(t)[i], bPatch.evaluate(t, 0)[i], 1e-10);
      EXPECT_NEAR(bCurve.evaluate(t)[i], bPatch.evaluate(t, 0.5)[i], 1e-10);
      EXPECT_NEAR(bCurve.evaluate(t)[i], bPatch.evaluate(t, 1.0)[i], 1e-10);
    }
}

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, tangent)
{
  SLIC_INFO("Testing bezier Patch tangent calculation");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType>;

  const int order_u = 3;
  const int order_v = 2;

  // clang-format off
  PointType controlPoints[(order_u + 1) * (order_v + 1)] = {
    PointType {0, 0, 0}, PointType {0, 4,  0}, PointType {0, 8, -3},
    PointType {2, 0, 6}, PointType {2, 4,  0}, PointType {2, 8,  0},
    PointType {4, 0, 0}, PointType {4, 4,  0}, PointType {4, 8,  3},
    PointType {6, 0, 0}, PointType {6, 4, -3}, PointType {6, 8,  0}};
  // clang-format on

  BezierPatchType bPatch(controlPoints, order_u, order_v);

  // Evaluate the tangent at the lower corfour corners,
  //  where the tangents are defined by the control nodes
  VectorType node00_u = VectorType(controlPoints[0], controlPoints[3]);
  VectorType node00_v = VectorType(controlPoints[0], controlPoints[1]);

  for(int i = 0; i < DIM; ++i)
  {
    EXPECT_NEAR(order_u * node00_u[i], bPatch.dt(0, 0, 0)[i], 1e-10);
    EXPECT_NEAR(order_v * node00_v[i], bPatch.dt(0, 0, 1)[i], 1e-10);
  }

  // Evaluate the tangent at some interior nodes
  VectorType interior1_u = {6.0, 0.0, 567.0 / 128.0};
  VectorType interior1_v = {0.0, 8.0, -159.0 / 64.0};

  VectorType interior2_u = {6.0, 0.0, -657.0 / 128.0};
  VectorType interior2_v = {0.0, 8.0, -123.0 / 64.0};

  for(int i = 0; i < DIM; ++i)
  {
    EXPECT_NEAR(interior1_u[i], bPatch.dt(0.25, 0.75, 0)[i], 1e-10);
    EXPECT_NEAR(interior1_v[i], bPatch.dt(0.25, 0.75, 1)[i], 1e-10);

    EXPECT_NEAR(interior2_u[i], bPatch.dt(0.75, 0.25, 0)[i], 1e-10);
    EXPECT_NEAR(interior2_v[i], bPatch.dt(0.75, 0.25, 1)[i], 1e-10);
  }
}

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, normal)
{
  SLIC_INFO("Testing bezier Patch normal calculation");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType>;

  const int order_u = 3;
  const int order_v = 2;

  // clang-format off
  PointType controlPoints[(order_u + 1) * (order_v + 1)] = {
    PointType {0, 0, 0}, PointType {0, 4,  0}, PointType {0, 8, -3},
    PointType {2, 0, 6}, PointType {2, 4,  0}, PointType {2, 8,  0},
    PointType {4, 0, 0}, PointType {4, 4,  0}, PointType {4, 8,  3},
    PointType {6, 0, 0}, PointType {6, 4, -3}, PointType {6, 8,  0}};
  // clang-format on

  BezierPatchType bPatch(controlPoints, order_u, order_v);

  // Test on some interior nodes
  VectorType interior_1 = {-279.0 / 16.0, 819.0 / 32.0, 48.0};  // (0.25, 0.25)
  VectorType interior_2 = {-567.0 / 16.0, 477.0 / 32.0, 48.0};  // (0.25, 0.75)
  VectorType interior_3 = {657.0 / 16.0, 369.0 / 32.0, 48.0};   // (0.75, 0.25)
  VectorType interior_4 = {369.0 / 16.0, -513.0 / 32.0, 48.0};  // (0.75, 0.75)

  for(int i = 0; i < DIM; ++i)
  {
    EXPECT_NEAR(interior_1[i], bPatch.normal(0.25, 0.25)[i], 1e-10);
    EXPECT_NEAR(interior_2[i], bPatch.normal(0.25, 0.75)[i], 1e-10);
    EXPECT_NEAR(interior_3[i], bPatch.normal(0.75, 0.25)[i], 1e-10);
    EXPECT_NEAR(interior_4[i], bPatch.normal(0.75, 0.75)[i], 1e-10);
  }
}

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, split_degenerate)
{
  SLIC_INFO("Testing bezier Patch splitting for order 0 surface");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType>;

  const int order_u = 0;
  const int order_v = 0;

  PointType data[(order_u + 1) * (order_v + 1)] = {
    PointType::make_point(0.6, 1.2, 0.2)};

  BezierPatchType p(data, 0, 0);

  BezierPatchType p1, p2, p3, p4;

  p.split(0.5, 0.5, p1, p2, p3, p4);

  for(int i = 0; i < DIM; i++)
  {
    EXPECT_DOUBLE_EQ(data[0][i], p1(0, 0)[i]);
    EXPECT_DOUBLE_EQ(data[0][i], p2(0, 0)[i]);
    EXPECT_DOUBLE_EQ(data[0][i], p3(0, 0)[i]);
    EXPECT_DOUBLE_EQ(data[0][i], p4(0, 0)[i]);
  }
}

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, split_curve)
{
  SLIC_INFO("Testing Bezier patch split with one degenerate axis");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierCurveType = primal::BezierCurve<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType>;

  const int order = 3;
  PointType data[order + 1] = {PointType {0.6, 1.2, 1.0},
                               PointType {1.3, 1.6, 1.8},
                               PointType {2.9, 2.4, 2.3},
                               PointType {3.2, 3.5, 3.0}};

  BezierCurveType bCurve(data, order);
  BezierPatchType bPatch(data, order, 0);

  BezierCurveType c1, c2;
  BezierPatchType p1, p2, p3, p4;

  bCurve.split(0.5, c1, c2);
  bPatch.split(0.5, 0.5, p1, p2, p3, p4);

  // Components split along the u-axis p1/p3 vs p2/p4 should equal curves
  // Components split along the v-axis p1/p2 vs p3/p4 should equal each other
  for(double u = 0; u <= 1; u += 0.1)
    for(double v = 0; v <= 1; v += 0.1)
      for(int i = 0; i < DIM; ++i)
      {
        EXPECT_NEAR(c1.evaluate(u)[i], p1.evaluate(u, v)[i], 1e-10);
        EXPECT_NEAR(c1.evaluate(u)[i], p3.evaluate(u, v)[i], 1e-10);

        EXPECT_NEAR(c2.evaluate(u)[i], p2.evaluate(u, v)[i], 1e-10);
        EXPECT_NEAR(c2.evaluate(u)[i], p4.evaluate(u, v)[i], 1e-10);
      }
}

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, split_plane)
{
  SLIC_INFO("Testing Bezier patch evaluation with one degenerate axis");
  // Should reduce to a Bezier curve along the nonempty dimension

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType>;

  const int order_u = 1;
  const int order_v = 1;

  // clang-format off
  PointType controlPoints[(order_u + 1) * (order_v + 1)] = {
    PointType {1.0, 0.0, 0.0}, PointType {0.0, 1.0, 0.0},
    PointType {0.0, 0.0, 1.0}, PointType {0.0, 0.0, 1.0}};
  // clang-format on

  BezierPatchType p0(controlPoints, order_u, order_v), p_arr[4];

  p0.split(0.5, 0.5, p_arr[0], p_arr[1], p_arr[2], p_arr[3]);

  // Ensure that after splitting, each point remains on the same hyperplane
  for(double u = 0.0; u <= 1.0; u += 0.1)
    for(double v = 0.0; v <= 1.0; v += 0.1)
      for(int n = 0; n < 4; ++n)
      {
        auto pt = p_arr[n].evaluate(u, v);
        EXPECT_NEAR(pt[0] + pt[1] + pt[2], 1.0, 1e-10);
      }
}

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, split_patch)
{
  SLIC_INFO("Testing Bezier patch evaluation with one degenerate axis");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType>;

  const int order_u = 3;
  const int order_v = 2;

  // clang-format off
  PointType controlPoints[(order_u + 1) * (order_v + 1)] = {
    PointType {0, 0, 0}, PointType {0, 4,  0}, PointType {0, 8, -3},
    PointType {2, 0, 6}, PointType {2, 4,  0}, PointType {2, 8,  0},
    PointType {4, 0, 0}, PointType {4, 4,  0}, PointType {4, 8,  3},
    PointType {6, 0, 0}, PointType {6, 4, -3}, PointType {6, 8,  0}};
  // clang-format on

  BezierPatchType p0(controlPoints, order_u, order_v), p1, p2, p3, p4;
  p0.split(0.5, 0.5, p1, p2, p3, p4);

  // Check that evaluation the appropriate edges line up
  for(double t = 0; t <= 1; t += 0.1)
  {
    for(int i = 0; i < DIM; i++)
    {
      EXPECT_NEAR(p1.evaluate(1, t)[i], p2.evaluate(0, t)[i], 1e-10);
      EXPECT_NEAR(p3.evaluate(1, t)[i], p4.evaluate(0, t)[i], 1e-10);

      EXPECT_NEAR(p1.evaluate(t, 1)[i], p3.evaluate(t, 0)[i], 1e-10);
      EXPECT_NEAR(p2.evaluate(t, 1)[i], p4.evaluate(t, 0)[i], 1e-10);
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, isPlanar)
{
  SLIC_INFO("Testing Bezier patch evaluation with one degenerate axis");
  // Should reduce to a Bezier curve along the nonempty dimension

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType>;

  // order (0, 0) -- always true
  {
    const int order_u = 0;
    const int order_v = 0;
    BezierPatchType patch(order_u, order_v);
    EXPECT_TRUE(patch.isPlanar());

    patch(0, 0) = PointType::make_point(1.0, 1.0, 1.0);
    EXPECT_TRUE(patch.isPlanar());
  }

  // order (0, 1), (1, 0) -- always true
  {
    PointType data[2] = {PointType {0.0, 1.0, 2.0}, PointType {2.0, 1.0, 0.0}};
    BezierPatchType patch1(data, 0, 1), patch2(data, 1, 0);

    EXPECT_TRUE(patch1.isPlanar());
    EXPECT_TRUE(patch2.isPlanar());
  }

  // order (1, 1) -- sometimes true
  {
    // clang-format off
    PointType data_planar[4] = {
      PointType {1.0, 0.0, 0.0}, PointType {0.0, 1.0, 0.0},
      PointType {0.0, 0.0, 1.0}, PointType {-1.0, 1.0, 1.0}};

    PointType data_nonplanar[4] = {
      PointType {1.0, 0.0, 0.0}, PointType {0.0, 1.0, 0.0},
      PointType {0.0, 0.0, 1.0}, PointType {0.0, 0.0, 2.0}};
    // clang-format on

    EXPECT_TRUE(BezierPatchType(data_planar, 1, 1).isPlanar());
    EXPECT_FALSE(BezierPatchType(data_nonplanar, 1, 1).isPlanar());
  }

  // order (2, 2)
  {
    // clang-format off
    PointType data[9] = {
      PointType {-1.0, -1.0, 0.0}, PointType {-1.0, 0.0, 0.0}, PointType {-1.0, 1.0, 0.0},  
      PointType { 0.0, -1.0, 0.0}, PointType { 0.0, 0.0, 0.0}, PointType { 0.0, 1.0, 0.0},
      PointType { 1.0, -1.0, 0.0}, PointType { 1.0, 0.0, 0.0}, PointType { 1.0, 1.0, 0.0}};
    // clang-format on

    const int ord_u = 2;
    const int ord_v = 2;

    BezierPatchType patch(data, ord_u, ord_v);

    // Begins as planar
    EXPECT_TRUE(patch.isPlanar());

    // Move middle point and check linearity with diffrent tolerances
    patch(1, 1).array() += 0.005 * VectorType({0.0, 0.0, 1.0}).array();

    // Linear for a coarse tolerance
    {
      const double tol = 0.01;
      const double tol_sq = tol * tol;
      EXPECT_TRUE(patch.isPlanar(tol_sq));
    }

    // Non-Linear for a finer tolerance
    {
      const double tol = 0.001;
      const double tol_sq = tol * tol;
      EXPECT_FALSE(patch.isPlanar(tol_sq));
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, reverse_orientation)
{
  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType>;

  SLIC_INFO("Testing reverseOrientation(axis) on Bezier patches");

  const int order_u = 3;
  const int order_v = 2;

  // clang-format off
  PointType controlPoints[(order_u + 1) * (order_v + 1)] = {
    PointType {0, 0, 0}, PointType {0, 4,  0}, PointType {0, 8, -3},
    PointType {2, 0, 6}, PointType {2, 4,  0}, PointType {2, 8,  0},
    PointType {4, 0, 0}, PointType {4, 4,  0}, PointType {4, 8,  3},
    PointType {6, 0, 0}, PointType {6, 4, -3}, PointType {6, 8,  0}};
  // clang-format on

  BezierPatchType original(controlPoints, order_u, order_v);
  BezierPatchType reversed(controlPoints, order_u, order_v);

  // Reverse along the u-axis
  reversed.reverseOrientation(0);
  for(double u = 0; u <= 1.0; u += 0.1)
    for(double v = 0; v <= 1.0; v += 0.1)
    {
      PointType o_pt = original.evaluate(u, v);
      PointType r_pt = reversed.evaluate(1 - u, v);

      for(int i = 0; i < DIM; ++i) EXPECT_NEAR(o_pt[i], r_pt[i], 1e-10);
    }

  // Reverse along the u-axis again, should return to original
  reversed.reverseOrientation(0);
  for(double u = 0; u <= 1.0; u += 0.1)
    for(double v = 0; v <= 1.0; v += 0.1)
    {
      PointType o_pt = original.evaluate(u, v);
      PointType r_pt = reversed.evaluate(u, v);

      for(int i = 0; i < DIM; ++i) EXPECT_NEAR(o_pt[i], r_pt[i], 1e-10);
    }

  // Reverse along the v-axis
  reversed.reverseOrientation(1);
  for(double u = 0; u <= 1.0; u += 0.1)
    for(double v = 0; v <= 1.0; v += 0.1)
    {
      PointType o_pt = original.evaluate(u, v);
      PointType r_pt = reversed.evaluate(u, 1 - v);

      for(int i = 0; i < DIM; ++i) EXPECT_NEAR(o_pt[i], r_pt[i], 1e-10);
    }

  // Reverse along the u-axis again
  reversed.reverseOrientation(0);
  for(double u = 0; u <= 1.0; u += 0.1)
    for(double v = 0; v <= 1.0; v += 0.1)
    {
      PointType o_pt = original.evaluate(u, v);
      PointType r_pt = reversed.evaluate(1 - u, 1 - v);

      for(int i = 0; i < DIM; ++i) EXPECT_NEAR(o_pt[i], r_pt[i], 1e-10);
    }
}

//------------------------------------------------------------------------------
TEST(primal_bezierpatch, rational_evaluation_split)
{
  // Test behavior with a spherical Bezier patch
  SLIC_INFO("Testing Bezier patch evaluation wirh rational weights");

  const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType>;

  const int ord_u = 3;
  const int ord_v = 3;

  // clang-format off
  PointType data[16] = {  
    PointType {0, 0,  1}, PointType {0, 0,  1}, PointType { 0, 0,  1}, PointType { 0, 0,  1},
    PointType {2, 0,  1}, PointType {2, 4,  1}, PointType {-2, 4,  1}, PointType {-2, 0,  1},
    PointType {2, 0, -1}, PointType {2, 4, -1}, PointType {-2, 4, -1}, PointType {-2, 0, -1},
    PointType {0, 0, -1}, PointType {0, 0, -1}, PointType { 0, 0, -1}, PointType { 0, 0, -1}};

  double weights[16] = {1.0,     1.0/3.0, 1.0/3.0, 1.0,
                        1.0/3.0, 1.0/9.0, 1.0/9.0, 1.0/3.0,
                        1.0/3.0, 1.0/9.0, 1.0/9.0, 1.0/3.0,
                        1.0,     1.0/3.0, 1.0/3.0, 1.0};
  // clang-format on

  BezierPatchType hemisphere(data, weights, ord_u, ord_v);

  // Verify that evaluation points are on the sphere
  for(double u = 0; u <= 1.0; u += 0.1)
    for(double v = 0; v <= 1.0; v += 0.1)
    {
      PointType pt = hemisphere.evaluate(u, v);
      EXPECT_NEAR(pt[0] * pt[0] + pt[1] * pt[1] + pt[2] * pt[2], 1.0, 1e-10);
    }

  BezierPatchType patches[4];
  hemisphere.split(0.5, 0.5, patches[0], patches[1], patches[2], patches[3]);

  // Verify that evaluation points are on still on the sphere for each subpatch
  for(double u = 0; u <= 1.0; u += 0.1)
    for(double v = 0; v <= 1.0; v += 0.1)
      for(int n = 0; n < 4; ++n)
      {
        PointType pt = patches[n].evaluate(u, v);
        EXPECT_NEAR(pt[0] * pt[0] + pt[1] * pt[1] + pt[2] * pt[2], 1.0, 1e-10);
      }

  // Do it again
  BezierPatchType sub_patches[4];
  patches[0].split(0.5,
                   0.5,
                   sub_patches[0],
                   sub_patches[1],
                   sub_patches[2],
                   sub_patches[3]);

  for(double u = 0; u <= 1.0; u += 0.1)
    for(double v = 0; v <= 1.0; v += 0.1)
      for(int n = 0; n < 4; ++n)
      {
        PointType pt = sub_patches[n].evaluate(u, v);
        EXPECT_NEAR(pt[0] * pt[0] + pt[1] * pt[1] + pt[2] * pt[2], 1.0, 1e-10);
      }
}

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}

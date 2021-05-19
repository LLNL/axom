// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/Triangle.hpp"

#include "axom/primal/operators/clip.hpp"

#include <limits>

namespace Primal3D
{
typedef axom::primal::Point<double, 3> PointType;
typedef axom::primal::Vector<double, 3> VectorType;
typedef axom::primal::BoundingBox<double, 3> BoundingBoxType;
typedef axom::primal::Triangle<double, 3> TriangleType;
typedef axom::primal::Polygon<double, 3> PolygonType;
}  // namespace Primal3D

TEST(primal_clip, simple_clip)
{
  using namespace Primal3D;
  BoundingBoxType bbox;
  bbox.addPoint(PointType::zero());
  bbox.addPoint(PointType::ones());

  PointType points[] = {
    PointType::make_point(2, 2, 2),
    PointType::make_point(2, 2, 4),
    PointType::make_point(2, 4, 2),
    PointType::make_point(-100, -100, 0.5),
    PointType::make_point(-100, 100, 0.5),
    PointType::make_point(100, 0, 0.5),
    PointType::make_point(0.25, 0.25, 0.5),
    PointType::make_point(0.75, 0.25, 0.5),
    PointType::make_point(0.66, 0.5, 0.5),
    PointType::make_point(1.5, 0.5, 0.5),
  };

  {
    TriangleType tri(points[0], points[1], points[2]);

    PolygonType poly = axom::primal::clip(tri, bbox);
    EXPECT_EQ(0, poly.numVertices());
  }

  {
    TriangleType tri(points[3], points[4], points[5]);

    PolygonType poly = axom::primal::clip(tri, bbox);
    EXPECT_EQ(4, poly.numVertices());

    SLIC_INFO("Intersection of triangle " << tri << " and bounding box " << bbox
                                          << " is polygon" << poly);
  }

  {
    TriangleType tri(points[3], points[4], points[5]);

    PolygonType poly = axom::primal::clip(tri, bbox);
    EXPECT_EQ(4, poly.numVertices());

    EXPECT_EQ(PointType(.5), poly.centroid());

    SLIC_INFO("Intersection of triangle " << tri << " and bounding box " << bbox
                                          << " is polygon" << poly);
  }

  {
    TriangleType tri(points[6], points[7], points[9]);

    PolygonType poly = axom::primal::clip(tri, bbox);
    EXPECT_EQ(4, poly.numVertices());

    SLIC_INFO("Intersection of triangle " << tri << " and bounding box " << bbox
                                          << " is polygon" << poly);
  }
}

TEST(primal_clip, unit_simplex)
{
  using namespace Primal3D;
  double delta = 1e-5;

  // Test the "unit simplex", and a jittered version
  PointType points[] = {PointType::make_point(1, 0, 0),
                        PointType::make_point(0, 1, 0),
                        PointType::make_point(0, 0, 1),
                        PointType::make_point(1 + delta, delta, delta),
                        PointType::make_point(delta, 1 + delta, delta),
                        PointType::make_point(delta, delta, 1 + delta)};

  BoundingBoxType bbox;
  bbox.addPoint(PointType::zero());
  bbox.addPoint(PointType(.75));

  // intersection of this triangle and cube is a hexagon
  {
    TriangleType tri(points[0], points[1], points[2]);

    PolygonType poly = axom::primal::clip(tri, bbox);
    EXPECT_EQ(6, poly.numVertices());

    SLIC_INFO("Intersection of triangle " << tri << " and bounding box " << bbox
                                          << " is polygon" << poly);
  }
  {
    TriangleType tri(points[3], points[4], points[5]);

    PolygonType poly = axom::primal::clip(tri, bbox);
    EXPECT_EQ(6, poly.numVertices());

    SLIC_INFO("Intersection of triangle " << tri << " and bounding box " << bbox
                                          << " is polygon" << poly);
  }
}

TEST(primal_clip, boundingBoxOptimization)
{
  using namespace Primal3D;

  SLIC_INFO("Checking correctness of optimization for skipping clipping "
            << " of planes that the triangle's bounding box doesn't cover");

  const double VAL1 = 3.;
  const double VAL2 = 2.;

  BoundingBoxType bbox;
  bbox.addPoint(PointType(-1.));
  bbox.addPoint(PointType(1.));

  PointType midpoint = PointType::zero();

  PointType points[] = {
    PointType::make_point(VAL1, VAL2, 0),
    PointType::make_point(-VAL1, VAL2, 0),
    PointType::make_point(VAL1, -VAL2, 0),
    PointType::make_point(-VAL1, -VAL2, 0),

    PointType::make_point(VAL1, 0, VAL2),
    PointType::make_point(-VAL1, 0, VAL2),
    PointType::make_point(VAL1, 0, -VAL2),
    PointType::make_point(-VAL1, 0, -VAL2),

    PointType::make_point(0, VAL2, VAL1),
    PointType::make_point(0, VAL2, -VAL1),
    PointType::make_point(0, -VAL2, VAL1),
    PointType::make_point(0, -VAL2, -VAL1),

    PointType::make_point(0, VAL1, VAL2),
    PointType::make_point(0, -VAL1, VAL2),
    PointType::make_point(0, VAL1, -VAL2),
    PointType::make_point(0, -VAL1, -VAL2),
  };

  for(int i = 0; i < 16; i += 2)
  {
    TriangleType tri(midpoint, points[i], points[i + 1]);
    PolygonType poly = axom::primal::clip(tri, bbox);
    SLIC_INFO(poly);
    EXPECT_EQ(5, poly.numVertices());
  }
}

TEST(primal_clip, experimentalData)
{
  using namespace Primal3D;

  const double EPS = 1e-8;

  // Triangle 248 from sphere mesh
  TriangleType tri(PointType::make_point(0.405431, 3.91921, 3.07821),
                   PointType::make_point(1.06511, 3.96325, 2.85626),
                   PointType::make_point(0.656002, 4.32465, 2.42221));

  // Block index {grid pt: (19,29,24); level: 5} from InOutOctree
  BoundingBoxType box12(PointType::make_point(0.937594, 4.06291, 2.50025),
                        PointType::make_point(1.25012, 4.37544, 2.81278));

  PolygonType poly = axom::primal::clip(tri, box12);
  EXPECT_EQ(3, poly.numVertices());

  SLIC_INFO("Intersection of triangle "
            << tri << " \n\t and bounding box " << box12 << " \n\t is polygon"
            << poly << " with centroid " << poly.centroid());

  // Check that the polygon vertices are on the triangle
  for(int i = 0; i < poly.numVertices(); ++i)
  {
    PointType bary = tri.physToBarycentric(poly[i]);
    PointType reconstructed = tri.baryToPhysical(bary);

    SLIC_INFO("Testing clipped polygon point "
              << poly[i] << "-- w/ barycentric coords " << bary
              << "\n\t-- reconstructed point is: " << reconstructed << "...\n");

    double barySum = bary[0] + bary[1] + bary[2];
    EXPECT_NEAR(1., barySum, EPS);

    for(int dim = 0; dim < 3; ++dim)
    {
      EXPECT_GE(bary[dim], -EPS);
      EXPECT_NEAR(poly[i][dim], reconstructed[dim], EPS);
    }

    EXPECT_TRUE(box12.contains(poly[i]));
  }

  // Check that the polygon centroid is on the triangle
  {
    PointType centroid = poly.centroid();
    PointType bary = tri.physToBarycentric(centroid);
    PointType reconstructed = tri.baryToPhysical(bary);

    SLIC_INFO("Testing clipped polygon centroid "
              << centroid << "-- w/ barycentric coords " << bary
              << "\n\t-- reconstructed point is: " << reconstructed << "...\n");

    double barySum = bary[0] + bary[1] + bary[2];
    EXPECT_NEAR(1., barySum, EPS);

    for(int dim = 0; dim < 3; ++dim)
    {
      EXPECT_GE(bary[dim], -EPS);
      EXPECT_NEAR(centroid[dim], reconstructed[dim], EPS);
    }
    EXPECT_TRUE(box12.contains(centroid));
  }
}

//------------------------------------------------------------------------------
#include "axom/slic/core/SimpleLogger.hpp"
using axom::slic::SimpleLogger;

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);

  SimpleLogger logger;  // create & initialize test logger,

  int result = RUN_ALL_TESTS();

  return result;
}

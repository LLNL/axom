// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "fmt/fmt.hpp"

#include "axom/config.hpp"
#include "axom/slic/interface/slic.hpp"

#include "axom/primal/geometry/OrientedBoundingBox.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Ray.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/Vector.hpp"

#include "axom/primal/operators/intersect.hpp"

#include <cmath>  // for sqrt

using namespace axom;

namespace
{
double randomDouble(double beg = 0., double end = 1.)
{
  double range = end - beg;

  if(range == 0)
  {
    range = 1.;
  }

  return beg + (rand() / (RAND_MAX / range));
}

template <int NDIMS>
primal::Point<double, NDIMS> randomPt(double beg, double end)
{
  primal::Point<double, NDIMS> pt;
  for(int i = 0; i < NDIMS; ++i)
  {
    pt[i] = randomDouble(beg, end);
  }

  return pt;
}

}  // namespace

template <int DIM>
primal::Triangle<double, DIM> roll(const primal::Triangle<double, DIM>& t,
                                   const int i)
{
  return primal::Triangle<double, DIM>(t[i % 3], t[(i + 1) % 3], t[(i + 2) % 3]);
}

template <int DIM>
void permuteCornersTest(const primal::Triangle<double, DIM>& a,
                        const primal::Triangle<double, DIM>& b,
                        const std::string& whattest,
                        const bool includeBdry,
                        const bool testtrue,
                        double EPS = 1E-8)
{
  SCOPED_TRACE(
    whattest +
    (includeBdry ? " (including boundary)" : " (NOT including boundary)"));

  int numFailures = 0;
  const int numTests = 3 * 3 * 4;

  for(int i = 0; i < 3; ++i)
  {
    for(int j = 0; j < 3; ++j)
    {
      auto t1 = roll(a, i);
      auto t2 = roll(b, j);
      if(primal::intersect(t1, t2, includeBdry, EPS) != testtrue)
      {
        ++numFailures;
      }
    }
  }

  const primal::Triangle<double, DIM> ap(a[0], a[2], a[1]);
  const primal::Triangle<double, DIM> bp(b[0], b[2], b[1]);

  for(int i = 0; i < 3; ++i)
  {
    for(int j = 0; j < 3; ++j)
    {
      auto t1 = roll(ap, i);
      auto t2 = roll(bp, j);
      if(primal::intersect(t1, t2, includeBdry, EPS) != testtrue)
      {
        ++numFailures;
      }
    }
  }

  // Now swap a for b

  for(int i = 0; i < 3; ++i)
  {
    for(int j = 0; j < 3; ++j)
    {
      auto t1 = roll(b, i);
      auto t2 = roll(a, j);
      if(primal::intersect(t1, t2, includeBdry, EPS) != testtrue)
      {
        ++numFailures;
      }
    }
  }

  for(int i = 0; i < 3; ++i)
  {
    for(int j = 0; j < 3; ++j)
    {
      auto t1 = roll(bp, i);
      auto t2 = roll(ap, j);
      if(primal::intersect(t1, t2, includeBdry, EPS) != testtrue)
      {
        ++numFailures;
      }
    }
  }

  if(numFailures == 0)
  {
    SUCCEED();
  }
  else
  {
    ADD_FAILURE()
      << (testtrue ? "Triangles should intersect but did not"
                   : "Triangles should not intersect but did")
      << "\n\t -- " << numFailures << " of " << numTests
      << " permutations failed when testing intersection of triangles "
      << "\n\t a: " << a << "\n\t b: " << b;
  }
}

TEST(primal_intersect, ray_segment_intersection)
{
  using PointType = primal::Point<double, 2>;
  using SegmentType = primal::Segment<double, 2>;
  using VectorType = primal::Vector<double, 2>;
  using RayType = primal::Ray<double, 2>;

  // STEP 0: construct segment
  PointType A {0, 0};
  PointType B {1, 1};
  SegmentType S(A, B);

  // STEP 1: construct ray
  PointType origin {0.5, -0.5};
  VectorType direction {0.0, 0.5};
  RayType R(origin, direction);

  // compute intersection
  {
    PointType ip;
    bool intersects1a = primal::intersect(R, S, ip);
    EXPECT_TRUE(intersects1a);

    double ray_param_b;
    bool intersects1b = primal::intersect(R, S, ray_param_b);
    EXPECT_TRUE(intersects1b);
    PointType interpolatedIP_b = R.at(ray_param_b);

    double ray_param_c;
    double seg_param_c;
    bool intersects1c = primal::intersect(R, S, ray_param_c, seg_param_c);
    EXPECT_TRUE(intersects1c);
    EXPECT_EQ(ray_param_b, ray_param_c);
    PointType interpolatedIP_c_ray = R.at(ray_param_c);
    PointType interpolatedIP_c_seg = S.at(seg_param_c);

    for(int i = 0; i < 2; ++i)
    {
      EXPECT_DOUBLE_EQ(0.5, ip[i]);
      EXPECT_DOUBLE_EQ(ip[1], interpolatedIP_b[i]);
      EXPECT_DOUBLE_EQ(ip[1], interpolatedIP_c_ray[i]);
      EXPECT_DOUBLE_EQ(ip[1], interpolatedIP_c_seg[i]);
    }
  }

  // construct a non-intersecting ray
  {
    PointType ip;
    double param;
    VectorType oppositeDirection {0, -1};
    RayType R2(origin, oppositeDirection);
    bool intersects2a = primal::intersect(R2, S, ip);
    EXPECT_FALSE(intersects2a);

    bool intersects2b = primal::intersect(R2, S, param);
    EXPECT_FALSE(intersects2b);
  }

  // construct another non-intersecting ray
  {
    PointType ip;
    double param;
    PointType origin3 {0.5, 1.};
    RayType R3(origin3, direction);
    bool intersects3a = primal::intersect(R3, S, ip);
    EXPECT_FALSE(intersects3a);

    bool intersects3b = primal::intersect(R3, S, param);
    EXPECT_FALSE(intersects3b);
  }
}

TEST(primal_intersect, more_ray_segment_intersection)
{
  constexpr int DIM = 2;
  constexpr double EPS = 1e-9;

  using PointType = primal::Point<double, DIM>;
  using SegmentType = primal::Segment<double, DIM>;
  using VectorType = primal::Vector<double, DIM>;
  using RayType = primal::Ray<double, DIM>;

  const int NVALS = 100;

  // Check several values for Ray aligned with x-axis
  for(int seg_pos = -10; seg_pos < 10; ++seg_pos)
  {
    for(int ray_height = 0; ray_height <= NVALS; ++ray_height)
    {
      for(int seg_height = 1; seg_height <= NVALS; ++seg_height)
      {
        const double ray_y = static_cast<double>(ray_height) / NVALS;
        PointType origin {0, ray_y};
        VectorType dir {1, 0};
        RayType ray(origin, dir);

        double x_val = static_cast<double>(seg_pos);
        double y_val = static_cast<double>(seg_height) / NVALS;
        SegmentType seg(PointType {x_val, 0}, PointType {x_val, y_val});

        double ray_param {0};
        double seg_param {0};
        bool expect_intersect = seg_pos >= 0 && ray_height <= seg_height;
        EXPECT_EQ(expect_intersect,
                  primal::intersect(ray, seg, ray_param, seg_param, EPS))
          << "Ray: " << ray << "; seg: " << seg << "; ray_param: " << ray_param
          << "; seg_param: " << seg_param << "; expect_intersect? "
          << expect_intersect;

        if(expect_intersect)
        {
          EXPECT_NEAR(x_val, ray_param, EPS);
          EXPECT_NEAR(ray_y, seg.at(seg_param)[1], EPS);
        }
      }
    }
  }

  // Check several values for ray aligned with y-axis
  for(int seg_pos = -10; seg_pos < 10; ++seg_pos)
  {
    for(int ray_offset = 0; ray_offset <= NVALS; ++ray_offset)
    {
      for(int seg_width = 1; seg_width <= NVALS; ++seg_width)
      {
        PointType origin {static_cast<double>(ray_offset) / NVALS, 0};
        // Note: Ray direction is normalized upon construction so the
        // following scale has no effect on parametrized intersection result
        VectorType dir {0, .5};
        RayType ray(origin, dir);

        double x_val = static_cast<double>(seg_width) / NVALS;
        double y_val = static_cast<double>(seg_pos);
        SegmentType seg(PointType {0, y_val}, PointType {x_val, y_val});

        double param {0};
        bool expect_intersect = seg_pos >= 0 && ray_offset <= seg_width;
        EXPECT_EQ(expect_intersect, primal::intersect(ray, seg, param))
          << "Ray: " << ray << "; seg: " << seg << "; param: " << param
          << "; expect_intersect? " << expect_intersect;
        if(expect_intersect)
        {
          EXPECT_NEAR(y_val, param, EPS);
        }
      }
    }
  }
}

TEST(primal_intersect, triangle_aabb_intersection)
{
  const int DIM = 3;
  typedef primal::Point<double, DIM> PointType;
  typedef primal::Triangle<double, DIM> TriangleType;
  typedef primal::BoundingBox<double, DIM> BoundingBoxType;

  double xArr[3] = {1., 0., 0.};
  double yArr[3] = {0., 1., 0.};
  double zArr[3] = {0., 0., 1.};

  PointType ptX(xArr);
  PointType ptY(yArr);
  PointType ptZ(zArr);

  TriangleType unitTri(ptX, ptY, ptZ);
  BoundingBoxType unitBB(PointType::zero(), PointType::ones());

  EXPECT_TRUE(primal::intersect(unitTri, unitBB));

  // Let's first move the bounding box around
  BoundingBoxType v0_BB(ptX);
  v0_BB.expand(.1);
  SLIC_INFO("Testing v0 bounding box: " << v0_BB << " against unit triangle");
  EXPECT_TRUE(v0_BB.contains(ptX));
  EXPECT_TRUE(primal::intersect(unitTri, v0_BB));

  BoundingBoxType v1_BB(ptY);
  v1_BB.expand(.1);
  SLIC_INFO("Testing v1 bounding box: " << v1_BB << " against unit triangle");
  EXPECT_TRUE(v1_BB.contains(ptY));
  EXPECT_TRUE(primal::intersect(unitTri, v1_BB));

  BoundingBoxType v2_BB(ptZ);
  v2_BB.expand(.1);
  SLIC_INFO("Testing v2 bounding box: " << v2_BB << " against unit triangle");
  EXPECT_TRUE(v2_BB.contains(ptZ));
  EXPECT_TRUE(primal::intersect(unitTri, v2_BB));

  BoundingBoxType mid_BB(PointType::zero());
  mid_BB.addPoint(PointType(0.9));
  SLIC_INFO(
    "Testing bounding box: "
    << mid_BB << " against unit triangle.  Note -- BB should intersect interior of triangle");
  EXPECT_TRUE(primal::intersect(unitTri, mid_BB));

  BoundingBoxType high_BB(PointType::ones());
  high_BB.addPoint(PointType(0.5));
  SLIC_INFO(
    "Testing bounding box: "
    << high_BB << " against unit triangle.  Note -- BB should not intersect interior of triangle");
  EXPECT_FALSE(primal::intersect(unitTri, high_BB));

  BoundingBoxType out_BB(PointType::ones());
  out_BB.addPoint(PointType(2));
  SLIC_INFO(
    "Testing bounding box: "
    << out_BB
    << " against unit triangle.  Note -- BB should not intersect triangle");
  EXPECT_FALSE(primal::intersect(unitTri, out_BB));

  BoundingBoxType negBB(PointType(-5), PointType(-10));
  SLIC_INFO(
    "Testing bounding box: "
    << negBB
    << " against unit triangle.  Note -- BB should not intersect triangle");
  EXPECT_FALSE(primal::intersect(unitTri, negBB));

  // Test new triangle whose edge crosses the BB
  double t2_0[3] = {10., 0., 0.};
  double t2_1[3] = {-10., 0., 0.};
  double t2_2[3] = {0., 100., 0};

  TriangleType xyTri(t2_0, t2_1, t2_2);
  BoundingBoxType bbOrigin(PointType::zero());
  bbOrigin.expand(1.);
  SLIC_INFO("Testing bounding box: "
            << bbOrigin << " against triangle " << xyTri
            << ".  Note -- BB should not intersect triangle");
  EXPECT_TRUE(primal::intersect(xyTri, bbOrigin));

  BoundingBoxType bbOrigin2(PointType::zero());
  bbOrigin.addPoint(PointType(-1.));
  bbOrigin.addPoint(PointType::make_point(-1., 1., 1.));
  SLIC_INFO("Testing bounding box: "
            << bbOrigin2 << " against triangle " << xyTri
            << ".  Note -- BB should not intersect triangle");
  EXPECT_TRUE(primal::intersect(xyTri, bbOrigin2));

  BoundingBoxType bbAbove(PointType::ones());
  bbAbove.addPoint(PointType(2.));
  SLIC_INFO("Testing bounding box: "
            << bbAbove << " against triangle " << xyTri
            << ".  Note -- BB should not intersect triangle");
  EXPECT_FALSE(primal::intersect(xyTri, bbAbove));

  BoundingBoxType bbBelow;
  bbBelow.addPoint(PointType(-1.));
  bbBelow.addPoint(PointType(-2.));
  SLIC_INFO("Testing bounding box: "
            << bbBelow << " against triangle " << xyTri
            << ".  Note -- BB should not intersect triangle");
  EXPECT_FALSE(primal::intersect(xyTri, bbBelow));

  BoundingBoxType bbPoint_OnTri;
  bbPoint_OnTri.addPoint(PointType::make_point(0., 1., 0.));
  SLIC_INFO("Testing point bounding box: "
            << bbPoint_OnTri << " against triangle " << xyTri
            << ".  Note -- BB is a point on triangle");
  EXPECT_TRUE(primal::intersect(xyTri, bbPoint_OnTri));

  BoundingBoxType bbPoint_OutsideTri;
  bbPoint_OutsideTri.addPoint(PointType::make_point(1., 1., 1.));
  SLIC_INFO("Testing point bounding box: "
            << bbPoint_OutsideTri << " against triangle " << xyTri
            << ".  Note -- BB is a point outside triangle");
  EXPECT_FALSE(primal::intersect(xyTri, bbPoint_OutsideTri));

  BoundingBoxType bbInvalid;
  SLIC_INFO("Testing point bounding box: "
            << bbInvalid << " against triangle " << xyTri
            << ".  Note -- BB is invalid (empty)");
  EXPECT_FALSE(primal::intersect(xyTri, bbInvalid));
}

TEST(primal_intersect, triangle_aabb_intersection_fromData)
{
  const int DIM = 3;
  typedef primal::Point<double, DIM> PointType;
  typedef primal::Triangle<double, DIM> TriangleType;
  typedef primal::BoundingBox<double, DIM> BoundingBoxType;

  PointType v0 = PointType::make_point(-31.015, 63.7756, 55.0043);
  PointType v1 = PointType::make_point(-29.0086, 59.2982, 58.0078);
  PointType v2 = PointType::make_point(-29.2009, 70.1039, 61.3229);

  TriangleType tri(v0, v1, v2);

  BoundingBoxType box0(PointType::make_point(-39.2793, 46.3735, 53.3791),
                       PointType::make_point(-26.1692, 60.1549, 57.0148));

  BoundingBoxType box1(PointType::make_point(-39.2793, 60.1549, 53.3791),
                       PointType::make_point(-26.1692, 73.9362, 57.0148));

  BoundingBoxType box2(PointType::make_point(-39.2793, 46.3735, 57.0148),
                       PointType::make_point(-26.1692, 60.1549, 60.6506));

  BoundingBoxType box3(PointType::make_point(-39.2793, 60.1549, 57.0148),
                       PointType::make_point(-26.1692, 73.9362, 60.6506));

  BoundingBoxType box4(PointType::make_point(-39.2793, 46.3735, 60.6506),
                       PointType::make_point(-26.1692, 60.1549, 64.2863));

  BoundingBoxType box5(PointType::make_point(-39.2793, 60.1549, 60.6506),
                       PointType::make_point(-26.1692, 73.9362, 64.2863));

  SLIC_INFO("Testing point bounding box: " << box0 << " against triangle "
                                           << tri);
  EXPECT_FALSE(primal::intersect(tri, box0));

  SLIC_INFO("Testing point bounding box: " << box1 << " against triangle "
                                           << tri);
  EXPECT_TRUE(primal::intersect(tri, box1));

  //
  axom::slic::setLoggingMsgLevel(axom::slic::message::Debug);

  SLIC_INFO("Testing point bounding box: " << box2 << " against triangle "
                                           << tri);
  EXPECT_TRUE(primal::intersect(tri, box2));

  axom::slic::setLoggingMsgLevel(axom::slic::message::Warning);

  SLIC_INFO("Testing point bounding box: " << box3 << " against triangle "
                                           << tri);
  EXPECT_TRUE(primal::intersect(tri, box3));

  SLIC_INFO("Testing point bounding box: " << box4 << " against triangle "
                                           << tri);
  EXPECT_FALSE(primal::intersect(tri, box4));

  SLIC_INFO("Testing point bounding box: " << box5 << " against triangle "
                                           << tri);
  EXPECT_TRUE(primal::intersect(tri, box5));
}

TEST(primal_intersect, triangle_aabb_intersection_fromData2)
{
  const int DIM = 3;
  typedef primal::Point<double, DIM> PointType;
  typedef primal::Triangle<double, DIM> TriangleType;
  typedef primal::BoundingBox<double, DIM> BoundingBoxType;

  // Triangle 569
  TriangleType tri(PointType::make_point(0, 5, 0),
                   PointType::make_point(-0.665356, 4.93844, -0.411212),
                   PointType::make_point(-0.665356, 4.93844, 0.411212));

  // {pt: (8,15,8); level: 4}
  BoundingBoxType box0(PointType::make_point(0, 4.375, 0),
                       PointType::make_point(0.625, 5, 0.625));

  // {pt: (6,15,7); level: 4}
  BoundingBoxType box1(PointType::make_point(-1.25, 4.375, -0.625),
                       PointType::make_point(-0.625, 5, 0));

  // {pt: (6,15,8); level: 4}
  BoundingBoxType box2(PointType::make_point(-1.25, 4.375, 0),
                       PointType::make_point(-0.625, 5, 0.625));

  // Block index {pt: (16,31,16); level: 5}
  BoundingBoxType box3(PointType::make_point(0, 4.6875, 0),
                       PointType::make_point(0.3125, 5, 0.3125));

  // Block index {pt: (8,15,8); level: 4}
  BoundingBoxType box4(PointType::make_point(0, 4.375, 0),
                       PointType::make_point(0.625, 5, 0.625));

  axom::slic::setLoggingMsgLevel(axom::slic::message::Info);

  SLIC_INFO("Testing point bounding box: " << box0 << " against triangle "
                                           << tri);
  EXPECT_TRUE(primal::intersect(tri, box0));

  SLIC_INFO("Testing point bounding box: " << box1 << " against triangle "
                                           << tri);
  EXPECT_TRUE(primal::intersect(tri, box1));

  SLIC_INFO("Testing point bounding box: " << box2 << " against triangle "
                                           << tri);
  EXPECT_TRUE(primal::intersect(tri, box2));

  SLIC_INFO("Testing point bounding box: " << box3 << " against triangle "
                                           << tri);
  EXPECT_TRUE(primal::intersect(tri, box3));

  SLIC_INFO("Testing point bounding box: " << box4 << " against triangle "
                                           << tri);
  EXPECT_TRUE(primal::intersect(tri, box4));

  axom::slic::setLoggingMsgLevel(axom::slic::message::Warning);
}

TEST(primal_intersect, 2D_triangle_triangle_intersection_barycentric)
{
  axom::slic::setLoggingMsgLevel(axom::slic::message::Info);

  using Triangle2 = primal::Triangle<double, 2>;
  using Point2 = primal::Point<double, 2>;
  using Bary = primal::Point<double, 3>;

  // Test several 2D triangle-triangle intersection cases
  // All cases are expected to intersect since one point is inside the triangle
  const bool expectIntersect = true;

  Triangle2 triA(Point2::make_point(0, 0),
                 Point2::make_point(1, 0),
                 Point2::make_point(1, 1));

  // Set first point to center of triA
  Point2 p0 = triA.baryToPhysical(Bary::make_point(1. / 3., 1. / 3., 1. / 3.));

  // Create some points based on barycentric coords
  std::vector<Bary> bary;

  // Add some barycentric coordinates w.r.t. the input triangle
  for(double x : {-0.2, -0.1, 0.0, 0.1, 0.2})
  {
    for(int j = -6; j <= 6; ++j)
    {
      Bary b;
      b[0] = x;
      b[1] = (1. - j / 3.) - x / 2.;
      b[2] = j / 3. - x / 2.;

      EXPECT_NEAR(1., b[0] + b[1] + b[2], 1E-8);  // check that they sum to one
      bary.push_back(b);
    }
  }

  bool includeBdry = true;

  int sz = bary.size();
  for(int i = 0; i < sz - 1; ++i)
  {
    for(int j = i; j < sz; ++j)
    {
      Triangle2 triB(p0,
                     triA.baryToPhysical(bary[i]),
                     triA.baryToPhysical(bary[j]));

      if(!triB.degenerate())
      {
        std::string str =
          fmt::format("Tri2D-Tri2D from barycenters. b1:{}, b2:{}",
                      bary[i],
                      bary[j]);
        permuteCornersTest(triA, triB, str, !includeBdry, expectIntersect);
        permuteCornersTest(triA, triB, str, includeBdry, expectIntersect);
      }
    }
  }

  axom::slic::setLoggingMsgLevel(axom::slic::message::Warning);
}

TEST(primal_intersect, 2D_triangle_triangle_intersection)
{
  typedef primal::Triangle<double, 2> Triangle2;
  typedef primal::Point<double, 2> Point2;

  // Triangle 569
  Triangle2 triA(Point2::make_point(0.0, 5.0),
                 Point2::make_point(5.0, 5.0),
                 Point2::make_point(0.0, 0.0));

  Triangle2 triB(Point2::make_point(0.0, 5.0),
                 Point2::make_point(5.0, 5.0),
                 Point2::make_point(0.0, 0.0));

  // axom::slic::setLoggingMsgLevel( axom::slic::message::Info);

  // Several intersection cases (and one non-intersection)

  permuteCornersTest(triA, triB, "identical 2D triangles", true, true);
  permuteCornersTest(triA, triB, "identical 2D triangles", false, true);

  Triangle2 triC(Point2::make_point(-1.0, -1.0),
                 Point2::make_point(-5.0, -5.0),
                 Point2::make_point(-7.0, -8.0));

  permuteCornersTest(triA, triC, "non-intersecting 2D triangles", true, false);
  permuteCornersTest(triA, triC, "non-intersecting 2D triangles", false, false);

  triA = Triangle2(Point2::make_point(4.3, 4.05),
                   Point2::make_point(-1.0, -0.06),
                   Point2::make_point(7.3, -1.3));

  triB = Triangle2(Point2::make_point(1.0, 0.0),
                   Point2::make_point(6.0, 0.5),
                   Point2::make_point(4.2, 2.1));

  permuteCornersTest(triA,
                     triB,
                     "2D tri B completely contained in tri A",
                     true,
                     true);
  permuteCornersTest(triA,
                     triB,
                     "2D tri B completely contained in tri A",
                     false,
                     true);

  triB = Triangle2(Point2::make_point(1.9, -2),
                   Point2::make_point(6.9, 2.1),
                   Point2::make_point(0.8, 5.1));

  permuteCornersTest(triA,
                     triB,
                     "intersecting 2D triangles, no corner in",
                     true,
                     true);
  permuteCornersTest(triA,
                     triB,
                     "intersecting 2D triangles, no corner in",
                     false,
                     true);

  triB = Triangle2(Point2::make_point(2.9, 1.6),
                   Point2::make_point(-1.5, 1.5),
                   Point2::make_point(0.8, 5.1));

  permuteCornersTest(triA,
                     triB,
                     "intersecting 2D triangles, one corner in",
                     true,
                     true);
  permuteCornersTest(triA,
                     triB,
                     "intersecting 2D triangles, one corner in",
                     false,
                     true);

  triB = Triangle2(Point2::make_point(2.9, 0),
                   Point2::make_point(2.1, 0.1),
                   Point2::make_point(0.8, 5.1));

  permuteCornersTest(triA,
                     triB,
                     "intersecting 2D triangles, two corners in",
                     true,
                     true);
  permuteCornersTest(triA,
                     triB,
                     "intersecting 2D triangles, two corners in",
                     false,
                     true);

  triB = Triangle2(Point2::make_point(2, -1),
                   Point2::make_point(-1.0, -0.06),
                   Point2::make_point(7.3, -1.3));

  permuteCornersTest(triA,
                     triB,
                     "2D t1 and t2 share a complete edge (and nothing else)",
                     true,
                     true);
  permuteCornersTest(triA,
                     triB,
                     "2D t1 and t2 share a complete edge (and nothing else)",
                     false,
                     false);

  Triangle2 triD(Point2::make_point(0, 0),
                 Point2::make_point(1, 0),
                 Point2::make_point(1, 1));

  Triangle2 triE(Point2::make_point(0, 0),
                 Point2::make_point(0.5, 0),
                 Point2::make_point(-1, -1));

  permuteCornersTest(triD,
                     triE,
                     "2D t1 edge is a subset of t2's, and they share a corner "
                     "(but nothing else)",
                     true,
                     true);
  permuteCornersTest(triD,
                     triE,
                     "2D t1 edge is a subset of t2's, and they share a corner "
                     "(but nothing else)",
                     false,
                     false);

  triE = Triangle2(Point2::make_point(0.5, 0),
                   Point2::make_point(1, 0),
                   Point2::make_point(-1, -1));

  permuteCornersTest(triD,
                     triE,
                     "2D t1 edge is a subset of t2's, and they share the other "
                     "corner (but nothing else)",
                     true,
                     true);
  permuteCornersTest(triD,
                     triE,
                     "2D t1 edge is a subset of t2's, and they share the other "
                     "corner (but nothing else)",
                     false,
                     false);

  triE = Triangle2(Point2::make_point(0.5, 0),
                   Point2::make_point(1.5, 0),
                   Point2::make_point(-1, -1));

  permuteCornersTest(triD,
                     triE,
                     "2D t1 edge overlaps t2 (no other intersection)",
                     true,
                     true);
  permuteCornersTest(triD,
                     triE,
                     "2D t1 edge overlaps t2 (no other intersection)",
                     false,
                     false);

  triE = Triangle2(Point2::make_point(-0.5, 0),
                   Point2::make_point(0.5, 0),
                   Point2::make_point(-1, -1));

  permuteCornersTest(
    triD,
    triE,
    "2D t1 edge overlaps t2 the other way (no other intersection)",
    true,
    true);
  permuteCornersTest(
    triD,
    triE,
    "2D t1 edge overlaps t2 the other way (no other intersection)",
    false,
    false);

  triE = Triangle2(Point2::make_point(-1, 0.5),
                   Point2::make_point(-1, -1),
                   Point2::make_point(2, -1));

  permuteCornersTest(triD,
                     triE,
                     "2D t1 point lands on t2 edge (no other intersection)",
                     true,
                     true);
  permuteCornersTest(triD,
                     triE,
                     "2D t1 point lands on t2 edge (no other intersection)",
                     false,
                     false);

  triE = Triangle2(Point2::make_point(0, 0),
                   Point2::make_point(-40, -0.7),
                   Point2::make_point(-23, 1.3));

  permuteCornersTest(triD,
                     triE,
                     "2D t1 point lands on t2 point (no other intersection)",
                     true,
                     true);
  permuteCornersTest(triD,
                     triE,
                     "2D t1 point lands on t2 point (no other intersection)",
                     false,
                     false);

  // Several non-intersection cases (and a few intersection)

  triE = Triangle2(Point2::make_point(0.2, -1e-3),
                   Point2::make_point(1, -1),
                   Point2::make_point(1.2, -1e-3));

  permuteCornersTest(triD, triE, "2D disjunct, close parallel sides", true, false);
  permuteCornersTest(triD, triE, "2D disjunct, close parallel sides", false, false);

  triE = Triangle2(Point2::make_point(0.2, -1e-3),
                   Point2::make_point(1, -1),
                   Point2::make_point(1, -1e-4));

  permuteCornersTest(triD, triE, "2D disjunct, close converging sides", true, false);
  permuteCornersTest(triD, triE, "2D disjunct, close converging sides", false, false);

  triE = Triangle2(Point2::make_point(10, 1),
                   Point2::make_point(2, 0),
                   Point2::make_point(11, -0.3));

  permuteCornersTest(triD, triE, "2D disjunct, fairly far-separated", true, false);
  permuteCornersTest(triD, triE, "2D disjunct, fairly far-separated", false, false);

  triE = Triangle2(Point2::make_point(0, 0.1),
                   Point2::make_point(-40, -0.7),
                   Point2::make_point(-23, 1.3));

  permuteCornersTest(triD, triE, "2D disjunct, point comes close", true, false);
  permuteCornersTest(triD, triE, "2D disjunct, point comes close", false, false);

  triE = Triangle2(Point2::make_point(-0.001, 0),
                   Point2::make_point(-40, -0.7),
                   Point2::make_point(-23, 1.3));

  permuteCornersTest(triD, triE, "2D disjunct, point comes close 2", true, false);
  permuteCornersTest(triD, triE, "2D disjunct, point comes close 2", false, false);

  triE = Triangle2(Point2::make_point(-0.5, 0),
                   Point2::make_point(-40, -0.7),
                   Point2::make_point(-23, 1.3));

  permuteCornersTest(triD, triE, "2D disjunct, point comes close 3", true, false);
  permuteCornersTest(triD, triE, "2D disjunct, point comes close 3", false, false);

  triE = Triangle2(Point2::make_point(-1.7, 0),
                   Point2::make_point(-40, -0.7),
                   Point2::make_point(-23, 1.3));

  permuteCornersTest(triD, triE, "2D disjunct, point comes close 4", true, false);
  permuteCornersTest(triD, triE, "2D disjunct, point comes close 4", false, false);

  triE = Triangle2(Point2::make_point(-5.1, 0),
                   Point2::make_point(-40, -0.7),
                   Point2::make_point(-23, 1.3));

  permuteCornersTest(triD, triE, "2D disjunct, point comes close 5", true, false);
  permuteCornersTest(triD, triE, "2D disjunct, point comes close 5", false, false);

  triE = Triangle2(Point2::make_point(0.5, 0.5),
                   Point2::make_point(-40, -0.7),
                   Point2::make_point(-23, 1.3));

  permuteCornersTest(triD, triE, "2D point lands on side 2", true, true);
  permuteCornersTest(triD, triE, "2D point lands on side 2", false, false);

  triE = Triangle2(Point2::make_point(0.49999, 0.5),
                   Point2::make_point(-40, -0.7),
                   Point2::make_point(-23, 1.3));

  permuteCornersTest(triD, triE, "2D point comes close to side", true, false);
  permuteCornersTest(triD, triE, "2D point comes close to side", false, false);

  triE = Triangle2(Point2::make_point(0.49, 0.5),
                   Point2::make_point(-40, -0.7),
                   Point2::make_point(-23, 1.3));

  permuteCornersTest(triD, triE, "2D point comes close to side 2", true, false);
  permuteCornersTest(triD, triE, "2D point comes close to side 2", false, false);

  triE = Triangle2(Point2::make_point(0.4, 0.5),
                   Point2::make_point(-40, -0.7),
                   Point2::make_point(-23, 1.3));

  permuteCornersTest(triD, triE, "2D point comes close to side 3", true, false);
  permuteCornersTest(triD, triE, "2D point comes close to side 3", false, false);

  triE = Triangle2(Point2::make_point(-0.1, 0.5),
                   Point2::make_point(-40, -0.7),
                   Point2::make_point(-23, 1.3));

  permuteCornersTest(triD, triE, "2D point comes close to side 4", true, false);
  permuteCornersTest(triD, triE, "2D point comes close to side 4", false, false);

  triE = Triangle2(Point2::make_point(-2.6, 2.5),
                   Point2::make_point(-40, -0.7),
                   Point2::make_point(-23, 1.3));

  permuteCornersTest(triD, triE, "2D point comes close to side 5", true, false);
  permuteCornersTest(triD, triE, "2D point comes close to side 5", false, false);

  triE = Triangle2(Point2::make_point(-6, 5),
                   Point2::make_point(-40, -0.7),
                   Point2::make_point(-23, 1.3));

  permuteCornersTest(triD, triE, "2D point comes close to side 6", true, false);
  permuteCornersTest(triD, triE, "2D point comes close to side 6", false, false);
}

bool makeTwoRandomIntersecting3DTriangles(primal::Triangle<double, 3>& l,
                                          primal::Triangle<double, 3>& r)
{
  typedef primal::Triangle<double, 3> Triangle3;
  typedef primal::Point<double, 3> Point3;
  typedef primal::Vector<double, 3> Vector3;

  //Step 1: Construct a random triangle
  Point3 A = randomPt<3>(0., 1.);
  Point3 B = randomPt<3>(0., 1.);
  Point3 C = randomPt<3>(0., 1.);
  l = Triangle3(A, B, C);

  //Step 2: Construct two random points on the triangle.
  Point3 P;
  Point3 Q;

  double a1 = randomDouble();
  double a2 = randomDouble();
  double a3 = randomDouble();

  double n1 = (a1 / (a1 + a2 + a3));
  double n2 = (a2 / (a1 + a2 + a3));
  double n3 = (a3 / (a1 + a2 + a3));

  double P_x = n1 * A[0] + n2 * B[0] + n3 * C[0];
  double P_y = n1 * A[1] + n2 * B[1] + n3 * C[1];
  double P_z = n1 * A[2] + n2 * B[2] + n3 * C[2];
  P = Point3::make_point(P_x, P_y, P_z);

  a1 = randomDouble();
  a2 = randomDouble();
  a3 = randomDouble();

  n1 = (a1 / (a1 + a2 + a3));
  n2 = (a2 / (a1 + a2 + a3));
  n3 = (a3 / (a1 + a2 + a3));

  double Q_x = n1 * A[0] + n2 * B[0] + n3 * C[0];
  double Q_y = n1 * A[1] + n2 * B[1] + n3 * C[1];
  double Q_z = n1 * A[2] + n2 * B[2] + n3 * C[2];

  Q = Point3::make_point(Q_x, Q_y, Q_z);

  /*PQ is so random segment on the triangle.  We create a vertex called vertex1
     and use it to create the triangle formed by P', Q' and vertex1. */

  //Step 3: choose some vertex away from the triangle
  Point3 vertex1 = randomPt<3>(0., 1.);

  //Step 4:
  // we scale the segments formed by both vertex 1 and P and by vertex 1 and Q
  // so that we now have a triangle whose base is not necessarily on the plane
  // formed by ABC
  Vector3 vertex2Direction = Vector3(Q, vertex1);
  Vector3 vertex3Direction = Vector3(P, vertex1);

  //construct the other two vertices of the triangle
  Point3 vertex2 = Point3::make_point(vertex1[0] - 2 * vertex2Direction[0],
                                      vertex1[1] - 2 * vertex2Direction[1],
                                      vertex1[2] - 2 * vertex2Direction[2]);
  Point3 vertex3 = Point3::make_point(vertex1[0] - 2 * vertex3Direction[0],
                                      vertex1[1] - 2 * vertex3Direction[1],
                                      vertex1[2] - 2 * vertex3Direction[2]);

  Triangle3 before = Triangle3(vertex1, Q, P);
  r = Triangle3(vertex1, vertex2, vertex3);

  return !l.degenerate() && !r.degenerate();
}

TEST(primal_intersect, 3D_triangle_triangle_intersection)
{
  typedef primal::Triangle<double, 3> Triangle3;
  typedef primal::Point<double, 3> Point3;

  Triangle3 tri3d_1(Point3::make_point(-1.0, -1.0, -1.0),
                    Point3::make_point(-2.0, -5.0, -5.0),
                    Point3::make_point(-4.0, -8.0, -8.0));

  Triangle3 tri3d_2(Point3::make_point(-1.0, -1.0, -1.0),
                    Point3::make_point(-2.0, -5.0, -5.0),
                    Point3::make_point(-4.0, -8.0, -8.0));

  permuteCornersTest(tri3d_1, tri3d_2, "3D identical triangles", true, true);
  permuteCornersTest(tri3d_1, tri3d_2, "3D identical triangles", false, true);

  Triangle3 tri3d_3(Point3::make_point(1.0, 1.0, 1.0),
                    Point3::make_point(5.0, 5.0, 5.0),
                    Point3::make_point(8.0, 7.0, 92.0));

  permuteCornersTest(tri3d_1, tri3d_3, "3D disjunct triangles", true, false);
  permuteCornersTest(tri3d_1, tri3d_3, "3D disjunct triangles", false, false);

  Triangle3 tri3A(Point3::make_point(0, 0, 0),
                  Point3::make_point(1, 0, 0),
                  Point3::make_point(0, 1.7, 2.3));

  Triangle3 tri3B(Point3::make_point(0, 0, 0),
                  Point3::make_point(1, 0, 0),
                  Point3::make_point(0, -2, 1.2));

  permuteCornersTest(tri3A, tri3B, "3D tris sharing a segment", true, true);
  permuteCornersTest(tri3A, tri3B, "3D tris sharing a segment", false, false);

  tri3B = Triangle3(Point3::make_point(-0.2, 0, 0),
                    Point3::make_point(0.7, 0, 0),
                    Point3::make_point(0, -2, 1.2));

  permuteCornersTest(tri3A, tri3B, "3D tris sharing part of a segment", true, true);
  permuteCornersTest(tri3A, tri3B, "3D tris sharing part of a segment", false, false);

  tri3B = Triangle3(Point3::make_point(-1, 0, 0),
                    Point3::make_point(0, 4.3, 6),
                    Point3::make_point(0, 1.7, 2.3));

  permuteCornersTest(tri3A, tri3B, "3D tris sharing a vertex", true, true);
  permuteCornersTest(tri3A, tri3B, "3D tris sharing a vertex", false, false);

  tri3B = Triangle3(Point3::make_point(0, -1, 0),
                    Point3::make_point(1, 1, 0),
                    Point3::make_point(0, 1.7, -2.3));

  permuteCornersTest(tri3A, tri3B, "3D tris, edges cross", true, true);
  permuteCornersTest(tri3A, tri3B, "3D tris, edges cross", false, false);

  tri3B = Triangle3(Point3::make_point(0, -1, -1),
                    Point3::make_point(0.5, 0, 0),
                    Point3::make_point(1, 1, -1));

  permuteCornersTest(tri3A, tri3B, "3D tris, B vertex lands on A's edge", true, true);
  permuteCornersTest(tri3A,
                     tri3B,
                     "3D tris, B vertex lands on A's edge",
                     false,
                     false);

  tri3B = Triangle3(Point3::make_point(0.5, -1, 0.1),
                    Point3::make_point(0.5, 1, 0.1),
                    Point3::make_point(1, 1, -1));

  permuteCornersTest(tri3A,
                     tri3B,
                     "3D tris intersect like two links in a chain",
                     true,
                     true);
  permuteCornersTest(tri3A,
                     tri3B,
                     "3D tris intersect like two links in a chain",
                     false,
                     true);

  tri3B = Triangle3(Point3::make_point(-1, -1, 1),
                    Point3::make_point(0, 2, 1),
                    Point3::make_point(5, 0, 1));

  permuteCornersTest(tri3A, tri3B, "3D tri A pokes through B", true, true);
  permuteCornersTest(tri3A, tri3B, "3D tri A pokes through B", false, true);

  tri3B = Triangle3(Point3::make_point(1, -1, 1),
                    Point3::make_point(1, 2, 1),
                    Point3::make_point(1, 0, -1));

  permuteCornersTest(tri3A, tri3B, "3D tri A vertex tangent on B", true, true);
  permuteCornersTest(tri3A, tri3B, "3D tri A vertex tangent on B", false, false);

  tri3B = Triangle3(Point3::make_point(1.00001, -1, 1),
                    Point3::make_point(1, 2, 1),
                    Point3::make_point(1, 0, -1));

  permuteCornersTest(tri3A,
                     tri3B,
                     "3D tri A vertex not quite tangent on B",
                     true,
                     false);
  permuteCornersTest(tri3A,
                     tri3B,
                     "3D tri A vertex not quite tangent on B",
                     false,
                     false);

  // 3D versions of 2D test cases (!)

  //future work: test triangle triangle with a bunch of random test cases

  srand(1);  // we want same random number sequence everytime to make sure our
             // tests don't differ on a case to case basis

  //Randomly generate a bunch of intersecting triangles (whose intersections
  // form segments) and test them How many tests are we actually performing
  // here?
  int rantests = 0;
  int skiptests = 0;
  for(int i = 0; i < 5000; i++)
  {
    Triangle3 randomTriangle, intersectingTriangle;

    if(makeTwoRandomIntersecting3DTriangles(randomTriangle, intersectingTriangle))
    {
      permuteCornersTest(randomTriangle, intersectingTriangle, "random", true, true);
      rantests += 1;
    }
    else
    {
      skiptests += 1;
    }
  }

  SLIC_INFO("Ran " << rantests << " and skipped " << skiptests
                   << " tests due to triangle degeneracy.");

  axom::slic::setLoggingMsgLevel(axom::slic::message::Warning);
}

TEST(primal_intersect, 3D_triangle_triangle_intersection_regression)
{
  axom::slic::setLoggingMsgLevel(axom::slic::message::Info);
  SLIC_INFO("Triangle intersection regression tests for discovered problems");

  typedef primal::Triangle<double, 3> Triangle3;
  typedef primal::Point<double, 3> Point3;

  {
    std::string msg = "Adjacent coplanar triangles with obtuse angles";
    Point3 vdata[4] = {Point3::make_point(-1.83697e-14, 62.5, 300),
                       Point3::make_point(16.17619, 60.37037, 300),
                       Point3::make_point(-5.790149e-16, 11.26926, 9.456031),
                       Point3::make_point(2.916699, 10.88527, 9.456031)};

    Triangle3 tri3d_1(vdata[0], vdata[1], vdata[2]);
    Triangle3 tri3d_2(vdata[2], vdata[1], vdata[3]);

    permuteCornersTest(tri3d_1, tri3d_2, msg, false, false);
    permuteCornersTest(tri3d_1, tri3d_2, msg, true, true);
  }

  {
    std::string msg = "Vertex-adjacent triangles";
    /* clang-format off */
    Point3 vdata[] = { Point3::make_point(-138.02488708496094,    -14398.0908203125, 111881.2421875),
                       Point3::make_point(   0.067092768847942352,-14407.21875,      111891.078125),
                       Point3::make_point(-136.77900695800781,    -14416.4912109375, 111891.078125),
                       Point3::make_point(   1.1611454486846924,  -14423.3466796875, 111904.359375),
                       Point3::make_point( 136.91319274902344,    -14397.947265625,  111891.078125)};
    /* clang-format on */

    Triangle3 tri3d_1(vdata[0], vdata[1], vdata[2]);
    Triangle3 tri3d_2(vdata[3], vdata[1], vdata[4]);

    permuteCornersTest(tri3d_1, tri3d_2, msg, false, false);
    permuteCornersTest(tri3d_1, tri3d_2, msg, true, true);
  }

  {
    std::string msg = "Vertex-adjacent triangles (simplified)";
    Point3 vdata[] = {Point3::make_point(-1, -1, -1),
                      Point3::make_point(0, 0, -0.005),
                      Point3::make_point(-1, 0, 0),
                      Point3::make_point(0, 1, -1),
                      Point3::make_point(1, 0, 0)};

    Triangle3 tri3d_1(vdata[0], vdata[1], vdata[2]);
    Triangle3 tri3d_2(vdata[3], vdata[1], vdata[4]);

    permuteCornersTest(tri3d_1, tri3d_2, msg, false, false);
    permuteCornersTest(tri3d_1, tri3d_2, msg, true, true);
  }

  {
    // https://github.com/LLNL/axom/issues/152
    std::string msg = "From Axom github issue #152";
    Point3 vdata[] = {Point3::make_point(76.648, 54.6752, 15.0012),
                      Point3::make_point(76.648, 54.6752, 14.5542),
                      Point3::make_point(76.582, 54.6752, 14.7879),
                      Point3::make_point(76.6252, 54.6752, 14.892),
                      Point3::make_point(76.5617, 54.6752, 14.7929)};

    Triangle3 tri3d_1(vdata[0], vdata[1], vdata[2]);
    Triangle3 tri3d_2(vdata[3], vdata[2], vdata[4]);

    // Output some debugging information about this configuration
    // (change logging level to Debug to see these messages)
    SLIC_DEBUG("Triangle 1: " << tri3d_1);
    SLIC_DEBUG("Triangle 2: " << tri3d_2);

    for(int i = 0; i < 3; ++i)
    {
      SLIC_DEBUG("t1 " << i << "\n\t-- distance " << tri3d_1[i]
                       << " to tri2: " << std::setprecision(17)
                       << sqrt(primal::squared_distance(tri3d_1[i], tri3d_2))
                       << "\n\t-- closest point: "
                       << primal::closest_point(tri3d_1[i], tri3d_2));
    }

    for(int i = 0; i < 3; ++i)
    {
      SLIC_DEBUG("t2 " << i << "\n\t-- distance " << tri3d_2[i]
                       << " to tri1: " << std::setprecision(17)
                       << sqrt(primal::squared_distance(tri3d_2[i], tri3d_1))
                       << "\n\t-- closest point: "
                       << primal::closest_point(tri3d_2[i], tri3d_1));
    }

    const bool expectIntersect = true;
    permuteCornersTest(tri3d_1, tri3d_2, msg, false, expectIntersect);
    permuteCornersTest(tri3d_1, tri3d_2, msg, true, expectIntersect);
  }

  {
    // https://github.com/LLNL/axom/issues/152
    std::string msg = "From Axom github issue #152 (simplified)";
    Point3 vdata[] = {Point3::make_point(0.066, 0, 0.2133),
                      Point3::make_point(0.066, 0, -0.2337),
                      Point3::make_point(0, 0, 0),
                      Point3::make_point(0.0432, 0, 0.1041),
                      Point3::make_point(-0.0203, 0, 0.005)};

    Triangle3 tri3d_1(vdata[0], vdata[1], vdata[2]);
    Triangle3 tri3d_2(vdata[3], vdata[2], vdata[4]);

    const bool expectIntersect = true;
    permuteCornersTest(tri3d_1, tri3d_2, msg, false, expectIntersect);
    permuteCornersTest(tri3d_1, tri3d_2, msg, true, expectIntersect);
  }

  {
    // https://github.com/LLNL/axom/issues/152
    std::string msg =
      "From Axom github issue #152 (simplified 2) -- one point inside other "
      "triangle";
    Point3 vdata[] = {Point3::make_point(1, 0, 0.5),
                      Point3::make_point(1, 0, -0.5),
                      Point3::make_point(0, 0, 0),
                      Point3::make_point(0.5, 0, 0.1),
                      Point3::make_point(-0.1, 0, -0.2)};

    Triangle3 tri3d_1(vdata[0], vdata[1], vdata[2]);
    Triangle3 tri3d_2(vdata[3], vdata[2], vdata[4]);

    permuteCornersTest(tri3d_1, tri3d_2, msg, false, true);
    permuteCornersTest(tri3d_1, tri3d_2, msg, true, true);
  }

  {
    // https://github.com/LLNL/axom/issues/152
    std::string msg =
      "From Axom github issue #152 (simplified 3) -- one point inside other "
      "triangle";
    Point3 vdata[] = {Point3::make_point(1, 0, 0.5),
                      Point3::make_point(1, 0, -0.5),
                      Point3::make_point(0, 0, 0),
                      Point3::make_point(0.5, 0, 0.1),
                      Point3::make_point(-0.1, 0, 0.05)};

    Triangle3 tri3d_1(vdata[0], vdata[1], vdata[2]);
    Triangle3 tri3d_2(vdata[3], vdata[2], vdata[4]);

    permuteCornersTest(tri3d_1, tri3d_2, msg, false, true);
    permuteCornersTest(tri3d_1, tri3d_2, msg, true, true);
  }

  {
    // https://github.com/LLNL/axom/issues/152
    std::string msg =
      "From Axom github issue #152 (simplified 4) -- one point inside other "
      "triangle";
    Point3 vdata[] = {Point3::make_point(1, 0, 0.5),
                      Point3::make_point(1, 0, -0.5),
                      Point3::make_point(0, 0, 0),
                      Point3::make_point(0.5, 0, 0.1),
                      Point3::make_point(-0.1, 0, 0.06)};

    Triangle3 tri3d_1(vdata[0], vdata[1], vdata[2]);
    Triangle3 tri3d_2(vdata[3], vdata[2], vdata[4]);

    permuteCornersTest(tri3d_1, tri3d_2, msg, false, true);
    permuteCornersTest(tri3d_1, tri3d_2, msg, true, true);
  }

  {
    // https://github.com/LLNL/axom/issues/152
    std::string msg =
      "From Axom github issue #152 (simplified 5) -- one point inside other "
      "triangle";
    Point3 vdata[] = {Point3::make_point(1, 0, 0.5),
                      Point3::make_point(1, 0, -0.5),
                      Point3::make_point(0, 0, 0),
                      Point3::make_point(0.5, 0, 0.1),
                      Point3::make_point(-0.1, 0, 0.04)};

    Triangle3 tri3d_1(vdata[0], vdata[1], vdata[2]);
    Triangle3 tri3d_2(vdata[3], vdata[2], vdata[4]);

    permuteCornersTest(tri3d_1, tri3d_2, msg, false, true);
    permuteCornersTest(tri3d_1, tri3d_2, msg, true, true);
  }

  // Revert logging level to Warning
  axom::slic::setLoggingMsgLevel(axom::slic::message::Warning);
}

TEST(primal_intersect, triangle_aabb_intersection_boundaryFace)
{
  const int DIM = 3;
  typedef primal::Point<double, DIM> PointType;
  typedef primal::Triangle<double, DIM> TriangleType;
  typedef primal::BoundingBox<double, DIM> BoundingBoxType;

  TriangleType tri(PointType::make_point(0, 5, 0),
                   PointType::make_point(0, 5, 5),
                   PointType::make_point(0, 5, 5));

  BoundingBoxType box0(PointType::make_point(-10, -10, -10),
                       PointType::make_point(0, 10, 10));

  BoundingBoxType box1(PointType::make_point(0, -10, -10),
                       PointType::make_point(10, 10, 10));

  axom::slic::setLoggingMsgLevel(axom::slic::message::Debug);

  SLIC_INFO("Testing point bounding box: " << box0 << " against triangle "
                                           << tri);
  EXPECT_TRUE(primal::intersect(tri, box0));

  SLIC_INFO("Testing point bounding box: " << box1 << " against triangle "
                                           << tri);
  EXPECT_TRUE(primal::intersect(tri, box1));

  // ---

  // Airfoil triangle 206
  TriangleType tri2(PointType::make_point(0.0340691, -1, 0.0236411),
                    PointType::make_point(0.028589, -1, 0.0221062),
                    PointType::make_point(0.0207793, -1, -0.0295674));

  // Block: (134,128,310) @ level 9
  BoundingBoxType box2(PointType::make_point(0.0230077, -1, -0.0208459),
                       PointType::make_point(0.0268708, -0.992188, -0.0201394));

  SLIC_INFO("Testing point bounding box: "
            << box2 << " against triangle " << tri2 << "\n\t -- intersects? "
            << (primal::intersect(tri2, box2) ? "yes" : "no")
            //<< "\n\t -- distance: " << (primal::distance(tri2, box2) ?
            // "yes":"no")
  );
  //EXPECT_TRUE( primal::intersect(tri, box1));

  axom::slic::setLoggingMsgLevel(axom::slic::message::Warning);
}

TEST(primal_intersect, ray_aabb_intersection_general3D)
{
  const int DIM = 3;
  typedef primal::Point<double, DIM> PointType;
  typedef primal::Ray<double, DIM> RayType;
  typedef primal::BoundingBox<double, DIM> BoundingBoxType;
  typedef primal::Vector<double, DIM> VectorType;

  // STEP 1: construct ray
  PointType origin = PointType::make_point(0.0, 0.0, 0.0);
  VectorType direction;
  direction[0] = 1.0;
  direction[1] = 1.0;
  direction[2] = 1.0;
  RayType R(origin, direction);

  BoundingBoxType box0(PointType::make_point(5.0, 5.0, 5.0),
                       PointType::make_point(10.0, 10.0, 10.0));

  BoundingBoxType box1(PointType::make_point(-5.0, -5.0, -5.0),
                       PointType::make_point(-1.0, -1.0, -1.0));

  axom::slic::setLoggingMsgLevel(axom::slic::message::Debug);
  PointType ip;

  bool intersects = primal::intersect(R, box0, ip);
  SLIC_INFO("Testing point bounding box: " << box0 << " against ray " << R);
  SLIC_INFO("Point at: " << ip);
  EXPECT_TRUE(intersects);

  intersects = primal::intersect(R, box1, ip);
  SLIC_INFO("Testing point bounding box: " << box1 << " against ray " << R);
  SLIC_INFO("Point at: " << ip);
  EXPECT_FALSE(intersects);

  //axom::slic::setLoggingMsgLevel( axom::slic::message::Warning);
}

TEST(primal_intersect, ray_aabb_intersection_tinyDirectionVector3D)
{
  const int DIM = 3;
  typedef primal::Point<double, DIM> PointType;
  typedef primal::Ray<double, DIM> RayType;
  typedef primal::BoundingBox<double, DIM> BoundingBoxType;
  typedef primal::Vector<double, DIM> VectorType;

  // STEP 1: construct ray
  PointType origin = PointType::make_point(11.0, 11.0, 11.0);
  VectorType direction;
  direction[0] = 0.0;
  direction[1] = 0.0;
  direction[2] = 0.0;
  RayType R(origin, direction);

  BoundingBoxType box0(PointType::make_point(5.0, 5.0, 5.0),
                       PointType::make_point(10.0, 10.0, 10.0));

  BoundingBoxType box1(PointType::make_point(-5.0, -5.0, -5.0),
                       PointType::make_point(-1.0, -1.0, -1.0));

  axom::slic::setLoggingMsgLevel(axom::slic::message::Debug);
  PointType ip;

  bool intersects = primal::intersect(R, box0, ip);
  SLIC_INFO("Testing point bounding box: " << box0 << " against ray " << R);
  SLIC_INFO("Point at: " << ip);
  EXPECT_FALSE(intersects);

  intersects = primal::intersect(R, box1, ip);
  SLIC_INFO("Testing point bounding box: " << box1 << " against ray " << R);
  SLIC_INFO("Point at: " << ip);
  EXPECT_FALSE(intersects);

  //axom::slic::setLoggingMsgLevel( axom::slic::message::Warning);
}

template <typename T, int DIM>
bool testPointsClose(const primal::Point<T, DIM>& lhp,
                     const primal::Point<T, DIM>& rhp,
                     double EPS = 1e-6)
{
  primal::Vector<T, DIM> v(lhp, rhp);
  return v.norm() < EPS;
}

template <typename T, int DIM>
void ensureTriPointMatchesSegPoint(const primal::Triangle<T, DIM>& tri,
                                   const primal::Point<T, 3>& bary,
                                   const primal::Segment<T, DIM>& seg,
                                   double t)
{
  primal::Point<T, DIM> tripoint = tri.at(bary);
  primal::Point<T, DIM> segpoint = seg.at(t);
  EXPECT_TRUE(testPointsClose(tripoint, segpoint));
}

template <typename T, int DIM>
void ensureTriPointMatchesRayPoint(const primal::Triangle<T, DIM>& tri,
                                   const primal::Point<T, 3>& bary,
                                   const primal::Ray<T, DIM>& ray,
                                   double t)
{
  primal::Point<T, DIM> tripoint = tri.at(bary);
  primal::Point<T, DIM> raypoint = ray.at(t);
  EXPECT_TRUE(testPointsClose(tripoint, raypoint));
}

template <int DIM>
void testRayIntersection(const primal::Triangle<double, DIM>& tri,
                         const primal::Ray<double, DIM>& ray,
                         const std::string& whattest,
                         const bool testtrue)
{
  SCOPED_TRACE(whattest);

  typedef primal::Point<double, DIM> PointType;

  PointType tip;
  double t = 0.;

  if(testtrue)
  {
    EXPECT_TRUE(intersect(tri, ray, t, tip));
    PointType tripoint = tri.baryToPhysical(tip);
    PointType raypoint = ray.at(t);
    EXPECT_TRUE(testPointsClose(tripoint, raypoint))
      << "Tripoint is " << tripoint << " and raypoint is " << raypoint;
  }
  else
  {
    EXPECT_FALSE(intersect(tri, ray, t, tip))
      << "Expected no intersection; Found one at point " << ray.at(t);
  }
}

template <int DIM>
void testTriSegBothEnds(const primal::Triangle<double, DIM>& tri,
                        const primal::Point<double, DIM>& p1,
                        const primal::Point<double, DIM>& p2,
                        const std::string& whattest,
                        const bool testtrue)
{
  SCOPED_TRACE(whattest);

  typedef primal::Point<double, DIM> PointType;
  typedef primal::Segment<double, DIM> SegmentType;

  PointType tip;
  double t = 0.;

  SegmentType seg1(p1, p2);
  SegmentType seg2(p2, p1);
  if(testtrue)
  {
    // Find the intersection of segment from p1 to p2
    double t1 = 0;
    PointType tip1;
    EXPECT_TRUE(intersect(tri, seg1, t1, tip1));
    PointType tripoint1 = tri.baryToPhysical(tip1);
    PointType segpoint1 = seg1.at(t1);
    EXPECT_TRUE(testPointsClose(tripoint1, segpoint1))
      << "Tripoint is " << tripoint1 << " and segpoint is " << segpoint1;

    // Find the intersection of segment from p1 to p2
    double t2 = 0;
    PointType tip2;
    EXPECT_TRUE(intersect(tri, seg2, t2, tip2));
    PointType tripoint2 = tri.baryToPhysical(tip2);
    PointType segpoint2 = seg2.at(t2);
    EXPECT_TRUE(testPointsClose(tripoint2, segpoint2));
  }
  else
  {
    EXPECT_FALSE(intersect(tri, seg1, t, tip))
      << "Expected no intersection; Found one at point " << seg1.at(t);

    EXPECT_FALSE(intersect(tri, seg2, t, tip))
      << "Expected no intersection; Found one at point " << seg2.at(t);
  }
}

TEST(primal_intersect, segment_aabb_2d_intersection)
{
  constexpr int DIM = 2;
  using PointType = primal::Point<double, DIM>;
  using SegmentType = primal::Segment<double, DIM>;
  using BoxType = primal::BoundingBox<double, DIM>;

  // Helper lambda for printing out intersection details
  auto print_details = [=](bool expected,
                           const BoxType& b,
                           const SegmentType& s,
                           const PointType& p) {
    if(expected)
    {
      SLIC_INFO("Found intersection between box " << b << " and line segment "
                                                  << s << " at point " << p);
    }
    else
    {
      SLIC_INFO("No expected intersection between box "
                << b << " and line segment " << s);
    }
  };

  // Simple intersection
  {
    BoxType box(PointType {3, 3}, PointType {5, 5});
    SegmentType seg(PointType {3.25, 0}, PointType {4.75, 6});
    PointType ipt;

    EXPECT_TRUE(primal::intersect(seg, box, ipt));
    print_details(true, box, seg, ipt);
  }

  // Reverse of first case
  {
    BoxType box(PointType {3, 3}, PointType {5, 5});
    SegmentType seg(PointType {4.75, 6}, PointType {3.25, 0});
    PointType ipt;

    EXPECT_TRUE(primal::intersect(seg, box, ipt));
    print_details(true, box, seg, ipt);
  }

  // No intersection
  {
    BoxType box(PointType {3, 3}, PointType {5, 5});
    SegmentType seg(PointType {-3.25, 0}, PointType {-4.75, 6});
    PointType ipt;

    EXPECT_FALSE(primal::intersect(seg, box, ipt));
    print_details(false, box, seg, ipt);
  }

  // Grazing intersection
  {
    BoxType box(PointType {3, 3}, PointType {5, 5});
    SegmentType seg(PointType {3, 4}, PointType {3, 6});
    PointType ipt;

    EXPECT_TRUE(primal::intersect(seg, box, ipt));
    print_details(true, box, seg, ipt);
  }

  // Find parametric coordinates of double intersection
  {
    const double EPS = 1e-10;

    BoxType box(PointType {3, 3}, PointType {5, 5});
    SegmentType seg(PointType {4, 2}, PointType {4, 10});
    double tmin, tmax;

    EXPECT_TRUE(primal::intersect(seg, box, tmin, tmax, EPS));
    EXPECT_NEAR(1. / 8., tmin, EPS);
    EXPECT_NEAR(3. / 8., tmax, EPS);

    PointType ipt = seg.at(tmin);
    print_details(true, box, seg, ipt);
  }

  // Find parametric coordinates of single intersection (left)
  {
    const double EPS = 1e-10;

    BoxType box(PointType {3, 3}, PointType {5, 5});
    SegmentType seg(PointType {4, 2}, PointType {4, 4});
    double tmin, tmax;

    EXPECT_TRUE(primal::intersect(seg, box, tmin, tmax, EPS));
    EXPECT_NEAR(1. / 2., tmin, EPS);
    EXPECT_NEAR(1., tmax, EPS);

    PointType ipt = seg.at(tmin);
    print_details(true, box, seg, ipt);
  }

  // Find parametric coordinates of single intersection (right)
  {
    const double EPS = 1e-10;

    BoxType box(PointType {3, 3}, PointType {5, 5});
    SegmentType seg(PointType {4, 4}, PointType {4, 10});
    double tmin, tmax;

    EXPECT_TRUE(primal::intersect(seg, box, tmin, tmax, EPS));
    EXPECT_NEAR(0, tmin, EPS);
    EXPECT_NEAR(1. / 6., tmax, EPS);

    PointType ipt = seg.at(tmin);
    print_details(true, box, seg, ipt);
  }

  // Find parametric coordinates of complete intersection
  {
    const double EPS = 1e-10;

    BoxType box(PointType {3, 3}, PointType {5, 5});
    SegmentType seg(PointType {3.5, 3.5}, PointType {4.5, 4.5});
    double tmin, tmax;

    EXPECT_TRUE(primal::intersect(seg, box, tmin, tmax, EPS));
    EXPECT_NEAR(0, tmin, EPS);
    EXPECT_NEAR(1., tmax, EPS);

    PointType ipt = seg.at(tmin);
    print_details(true, box, seg, ipt);
  }
}

TEST(primal_intersect, segment_aabb_3d_intersection)
{
  constexpr int DIM = 3;
  using PointType = primal::Point<double, DIM>;
  using SegmentType = primal::Segment<double, DIM>;
  using BoxType = primal::BoundingBox<double, DIM>;

  // Helper lambda for printing out intersection details
  auto print_details = [=](bool expected,
                           const BoxType& b,
                           const SegmentType& s,
                           const PointType& p) {
    if(expected)
    {
      SLIC_INFO("Found intersection between box " << b << " and line segment "
                                                  << s << " at point " << p);
    }
    else
    {
      SLIC_INFO("No expected intersection between box "
                << b << " and line segment " << s);
    }
  };

  // Simple intersection
  {
    BoxType box(PointType {3, 3, 3}, PointType {5, 5, 5});
    SegmentType seg(PointType {3.25, 0, 4}, PointType {4.75, 6, 4});
    PointType ipt;

    EXPECT_TRUE(primal::intersect(seg, box, ipt));
    print_details(true, box, seg, ipt);
  }

  // Reverse of first case
  {
    BoxType box(PointType {3, 3, 3}, PointType {5, 5, 5});
    SegmentType seg(PointType {4.75, 6, 4}, PointType {3.25, 0, 4});
    PointType ipt;

    EXPECT_TRUE(primal::intersect(seg, box, ipt));
    print_details(true, box, seg, ipt);
  }

  // No intersection
  {
    BoxType box(PointType {3, 3, 3}, PointType {5, 5, 5});
    SegmentType seg(PointType {-3.25, 0, 2}, PointType {-4.75, 6, -2});
    PointType ipt;

    EXPECT_FALSE(primal::intersect(seg, box, ipt));
    print_details(false, box, seg, ipt);
  }

  // Grazing intersection
  {
    BoxType box(PointType {3, 3, 3}, PointType {5, 5, 5});
    SegmentType seg(PointType {3, 4, 3}, PointType {3, 6, 3});
    PointType ipt;

    EXPECT_TRUE(primal::intersect(seg, box, ipt));
    print_details(true, box, seg, ipt);
  }

  // Find parametric coordinates of double intersection
  {
    const double EPS = 1e-10;

    BoxType box(PointType {3, 3, 3}, PointType {5, 5, 5});
    SegmentType seg(PointType {4, 2, 4}, PointType {4, 10, 4});
    double tmin, tmax;

    EXPECT_TRUE(primal::intersect(seg, box, tmin, tmax, EPS));
    EXPECT_NEAR(1. / 8., tmin, EPS);
    EXPECT_NEAR(3. / 8., tmax, EPS);

    PointType ipt = seg.at(tmin);
    print_details(true, box, seg, ipt);
  }
}

TEST(primal_intersect, triangle_segment_intersection)
{
  const int DIM = 3;
  typedef primal::Point<double, DIM> PointType;
  typedef primal::Triangle<double, DIM> TriangleType;
  typedef primal::Segment<double, DIM> SegmentType;

  double xArr[3] = {1., 0., 0.};
  double yArr[3] = {0., 1., 0.};
  double zArr[3] = {0., 0., 1.};
  double mArr[3] = {1. / 3., 1. / 3., 1. / 3.};

  PointType ptX(xArr);
  PointType ptY(yArr);
  PointType ptZ(zArr);
  PointType ptM(mArr);
  PointType r0 = PointType::make_point(5., 5., 5.);
  PointType testp = PointType::make_point(6., 5., 5.);

  TriangleType tri(ptX, ptY, ptZ);
  SegmentType testSeg(r0, ptX);

  // Clear miss
  testTriSegBothEnds(tri, r0, testp, "clear miss", false);

  // Succession of misses
  // Copied from ray test, and testing both orders for segment (AB and BA)
  testTriSegBothEnds(tri, r0, testp, "miss 1", false);
  testp = PointType::make_point(0., .5, .6);
  testTriSegBothEnds(tri, r0, testp, "miss 2", false);
  testp = PointType::make_point(0., .85, .16);
  testTriSegBothEnds(tri, r0, testp, "miss 3", false);
  testp = PointType::make_point(.4, 1.2, 0);
  testTriSegBothEnds(tri, r0, testp, "miss 4", false);
  testp = PointType::make_point(1., 0.000001, 0);
  testTriSegBothEnds(tri, r0, testp, "miss 5", false);
  testp = PointType::make_point(0.4, 0, 0.7);
  testTriSegBothEnds(tri, r0, testp, "miss 6", false);
  testp = PointType::make_point(0.3, 0.4, 0.5);
  testTriSegBothEnds(tri, r0, testp, "miss 7", false);
  testp = PointType::make_point(0.4, 0.4, 0.4);
  testTriSegBothEnds(tri, r0, testp, "miss 8", false);

  // Some hits
  testp = PointType::make_point(0.78, -0.2, -0.2);
  testTriSegBothEnds(tri, r0, testp, "hit 1", true);
  testp = PointType::make_point(0.4, 0.3, 0.2);
  testTriSegBothEnds(tri, r0, testp, "hit 2", true);
  testp = PointType::make_point(0.2, 0.2, 0.2);
  testTriSegBothEnds(tri, r0, testp, "hit 3", true);

  // End points, triangle boundaries
  PointType testp2 = PointType::make_point(1., 1., 1.);
  testp = PointType::make_point(1., .1, .1);
  testTriSegBothEnds(tri, testp, testp2, "shy of corner", false);
  testp = PointType::make_point(1., -.1, -.1);
  testTriSegBothEnds(tri, testp, testp2, "beyond corner", true);
  testTriSegBothEnds(tri, testp, ptX, "beyond corner 2", true);

  testp2 = PointType::make_point(0, 1, 1);
  testp = PointType::make_point(0, .4, .7);
  testTriSegBothEnds(tri, testp, testp2, "shy of edge", false);
  testp = PointType::make_point(0, .6, .3);
  testTriSegBothEnds(tri, testp, testp2, "beyond edge", true);
  testp = PointType::make_point(0, .7, .3);
  testTriSegBothEnds(tri, testp2, ptM, "segment endpoint 1", true);
  testTriSegBothEnds(tri, ptM, r0, "segment endpoint 2", true);
}

TEST(primal_intersect, triangle_ray_intersection)
{
  const int DIM = 3;
  typedef primal::Point<double, DIM> PointType;
  typedef primal::Triangle<double, DIM> TriangleType;
  typedef primal::Ray<double, DIM> RayType;
  typedef primal::Segment<double, DIM> SegmentType;

  double xArr[3] = {1., 0., 0.};
  double yArr[3] = {0., 1., 0.};
  double zArr[3] = {0., 0., 1.};
  double mArr[3] = {1. / 3., 1. / 3., 1. / 3.};

  double nxArr[3] = {-1., 2., 2.};
  double nyArr[3] = {2., -1., 2.};
  double nzArr[3] = {2., 2., -1.};

  PointType ptX(xArr);
  PointType ptY(yArr);
  PointType ptZ(zArr);
  PointType ptnX(nxArr);
  PointType ptnY(nyArr);
  PointType ptnZ(nzArr);
  PointType ptM(mArr);
  PointType r0 = PointType::make_point(5., 5., 5.);
  PointType o = PointType::make_point(0, 0, 0);
  PointType ox = PointType::make_point(1, 0, 0);
  PointType oy = PointType::make_point(0, 1, 0);

  TriangleType tri(ptX, ptY, ptZ);
  TriangleType tri2(o, ox, oy);
  RayType testRay(SegmentType(ptX, ptY));

  // Clear miss
  testRay = RayType(SegmentType(r0, PointType::make_point(6., 5., 5.)));
  testRayIntersection(tri, testRay, "clear miss", false);

  // More misses
  testRay = RayType(SegmentType(r0, PointType::make_point(0., 1., .6)));
  testRayIntersection(tri, testRay, "miss 1", false);
  testRay = RayType(SegmentType(r0, PointType::make_point(0., .5, .6)));
  testRayIntersection(tri, testRay, "miss 2", false);
  testRay = RayType(SegmentType(r0, PointType::make_point(0., .85, .16)));
  testRayIntersection(tri, testRay, "miss 3", false);
  testRay = RayType(SegmentType(r0, PointType::make_point(.4, 1.2, 0)));
  testRayIntersection(tri, testRay, "miss 4", false);
  testRay = RayType(SegmentType(r0, PointType::make_point(1., 0.000001, 0)));
  testRayIntersection(tri, testRay, "miss 5", false);
  testRay = RayType(SegmentType(r0, PointType::make_point(0.4, 0, 0.7)));
  testRayIntersection(tri, testRay, "miss 6", false);

  // Edge intersections should be reported as hits
  PointType tripoint;
  testRay = RayType(SegmentType(r0, ptX));
  testRayIntersection(tri, testRay, "edge hit 1", true);
  testRay = RayType(SegmentType(r0, ptY));
  testRayIntersection(tri, testRay, "edge hit 2", true);
  testRay = RayType(SegmentType(r0, ptZ));
  testRayIntersection(tri, testRay, "edge hit 3", true);
  testRay = RayType(SegmentType(r0, PointType::make_point(0., 0.7, 0.3)));
  testRayIntersection(tri, testRay, "edge hit 4", true);
  testRay = RayType(SegmentType(r0, PointType::make_point(0.7, 0.3, 0.)));
  testRayIntersection(tri, testRay, "edge hit 5", true);
  testRay = RayType(SegmentType(o, PointType::make_point(0.2, 0., 0.8)));
  testRayIntersection(tri, testRay, "edge hit 6", true);

  // Hits
  testRay = RayType(SegmentType(r0, PointType::make_point(0.2, 0., 0.2)));
  testRayIntersection(tri, testRay, "hit 1", true);
  testRay = RayType(SegmentType(r0, PointType::make_point(0., 0., 0.)));
  testRayIntersection(tri, testRay, "hit 2", true);
  testRay = RayType(SegmentType(r0, PointType::make_point(0.1, 0.6, 0.)));
  testRayIntersection(tri, testRay, "hit 3", true);

  // Coplanar miss
  testRay = RayType(SegmentType(PointType::make_point(-0.1, 1.1, 0.),
                                PointType::make_point(-0.1, 0., 1.1)));
  testRayIntersection(tri, testRay, "coplanar miss", false);

  // Coplanar intersection (reported as miss by function)
  testRay = RayType(SegmentType(PointType::make_point(1, 0.5, 0),
                                PointType::make_point(-1, 0.5, 0)));
  testRayIntersection(tri2,
                      testRay,
                      "coplanar intersection, reported as miss by design",
                      false);

  // Coplanar, interior ray origin (reported as miss by function)
  testRay = RayType(SegmentType(ptM, PointType::make_point(0.5, 0., 0.5)));
  testRayIntersection(
    tri,
    testRay,
    "coplanar interior ray origin, reported as miss by design",
    false);

  // Not coplanar, interior ray origin (reported as miss by function)
  testRay = RayType(SegmentType(PointType::make_point(0.2, 0.18, 0),
                                PointType::make_point(0., 0., 0.5)));
  testRayIntersection(
    tri2,
    testRay,
    "non-coplanar interior ray origin, reported as miss by design",
    false);
}

TEST(primal_intersect, triangle_ray_intersection_unit_ray)
{
  const int DIM = 3;
  typedef primal::Point<double, DIM> PointType;
  typedef primal::Vector<double, DIM> VectorType;
  typedef primal::Triangle<double, DIM> TriangleType;
  typedef primal::Ray<double, DIM> RayType;

  PointType o = PointType::zero();
  VectorType v = VectorType::make_vector(0, 0, 1);
  RayType r(o, v);

  TriangleType t(PointType::make_point(2, 2, 2),
                 PointType::make_point(2, 4, 2),
                 PointType::make_point(3, 3, 2));

  EXPECT_FALSE(axom::primal::intersect(t, r));

  PointType intBary;
  double intersectionParam = 0.;
  TriangleType t2(PointType::make_point(-1, -1, 2),
                  PointType::make_point(-1, 1, 2),
                  PointType::make_point(2, 0, 2));
  EXPECT_TRUE(axom::primal::intersect(t2, r, intersectionParam, intBary));

  // Here, intersectionParam is the distance along the ray, with the source
  // being 0.
  // This is different from a segment's intersection parameter (see below).
  EXPECT_DOUBLE_EQ(2.0, intersectionParam);

  PointType intersectionPoint = r.at(intersectionParam);
  PointType triIntersectionPoint = t2.baryToPhysical(intBary);
  SLIC_INFO("Intersection (unscaled barycentric) is " << intBary);
  SLIC_INFO("Intersection param is " << intersectionParam);
  SLIC_INFO("Intersection point is " << intersectionPoint);
  EXPECT_TRUE(testPointsClose(intersectionPoint, triIntersectionPoint));
}

TEST(primal_intersect, triangle_ray_intersection_fpe_regression)
{
  const int DIM = 3;
  typedef primal::Point<double, DIM> PointType;
  typedef primal::Vector<double, DIM> VectorType;
  typedef primal::Triangle<double, DIM> TriangleType;
  typedef primal::Ray<double, DIM> RayType;

  // This regression test triggers a division by zero in
  // a previous implementation of primal::intersect()
  // when normalizing the barycentric coordinates

  double t = 0.;
  PointType bary;

  double val = std::sqrt(2) / 2;
  RayType ray(PointType::make_point(-0.1, 1.1, 0.),
              VectorType::make_vector(0, -val, val));

  TriangleType tri(PointType::make_point(1, 0, 0),
                   PointType::make_point(0, 1, 0),
                   PointType::make_point(0, 0, 1));

  EXPECT_FALSE(axom::primal::intersect(tri, ray, t, bary));
}

TEST(primal_intersect, triangle_ray_intersection_unit_seg)
{
  const int DIM = 3;
  typedef primal::Point<double, DIM> PointType;
  typedef primal::Triangle<double, DIM> TriangleType;
  typedef primal::Segment<double, DIM> SegmentType;

  PointType o = PointType::zero();
  PointType d = PointType::make_point(0, 0, 4);
  SegmentType s(o, d);

  TriangleType t(PointType::make_point(-1, -1, 2),
                 PointType::make_point(-1, 1, 2),
                 PointType::make_point(2, 0, 2));

  PointType intBary;
  double intersectionParam = 0.;
  EXPECT_TRUE(axom::primal::intersect(t, s, intersectionParam, intBary));

  // Here, intersectionParam is the distance along the segment, with the source
  // being 0
  // and the target being 1.  This is different from a ray's intersection
  // parameter
  // (see above).
  EXPECT_DOUBLE_EQ(0.5, intersectionParam);

  PointType intersectionPoint = s.at(intersectionParam);
  PointType triIntersectionPoint = t.baryToPhysical(intBary);
  SLIC_INFO("Intersection (unscaled barycentric) is " << intBary);
  SLIC_INFO("Intersection param is " << intersectionParam);
  SLIC_INFO("Intersection point is " << intersectionPoint);
  EXPECT_TRUE(testPointsClose(intersectionPoint, triIntersectionPoint));
}

//------------------------------------------------------------------------------
TEST(primal_intersect, obb_obb_test_intersection2D)
{
  const int DIM = 2;
  typedef double CoordType;
  typedef primal::Point<CoordType, DIM> QPoint;
  typedef primal::Vector<CoordType, DIM> QVector;
  typedef primal::OrientedBoundingBox<CoordType, DIM> QOBBox;

  QPoint pt1;       // origin
  QVector u1[DIM];  // make standard axes
  QVector u2[DIM];  // axes rotated by 90 degrees
  for(int i = 0; i < DIM; i++)
  {
    u1[i] = QVector();
    u2[i] = QVector();
    u1[i][i] = 1.;
    u2[i][0] = 1.;
    u2[i][1] = 1. - 2 * i;
  }

  QVector e1(1.);
  QVector e2(1.42);
  QPoint pt2(2.);

  // Pass a box which just kisses on an edge initially through the unit box
  QOBBox obbox1(pt1, u1, e1);
  QOBBox obbox2(pt2, u2, e2);

  EXPECT_TRUE(axom::primal::intersect(obbox1, obbox2));

  QVector shift(-0.5);
  obbox2.shift(shift);

  EXPECT_TRUE(axom::primal::intersect(obbox1, obbox2));

  obbox2.shift(2. * shift);
  EXPECT_TRUE(axom::primal::intersect(obbox1, obbox2));

  obbox2.shift(100. * shift);
  EXPECT_FALSE(axom::primal::intersect(obbox1, obbox2));

  // check edge-edge intersection
  QOBBox obbox3;
  QOBBox obbox4;

  obbox1.bisect(obbox3, obbox4);
  EXPECT_TRUE(axom::primal::intersect(obbox3, obbox4));

  // check vertex-vertex intersection
  QOBBox obbox5(obbox1);
  obbox5.shift(-2. * shift);
  EXPECT_TRUE(axom::primal::intersect(obbox4, obbox5));
}

//------------------------------------------------------------------------------
TEST(primal_intersect, obb_obb_test_intersection3D)
{
  const int DIM = 3;
  typedef double CoordType;
  typedef primal::Point<CoordType, DIM> QPoint;
  typedef primal::Vector<CoordType, DIM> QVector;
  typedef primal::OrientedBoundingBox<CoordType, DIM> QOBBox;

  QPoint pt1;       // origin
  QVector u1[DIM];  // make standard axes
  QVector u2[DIM];  // first two axes rotated by 90 degrees
  QVector u3[DIM];  // all axes rotated
  for(int i = 0; i < DIM; i++)
  {
    u1[i] = QVector();
    u2[i] = QVector();
    u3[i] = QVector();
    u1[i][i] = 1.;
  }
  u2[0][0] = 1.;
  u2[0][1] = 1.;
  u2[1][0] = 1.;
  u2[1][1] = -1.;
  u2[2][2] = 1.;

  QVector e1(1.);
  QVector e2(1.41422);
  QPoint pt2(2.);

  QOBBox obbox1(pt1, u1, e1);
  QOBBox obbox2(pt2, u2, e2);

  // test edge-face
  EXPECT_TRUE(axom::primal::intersect(obbox1, obbox2));

  // pass obbox2 through obbox1 to test body-body
  QVector shift(-0.5);

  obbox2.shift(2. * shift);
  EXPECT_TRUE(axom::primal::intersect(obbox1, obbox2));

  obbox2.shift(100. * shift);
  EXPECT_FALSE(axom::primal::intersect(obbox1, obbox2));

  // test vertex-vertex
  QOBBox obbox3(obbox1);
  obbox3.shift(-2. * shift);
  EXPECT_TRUE(axom::primal::intersect(obbox1, obbox3));

  QVector shift2;
  shift2[0] = 1.;
  shift2[1] = 1.;
  QOBBox obbox4(obbox1);
  obbox4.shift(shift2);
  // edge-edge
  EXPECT_TRUE(axom::primal::intersect(obbox1, obbox4));

  QOBBox obbox5;
  QOBBox obbox6;
  obbox1.bisect(obbox5, obbox6);
  // face-face
  EXPECT_TRUE(axom::primal::intersect(obbox5, obbox6));

  // now for vertex-edge and vertex-face
  //
  // NOTE: here we are applying an arbitrary rotation to create obbox7. The
  // only essential point is that it is rotated so that it has a unique vertex
  // with smallest z coordinate, thus (once we shift it appropriately) we can
  // look at cases where one vertex (this smallest one) "kisses" parts of the
  // other bounding box.
  u3[0][0] = 0.7071;
  u3[0][2] = 0.7071;
  u3[1][0] = 0.5;
  u3[1][1] = 0.7071;
  u3[1][2] = -0.5;
  u3[2][0] = -0.5;
  u3[2][1] = 0.7071;
  u3[2][2] = 0.5;

  QVector e3(1.);
  QOBBox obbox7(pt1, u3, e3);
  QVector shift3;
  // shift so the lowest vertex is brushing the top face of obbox1 at (0,0,1)
  shift3[0] = -0.292893;
  shift3[1] = 0;
  shift3[2] = 1.70711 + 1.;
  obbox7.shift(shift3);
  // vertex-face
  EXPECT_TRUE(axom::primal::intersect(obbox1, obbox7));

  shift3[2] = 0.;
  shift3[1] = 1.;
  obbox7.shift(shift3);
  // vertex-edge
  EXPECT_TRUE(axom::primal::intersect(obbox1, obbox7));

  obbox7.shift(100. * shift3);
  EXPECT_FALSE(axom::primal::intersect(obbox1, obbox7));
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
#include "axom/slic/core/SimpleLogger.hpp"
using axom::slic::SimpleLogger;

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);

  SimpleLogger logger;  // create & initialize test logger,
  axom::slic::setLoggingMsgLevel(axom::slic::message::Warning);

  int result = RUN_ALL_TESTS();
  return result;
}

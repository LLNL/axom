// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file primal_surface_intersect.cpp
 * \brief This file tests surface (bilinear/patch/NURBS) intersection routines
 */

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/slic.hpp"

#include "axom/primal/geometry/BezierPatch.hpp"
#include "axom/primal/operators/intersect.hpp"

#include <cmath>
#include <iomanip>
#include <algorithm>
#include <sstream>

namespace primal = axom::primal;

/**
 * Helper function to compute the intersections of a Bezier patch and a ray 
 * and check that their intersection points match our expectations.
 * Patch parameters are stored in \a exp_u, \a exp_v and \a exp_t.
 * Intersections are computed within tolerance \a eps and our checks use \a test_eps.
 *
 * Param \a shouldPrintIntersections is used for debugging and for generating
 * the initial array of expected intersections.
 */
template <typename CoordType>
void checkIntersections(const primal::Ray<CoordType, 3>& ray,
                        const primal::BezierPatch<CoordType, 3>& patch,
                        const axom::Array<CoordType>& exp_u,
                        const axom::Array<CoordType>& exp_v,
                        const axom::Array<CoordType>& exp_t,
                        double eps,
                        double test_eps,
                        bool shouldPrintIntersections = false)
{
  constexpr int DIM = 3;
  using Array = axom::Array<CoordType>;

  // Check validity of input data.
  // They should have the same size
  EXPECT_EQ(exp_u.size(), exp_v.size());
  EXPECT_EQ(exp_u.size(), exp_t.size());

  const int num_exp_intersections = static_cast<int>(exp_u.size());
  const bool exp_intersect = (num_exp_intersections > 0);

  // Intersect the ray and the patch, intersection parameters will be
  // in arrays (u, v) and t, for the patch and ray, respectively
  Array u, v, t;
  bool ray_intersects = intersect(ray, patch, u, v, t);
  EXPECT_EQ(exp_intersect, ray_intersects);
  EXPECT_EQ(u.size(), v.size());
  EXPECT_EQ(u.size(), t.size());

  // check that we found the expected number of intersection points
  const int num_actual_intersections = static_cast<int>(u.size());
  EXPECT_EQ(num_exp_intersections, num_actual_intersections);

  // check that the evaluated intersection points are identical
  for(int i = 0; i < num_actual_intersections; ++i)
  {
    auto p1 = ray.at(t[i]);
    auto p2 = patch.evaluate(u[i], v[i]);

    EXPECT_NEAR(0., primal::squared_distance(p1, p2), test_eps);

    for(int d = 0; d < DIM; ++d)
    {
      EXPECT_NEAR(p1[d], p2[d], test_eps);
    }
  }

  if(shouldPrintIntersections)
  {
    std::stringstream sstr;

    sstr << "Intersections for ray and patch: "
         << "\n\t" << ray << "\n\t" << patch;

    sstr << "\ns (" << u.size() << "): ";
    for(auto i = 0u; i < u.size(); ++i)
    {
      sstr << std::setprecision(16) << "(" << u[i] << "," << v[i] << "),";
    }

    sstr << "\nt (" << t.size() << "): ";
    for(auto i = 0u; i < t.size(); ++i)
    {
      sstr << std::setprecision(16) << t[i] << ",";
    }

    SLIC_INFO(sstr.str());
  }

  for(int i = 0; i < num_actual_intersections; ++i)
  {
    EXPECT_NEAR(exp_u[i], u[i], test_eps);
    EXPECT_NEAR(exp_v[i], v[i], test_eps);
    EXPECT_NEAR(exp_t[i], t[i], test_eps);

    if(shouldPrintIntersections)
    {
      SLIC_INFO("\t" << i << ": {u:" << u[i] << ", v:" << v[i]
                     << std::setprecision(16) << ", t:" << t[i]
                     << ", u_actual:" << exp_u[i] << ", v_actual:" << exp_v[i]
                     << ", t_actual:" << exp_t[i] << "}");
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_surface_inter, bilinear_intersect)
{
  static const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType, DIM>;
  using RayType = primal::Ray<CoordType, DIM>;

  const double eps = 1E-16;
  const double eps_test = 1E-10;

  SLIC_INFO("primal: testing bilinear patch intersection");

  // Set control points
  BezierPatchType bilinear_patch(1, 1);
  bilinear_patch(0, 0) = PointType({-1.0, 1.0, 1.0});
  bilinear_patch(1, 0) = PointType({-1.0, -1.0, 2.0});
  bilinear_patch(1, 1) = PointType({1.0, -1.0, 1.0});
  bilinear_patch(0, 1) = PointType({1.0, 1.0, 2.0});

  // Ray with single intersection
  PointType ray_origin({0.0, 0.0, 1.75});
  VectorType ray_direction({1.0, 1.0, 0.0});
  RayType ray(ray_origin, ray_direction);

  checkIntersections(ray,
                     bilinear_patch,
                     {0.146446609407},
                     {0.853553390593},
                     {1.0},
                     eps,
                     eps_test);

  // Ray with no intersections
  ray_direction = VectorType({1.0, -1.0, 0.0});
  ray = RayType(ray_origin, ray_direction);
  checkIntersections(ray, bilinear_patch, {}, {}, {}, eps, eps_test);

  // Ray with no intersections, but in a way that is difficult for
  //  the standard GARP implementation
  ray_direction = VectorType({1.0, 0.0, 0.0});
  ray = RayType(ray_origin, ray_direction);
  checkIntersections(ray, bilinear_patch, {}, {}, {}, eps, eps_test);

  // Ray with two intersections
  ray_origin = PointType({-1.0, -1.0, 1.75});
  ray_direction = VectorType({1.0, 1.0, 0.0});
  ray = RayType(ray_origin, ray_direction);

  checkIntersections(ray,
                     bilinear_patch,
                     {0.853553390593, 0.146446609407},
                     {0.146446609407, 0.853553390593},
                     {0.414213562373, 2.41421356237},
                     eps,
                     eps_test);

  // Ray with infinitely many intersections
  ray_origin = PointType({-2.0, 0.0, 1.5});
  ray_direction = VectorType({1.0, 0.0, 0.0});
  ray = RayType(ray_origin, ray_direction);

  checkIntersections(ray, bilinear_patch, {0.5}, {0.0}, {1.0}, eps, eps_test);

  // Ray with no intersections on line with infinitely many intersections
  ray_origin = PointType({2.0, 0.0, 1.5});
  ray_direction = VectorType({1.0, 0.0, 0.0});
  ray = RayType(ray_origin, ray_direction);

  checkIntersections(ray, bilinear_patch, {}, {}, {}, eps, eps_test);

  // Ray with infinitely many intersections, but the origin is on the patch
  ray_origin = PointType({0.4, 0.0, 1.5});
  ray_direction = VectorType({1.0, 0.0, 0.0});
  ray = RayType(ray_origin, ray_direction);

  checkIntersections(ray, bilinear_patch, {0.5}, {0.7}, {0.0}, eps, eps_test);
}

//------------------------------------------------------------------------------
TEST(primal_surface_inter, difficult_garp_case)
{
  static const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType, DIM>;
  using RayType = primal::Ray<CoordType, DIM>;

  const double eps = 1E-16;
  const double eps_test = 1E-10;

  SLIC_INFO("primal: testing bilinear patch intersection");

  // Set control points
  BezierPatchType bilinear_patch(1, 1);
  bilinear_patch(0, 0) = PointType({-2.0, 1.0, 2.0});
  bilinear_patch(1, 0) = PointType({-1.0, -1.0, 1.0});
  bilinear_patch(1, 1) = PointType({1.0, -1.0, 1.0});
  bilinear_patch(0, 1) = PointType({2.0, 1.0, 1.0});

  VectorType ray_direction({-1.0, -2.0, 0.0});
  // The first step of the GARP algorithm is to solve a quadratic equation a + bt + ct^2 = 0,
  //   and this configuration of patch + ray direction is such that c = 0
  
  // Ray with single intersection
  PointType ray_origin({0.0, 1.0, 1.25});
  RayType ray(ray_origin, ray_direction);
  checkIntersections(ray, bilinear_patch, {2./3.}, {0.25}, {sqrt(20. / 9.)}, eps, eps_test);

  // Ray with no intersections
  ray_origin = PointType({0.0, 1.0, 1.75});
  ray = RayType(ray_origin, ray_direction);
  checkIntersections(ray, bilinear_patch, {}, {}, {}, eps, eps_test);
}

//------------------------------------------------------------------------------
TEST(primal_surface_inter, flat_bilinear_intersect)
{
  static const int DIM = 3;
  using CoordType = double;
  using PointType = primal::Point<CoordType, DIM>;
  using VectorType = primal::Vector<CoordType, DIM>;
  using BezierPatchType = primal::BezierPatch<CoordType, DIM>;
  using RayType = primal::Ray<CoordType, DIM>;

  const double eps = 1E-16;
  const double eps_test = 1E-10;

  SLIC_INFO("primal: testing bilinear patch intersection");

  // Set control points
  BezierPatchType bilinear_patch(1, 1);
  bilinear_patch(0, 0) = PointType({-2.0, 1.0, 1.0});
  bilinear_patch(1, 0) = PointType({-1.0, -1.0, 1.0});
  bilinear_patch(1, 1) = PointType({1.0, -1.0, 1.0});
  bilinear_patch(0, 1) = PointType({2.0, 1.0, 1.0});

  // Ray with single intersection
  PointType ray_origin({0.0, 0.0, 2.0});
  VectorType ray_direction({1, 0.5, -1.0});
  RayType ray(ray_origin, ray_direction);

  checkIntersections(ray, bilinear_patch, {0.25}, {11. / 14.}, {1.5}, eps, eps_test);

  // Ray with a single intersection that is coplanar with an isocurve
  ray_direction = VectorType({1.0, 0.0, -1.0});
  ray = RayType(ray_origin, ray_direction);
  checkIntersections(ray, bilinear_patch, {0.5}, {5. / 6.}, {sqrt(2)}, eps, eps_test);

  // Ray with no intersections
  ray_direction = VectorType({1.0, -1.0, -0.5});
  ray = RayType(ray_origin, ray_direction);
  checkIntersections(ray, bilinear_patch, {}, {}, {}, eps, eps_test);

  // Ray with no intersections that is parallel to the patch
  ray_direction = VectorType({1.0, 1.0, 0.0});
  ray = RayType(ray_origin, ray_direction);
  checkIntersections(ray, bilinear_patch, {}, {}, {}, eps, eps_test);

  // Ray with no intersections that is parallel to the patch along an isocurve
  ray_direction = VectorType({1.0, 0.0, 0.0});
  ray = RayType(ray_origin, ray_direction);
  checkIntersections(ray, bilinear_patch, {}, {}, {}, eps, eps_test);

  // Ray with no intersections on a line with infinitely many intersections
  ray_origin = PointType({-2.0, 0.5, 1});
  ray_direction = VectorType({-1.0, 0.0, 0.0});
  ray = RayType(ray_origin, ray_direction);

  checkIntersections(ray, bilinear_patch, {}, {}, {}, eps, eps_test);

  // Ray with infinitely many intersections
  ray_origin = PointType({-2.0, 0.0, 1.0});
  ray_direction = VectorType({1.0, 0.0, 0.0});
  ray = RayType(ray_origin, ray_direction);

  // checkIntersections(ray, bilinear_patch, {0.25}, {0.0}, {0.25}, eps, eps_test);
}

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;  // create & initialize test logger,

  result = RUN_ALL_TESTS();

  return result;
}

// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/slic.hpp"

#include "axom/core/Types.hpp"
#include "axom/core/execution/for_all.hpp"
#include "axom/core/memory_management.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/Plane.hpp"

#include "axom/primal/operators/clip.hpp"
#include "axom/primal/operators/intersection_volume.hpp"
#include "axom/primal/operators/split.hpp"

#include <limits>

namespace Primal3D
{
using PointType = axom::primal::Point<double, 3>;
using VectorType = axom::primal::Vector<double, 3>;
using BoundingBoxType = axom::primal::BoundingBox<double, 3>;
using HexahedronType = axom::primal::Hexahedron<double, 3>;
using TriangleType = axom::primal::Triangle<double, 3>;
using TetrahedronType = axom::primal::Tetrahedron<double, 3>;
using OctahedronType = axom::primal::Octahedron<double, 3>;
using PolyhedronType = axom::primal::Polyhedron<double, 3>;
using PolygonType = axom::primal::Polygon<double, 3>;
using PlaneType = axom::primal::Plane<double, 3>;
using PolyhedronType = axom::primal::Polyhedron<double, 3>;
}  // namespace Primal3D

TEST(primal_clip, simple_clip)
{
  using namespace Primal3D;
  BoundingBoxType bbox;
  bbox.addPoint(PointType::zero());
  bbox.addPoint(PointType::ones());

  PointType points[] = {
    PointType {2, 2, 2},
    PointType {2, 2, 4},
    PointType {2, 4, 2},

    PointType {-100, -100, 0.5},
    PointType {-100, 100, 0.5},
    PointType {100, 0, 0.5},

    PointType {0.25, 0.25, 0.5},
    PointType {0.75, 0.25, 0.5},
    PointType {0.66, 0.5, 0.5},
    PointType {1.5, 0.5, 0.5},
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

    EXPECT_EQ(PointType(.5), poly.vertexMean());

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
  constexpr double delta = 1e-5;

  // Test the "unit simplex", and a jittered version
  PointType points[] = {PointType {1, 0, 0},
                        PointType {0, 1, 0},
                        PointType {0, 0, 1},
                        PointType {1 + delta, delta, delta},
                        PointType {delta, 1 + delta, delta},
                        PointType {delta, delta, 1 + delta}};

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

  constexpr double VAL1 = 3.;
  constexpr double VAL2 = 2.;

  BoundingBoxType bbox;
  bbox.addPoint(PointType(-1.));
  bbox.addPoint(PointType(1.));

  PointType midpoint = PointType::zero();

  PointType points[] = {
    PointType {VAL1, VAL2, 0},
    PointType {-VAL1, VAL2, 0},
    PointType {VAL1, -VAL2, 0},
    PointType {-VAL1, -VAL2, 0},

    PointType {VAL1, 0, VAL2},
    PointType {-VAL1, 0, VAL2},
    PointType {VAL1, 0, -VAL2},
    PointType {-VAL1, 0, -VAL2},

    PointType {0, VAL2, VAL1},
    PointType {0, VAL2, -VAL1},
    PointType {0, -VAL2, VAL1},
    PointType {0, -VAL2, -VAL1},

    PointType {0, VAL1, VAL2},
    PointType {0, -VAL1, VAL2},
    PointType {0, VAL1, -VAL2},
    PointType {0, -VAL1, -VAL2},
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

  constexpr double EPS = 1e-8;

  // Triangle 248 from sphere mesh
  TriangleType tri(PointType {0.405431, 3.91921, 3.07821},
                   PointType {1.06511, 3.96325, 2.85626},
                   PointType {0.656002, 4.32465, 2.42221});

  // Block index {grid pt: (19,29,24); level: 5} from InOutOctree
  BoundingBoxType box12(PointType {0.937594, 4.06291, 2.50025},
                        PointType {1.25012, 4.37544, 2.81278});

  PolygonType poly = axom::primal::clip(tri, box12);
  EXPECT_EQ(3, poly.numVertices());

  SLIC_INFO("Intersection of triangle "
            << tri << " \n\t and bounding box " << box12 << " \n\t is polygon"
            << poly << " with centroid " << poly.vertexMean());

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
    PointType centroid = poly.vertexMean();
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

template <typename ExecPolicy>
void unit_check_poly_clip()
{
  using namespace Primal3D;

  constexpr double EPS = 1.e-24;

  // Create a square in 3D space as our "polyhedron", and a plane splitting it
  // between the middle. The polyhedron after a call to poly_clip_vertices
  // should have two new vertices where the plane crosses the 0-1 and 2-3 edges:
  //      |
  // 3 ---|--- 2        3 -- 2' -- 2
  // |    |    |        |          |
  // |    p->  |   -->  |          |
  // |    |    |        |          |
  // 0 ---|--- 1        0 -- 1' -- 1
  //      |
  // In addition, vertices 0 and 3 should be marked as clipped.

  PolyhedronType square;
  square.addVertex({0.0, 0.0, 0.0});
  square.addVertex({1.0, 0.0, 0.0});
  square.addVertex({1.0, 1.0, 0.0});
  square.addVertex({0.0, 1.0, 0.0});

  square.addNeighbors(0, {1, 3});
  square.addNeighbors(1, {0, 2});
  square.addNeighbors(2, {1, 3});
  square.addNeighbors(3, {0, 2});

  PlaneType plane(VectorType {1, 0, 0}, PointType {0.5, 0.0, 0.0});

  const int current_allocator = axom::getDefaultAllocatorID();
  axom::setDefaultAllocator(axom::execution_space<ExecPolicy>::allocatorID());

  PolyhedronType* out_square = axom::allocate<PolyhedronType>(1);
  out_square[0] = square;

  unsigned int* out_clipped = axom::allocate<unsigned int>(1);
  out_clipped[0] = 0;

  axom::for_all<ExecPolicy>(
    0,
    1,
    AXOM_LAMBDA(int /* idx */) {
      axom::primal::detail::poly_clip_vertices(out_square[0],
                                               plane,
                                               EPS,
                                               out_clipped[0]);
    });

  const PolyhedronType& clippedSquare = out_square[0];

  EXPECT_EQ(clippedSquare.numVertices(), 6);
  // Check that existing vertices were not modified
  EXPECT_EQ(clippedSquare[0], square[0]);
  EXPECT_EQ(clippedSquare[1], square[1]);
  EXPECT_EQ(clippedSquare[2], square[2]);
  EXPECT_EQ(clippedSquare[3], square[3]);

  int lo_idx = -1, hi_idx = -1;
  for(int i = 4; i < 6; i++)
  {
    if(clippedSquare[i] == PointType {0.5, 0.0, 0.0})
    {
      lo_idx = i;
    }
    else if(clippedSquare[i] == PointType {0.5, 1.0, 0.0})
    {
      hi_idx = i;
    }
  }

  // Check that plane-crossing vertices were generated
  EXPECT_NE(lo_idx, -1);
  EXPECT_NE(hi_idx, -1);

  // Check that vertices outside the plane were marked as clipped
  EXPECT_EQ(out_clipped[0], (1 | (1 << 3)));

  // Generate sets of expected neighbors
  std::vector<std::set<int>> expectedNbrs(6);
  expectedNbrs[0] = {3, lo_idx};
  expectedNbrs[1] = {2, lo_idx};
  expectedNbrs[2] = {1, hi_idx};
  expectedNbrs[3] = {0, hi_idx};
  expectedNbrs[lo_idx] = {0, 1};
  expectedNbrs[hi_idx] = {2, 3};

  // Check connectivity
  for(int iv = 0; iv < 6; iv++)
  {
    // Each vertex should only have two neighbors
    EXPECT_EQ(clippedSquare.getNumNeighbors(iv), 2);
    std::set<int> resultNbrs {clippedSquare.getNeighbors(iv)[0],
                              clippedSquare.getNeighbors(iv)[1]};
    // Check that neighbors match expected connectivity
    EXPECT_EQ(expectedNbrs[iv], resultNbrs);
  }

  axom::deallocate(out_square);
  axom::deallocate(out_clipped);
  axom::setDefaultAllocator(current_allocator);
}

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
template <typename ExecPolicy>
void check_hex_tet_clip(double EPS)
{
  using namespace Primal3D;

  constexpr bool CHECK_SIGN = true;

  // Save current/default allocator
  const int current_allocator = axom::getDefaultAllocatorID();

  // Set new default to device if available
  axom::setDefaultAllocator(axom::execution_space<ExecPolicy>::allocatorID());

  // Allocate memory for shapes
  TetrahedronType* tet = axom::allocate<TetrahedronType>(1);
  HexahedronType* hex = axom::allocate<HexahedronType>(1);
  PolyhedronType* res = axom::allocate<PolyhedronType>(1);

  tet[0] = TetrahedronType(PointType {1, 0, 0},
                           PointType {1, 1, 0},
                           PointType {0, 1, 0},
                           PointType {1, 0, 1});

  hex[0] = HexahedronType(PointType {0, 0, 0},
                          PointType {1, 0, 0},
                          PointType {1, 1, 0},
                          PointType {0, 1, 0},
                          PointType {0, 0, 1},
                          PointType {1, 0, 1},
                          PointType {1, 1, 1},
                          PointType {0, 1, 1});
  axom::for_all<ExecPolicy>(
    1,
    AXOM_LAMBDA(int i) { res[i] = axom::primal::clip(hex[i], tet[i]); });

  EXPECT_NEAR(0.1666, res[0].volume(), EPS);
  EXPECT_NEAR(0.1666,
              axom::primal::intersection_volume<double>(hex[0], tet[0]),
              EPS);

  // Test checkSign optional parameter using shapes with negative volumes
  axom::utilities::swap<PointType>(tet[0][1], tet[0][2]);
  axom::utilities::swap<PointType>(hex[0][1], hex[0][3]);
  axom::utilities::swap<PointType>(hex[0][5], hex[0][7]);

  EXPECT_LT(tet[0].signedVolume(), 0.0);
  EXPECT_LT(hex[0].signedVolume(), 0.0);

  axom::for_all<ExecPolicy>(
    1,
    AXOM_LAMBDA(int i) {
      res[i] = axom::primal::clip(hex[i], tet[i], EPS, CHECK_SIGN);
    });

  EXPECT_NEAR(0.1666, res[0].volume(), EPS);
  EXPECT_NEAR(
    0.1666,
    axom::primal::intersection_volume<double>(hex[0], tet[0], EPS, CHECK_SIGN),
    EPS);

  axom::deallocate(tet);
  axom::deallocate(hex);
  axom::deallocate(res);

  axom::setDefaultAllocator(current_allocator);
}

template <typename ExecPolicy>
void check_oct_tet_clip(double EPS)
{
  using namespace Primal3D;

  constexpr bool CHECK_SIGN = true;

  // Save current/default allocator
  const int current_allocator = axom::getDefaultAllocatorID();

  // Set new default to device if available
  axom::setDefaultAllocator(axom::execution_space<ExecPolicy>::allocatorID());

  // Allocate memory for shapes
  TetrahedronType* tet = axom::allocate<TetrahedronType>(1);
  OctahedronType* oct = axom::allocate<OctahedronType>(1);
  PolyhedronType* res = axom::allocate<PolyhedronType>(1);

  tet[0] = TetrahedronType(PointType {1, 0, 0},
                           PointType {1, 1, 0},
                           PointType {0, 1, 0},
                           PointType {1, 0, 1});

  oct[0] = OctahedronType(PointType {1, 0, 0},
                          PointType {1, 1, 0},
                          PointType {0, 1, 0},
                          PointType {0, 1, 1},
                          PointType {0, 0, 1},
                          PointType {1, 0, 1});

  axom::for_all<ExecPolicy>(
    1,
    AXOM_LAMBDA(int i) { res[i] = axom::primal::clip(oct[i], tet[i]); });

  EXPECT_NEAR(0.1666, res[0].volume(), EPS);
  EXPECT_NEAR(0.1666,
              axom::primal::intersection_volume<double>(oct[0], tet[0]),
              EPS);

  // Test checkSign optional parameter using shapes with negative volumes
  axom::utilities::swap<PointType>(tet[0][1], tet[0][2]);
  axom::utilities::swap<PointType>(oct[0][1], oct[0][2]);
  axom::utilities::swap<PointType>(oct[0][4], oct[0][5]);

  EXPECT_LT(tet[0].signedVolume(), 0.0);

  axom::for_all<ExecPolicy>(
    1,
    AXOM_LAMBDA(int i) {
      res[i] = axom::primal::clip(oct[i], tet[i], EPS, CHECK_SIGN);
    });

  EXPECT_NEAR(0.1666, res[0].volume(), EPS);
  EXPECT_NEAR(
    0.1666,
    axom::primal::intersection_volume<double>(oct[0], tet[0], EPS, CHECK_SIGN),
    EPS);

  axom::deallocate(tet);
  axom::deallocate(oct);
  axom::deallocate(res);

  axom::setDefaultAllocator(current_allocator);
}

template <typename ExecPolicy>
void check_tet_tet_clip(double EPS)
{
  using namespace Primal3D;

  constexpr bool CHECK_SIGN = true;

  // Save current/default allocator
  const int current_allocator = axom::getDefaultAllocatorID();

  // Set new default to device if available
  axom::setDefaultAllocator(axom::execution_space<ExecPolicy>::allocatorID());

  // Allocate memory for shapes
  TetrahedronType* tet1 = axom::allocate<TetrahedronType>(1);
  TetrahedronType* tet2 = axom::allocate<TetrahedronType>(1);
  PolyhedronType* res = axom::allocate<PolyhedronType>(1);

  tet1[0] = TetrahedronType(PointType {1, 0, 0},
                            PointType {1, 1, 0},
                            PointType {0, 1, 0},
                            PointType {1, 1, 1});

  tet2[0] = TetrahedronType(PointType {0, 0, 0},
                            PointType {1, 0, 0},
                            PointType {1, 1, 0},
                            PointType {1, 1, 1});

  axom::for_all<ExecPolicy>(
    1,
    AXOM_LAMBDA(int i) { res[i] = axom::primal::clip(tet1[i], tet2[i]); });

  EXPECT_NEAR(0.0833, res[0].volume(), EPS);
  EXPECT_NEAR(0.0833,
              axom::primal::intersection_volume<double>(tet1[0], tet2[0]),
              EPS);

  // Test checkSign optional parameter using shapes with negative volumes
  axom::utilities::swap<PointType>(tet1[0][1], tet1[0][2]);
  axom::utilities::swap<PointType>(tet2[0][1], tet2[0][2]);

  EXPECT_LT(tet1[0].signedVolume(), 0.0);
  EXPECT_LT(tet2[0].signedVolume(), 0.0);

  axom::for_all<ExecPolicy>(
    1,
    AXOM_LAMBDA(int i) {
      res[i] = axom::primal::clip(tet1[i], tet2[i], EPS, CHECK_SIGN);
    });

  EXPECT_NEAR(0.0833, res[0].volume(), EPS);
  EXPECT_NEAR(
    0.0833,
    axom::primal::intersection_volume<double>(tet1[0], tet2[0], EPS, CHECK_SIGN),
    EPS);

  axom::deallocate(tet1);
  axom::deallocate(tet2);
  axom::deallocate(res);

  axom::setDefaultAllocator(current_allocator);
}

TEST(primal_clip, unit_poly_clip_vertices_sequential)
{
  unit_check_poly_clip<axom::SEQ_EXEC>();
}

TEST(primal_clip, clip_hex_tet_sequential)
{
  constexpr double EPS = 1e-4;
  check_hex_tet_clip<axom::SEQ_EXEC>(EPS);
}

TEST(primal_clip, clip_oct_tet_sequential)
{
  constexpr double EPS = 1e-4;
  check_oct_tet_clip<axom::SEQ_EXEC>(EPS);
}

TEST(primal_clip, clip_tet_tet_sequential)
{
  constexpr double EPS = 1e-4;
  check_tet_tet_clip<axom::SEQ_EXEC>(EPS);
}

  #ifdef AXOM_USE_OPENMP
TEST(primal_clip, unit_poly_clip_vertices_omp)
{
  unit_check_poly_clip<axom::OMP_EXEC>();
}

TEST(primal_clip, clip_hex_tet_omp)
{
  constexpr double EPS = 1e-4;
  check_hex_tet_clip<axom::OMP_EXEC>(EPS);
}

TEST(primal_clip, clip_oct_tet_omp)
{
  constexpr double EPS = 1e-4;
  check_oct_tet_clip<axom::OMP_EXEC>(EPS);
}

TEST(primal_clip, clip_tet_tet_omp)
{
  constexpr double EPS = 1e-4;
  check_tet_tet_clip<axom::OMP_EXEC>(EPS);
}
  #endif /* AXOM_USE_OPENMP */

  #if defined(AXOM_USE_CUDA)
AXOM_CUDA_TEST(primal_clip, unit_poly_clip_vertices_cuda)
{
  unit_check_poly_clip<axom::CUDA_EXEC<256>>();
}

AXOM_CUDA_TEST(primal_clip, clip_hex_tet_cuda)
{
  constexpr double EPS = 1e-4;
  check_hex_tet_clip<axom::CUDA_EXEC<256>>(EPS);
}

AXOM_CUDA_TEST(primal_clip, clip_oct_tet_cuda)
{
  constexpr double EPS = 1e-4;
  check_oct_tet_clip<axom::CUDA_EXEC<256>>(EPS);
}

TEST(primal_clip, clip_tet_tet_cuda)
{
  constexpr double EPS = 1e-4;
  check_tet_tet_clip<axom::CUDA_EXEC<256>>(EPS);
}
  #endif /* AXOM_USE_CUDA */

  #if defined(AXOM_USE_HIP)
TEST(primal_clip, unit_poly_clip_vertices_hip)
{
  unit_check_poly_clip<axom::HIP_EXEC<256>>();
}

TEST(primal_clip, clip_hex_tet_hip)
{
  constexpr double EPS = 1e-4;
  check_hex_tet_clip<axom::HIP_EXEC<256>>(EPS);
}

TEST(primal_clip, clip_oct_tet_hip)
{
  constexpr double EPS = 1e-4;
  check_oct_tet_clip<axom::HIP_EXEC<256>>(EPS);
}

TEST(primal_clip, clip_tet_tet_hip)
{
  constexpr double EPS = 1e-4;
  check_tet_tet_clip<axom::HIP_EXEC<256>>(EPS);
}
  #endif /* AXOM_USE_HIP */

#endif /* AXOM_USE_RAJA && AXOM_USE_UMPIRE */

// Tetrahedron does not clip hexahedron.
TEST(primal_clip, hex_tet_clip_nonintersect)
{
  using namespace Primal3D;

  TetrahedronType tet(PointType {-1, -1, -1},
                      PointType {-1, 0, 0},
                      PointType {-1, -1, 0},
                      PointType {0, 0, 0});
  HexahedronType hex(PointType {1, 0, 0},
                     PointType {2, 0, 0},
                     PointType {2, 1, 0},
                     PointType {1, 1, 0},
                     PointType {1, 0, 1},
                     PointType {2, 0, 1},
                     PointType {2, 1, 1},
                     PointType {1, 1, 1});

  PolyhedronType poly = axom::primal::clip(hex, tet);
  EXPECT_EQ(0.0, poly.volume());
  EXPECT_EQ(0.0, axom::primal::intersection_volume<double>(hex, tet));
  EXPECT_EQ(0.0, axom::primal::intersection_volume<double>(tet, hex));
  EXPECT_EQ(axom::primal::clip(tet, hex).volume(), poly.volume());
}

// Tetrahedron is encapsulated by the hexahedron
TEST(primal_clip, hex_tet_clip_encapsulate)
{
  using namespace Primal3D;
  constexpr double EPS = 1e-4;

  TetrahedronType tet(PointType {1, 0, 0},
                      PointType {1, 1, 0},
                      PointType {0, 1, 0},
                      PointType {1, 0, 1});
  HexahedronType hex(PointType {0, 0, 0},
                     PointType {1, 0, 0},
                     PointType {1, 1, 0},
                     PointType {0, 1, 0},
                     PointType {0, 0, 1},
                     PointType {1, 0, 1},
                     PointType {1, 1, 1},
                     PointType {0, 1, 1});

  PolyhedronType poly = axom::primal::clip(hex, tet);

  // Expected result should be 0.666 / 4 = 0.1666, volume of tet.
  EXPECT_NEAR(0.1666, poly.volume(), EPS);
  EXPECT_NEAR(0.1666, axom::primal::intersection_volume<double>(hex, tet), EPS);
  EXPECT_NEAR(0.1666, axom::primal::intersection_volume<double>(tet, hex), EPS);
  EXPECT_NEAR(axom::primal::clip(tet, hex).volume(), poly.volume(), EPS);
}

// Hexahedron is encapsulated inside the tetrahedron
TEST(primal_clip, hex_tet_clip_encapsulate_inv)
{
  using namespace Primal3D;
  constexpr double EPS = 1e-4;

  TetrahedronType tet(PointType {0, 0, 0},
                      PointType {0, 3, 0},
                      PointType {0, 0, 3},
                      PointType {3, 0, 0});
  HexahedronType hex(PointType {0, 0, 0},
                     PointType {1, 0, 0},
                     PointType {1, 1, 0},
                     PointType {0, 1, 0},
                     PointType {0, 0, 1},
                     PointType {1, 0, 1},
                     PointType {1, 1, 1},
                     PointType {0, 1, 1});

  PolyhedronType poly = axom::primal::clip(hex, tet);

  // Expected result should be 1.0, volume of hex.
  EXPECT_NEAR(1.0, poly.volume(), EPS);
  EXPECT_NEAR(1.0, axom::primal::intersection_volume<double>(hex, tet), EPS);
  EXPECT_NEAR(1.0, axom::primal::intersection_volume<double>(tet, hex), EPS);
  EXPECT_NEAR(axom::primal::clip(tet, hex).volume(), poly.volume(), EPS);
}

// Half of the hexahedron is clipped by the tetrahedron
TEST(primal_clip, hex_tet_clip_half)
{
  using namespace Primal3D;
  constexpr double EPS = 1e-4;

  TetrahedronType tet(PointType {0, 0, 0},
                      PointType {0, 3, 0},
                      PointType {0, 0, 3},
                      PointType {3, 0, 0});
  HexahedronType hex(PointType {-1, 0, 0},
                     PointType {1, 0, 0},
                     PointType {1, 1, 0},
                     PointType {-1, 1, 0},
                     PointType {-1, 0, 1},
                     PointType {1, 0, 1},
                     PointType {1, 1, 1},
                     PointType {-1, 1, 1});

  PolyhedronType poly = axom::primal::clip(hex, tet);

  // Expected result should be 1.0, half the volume of hex.
  EXPECT_NEAR(1.0, poly.volume(), EPS);
  EXPECT_NEAR(1.0, axom::primal::intersection_volume<double>(hex, tet), EPS);
  EXPECT_NEAR(1.0, axom::primal::intersection_volume<double>(tet, hex), EPS);
  EXPECT_NEAR(axom::primal::clip(tet, hex).volume(), poly.volume(), EPS);
}

// Hexahedron is adjacent to tetrahedron
TEST(primal_clip, hex_tet_clip_adjacent)
{
  using namespace Primal3D;

  TetrahedronType tet(PointType {0, -1, 0},
                      PointType {0, 0, 1},
                      PointType {1, 0, 1},
                      PointType {1, 0, 0});
  HexahedronType hex(PointType {0, 0, 0},
                     PointType {1, 0, 0},
                     PointType {1, 1, 0},
                     PointType {0, 1, 0},
                     PointType {0, 0, 1},
                     PointType {1, 0, 1},
                     PointType {1, 1, 1},
                     PointType {0, 1, 1});

  PolyhedronType poly = axom::primal::clip(hex, tet);

  EXPECT_EQ(0.0, poly.volume());
  EXPECT_EQ(0.0, axom::primal::intersection_volume<double>(hex, tet));
  EXPECT_EQ(0.0, axom::primal::intersection_volume<double>(tet, hex));
  EXPECT_EQ(axom::primal::clip(tet, hex).volume(), poly.volume());
}

// Tetrahedron clips hexahedron at a single vertex
TEST(primal_clip, hex_tet_clip_point)
{
  using namespace Primal3D;

  TetrahedronType tet(PointType {0, -1, 0},
                      PointType {0, 0, 1},
                      PointType {1, 0, 1},
                      PointType {1, 0, 0});
  HexahedronType hex(PointType {-1, -2, -1},
                     PointType {1, -2, -1},
                     PointType {1, -1, -1},
                     PointType {-1, -1, -1},
                     PointType {-1, -2, 1},
                     PointType {1, -2, 1},
                     PointType {1, -1, 1},
                     PointType {-1, -1, 1});

  PolyhedronType poly = axom::primal::clip(hex, tet);

  EXPECT_EQ(0.0, poly.volume());
  EXPECT_EQ(0.0, axom::primal::intersection_volume<double>(hex, tet));
  EXPECT_EQ(0.0, axom::primal::intersection_volume<double>(tet, hex));
  EXPECT_EQ(axom::primal::clip(tet, hex).volume(), poly.volume());
}

// Tetrahedron does not clip octahedron.
TEST(primal_clip, oct_tet_clip_nonintersect)
{
  using namespace Primal3D;

  TetrahedronType tet(PointType {-1, -1, -1},
                      PointType {-1, 0, 0},
                      PointType {-1, -1, 0},
                      PointType {0, 0, 0});
  OctahedronType oct(PointType {1, 0, 0},
                     PointType {1, 1, 0},
                     PointType {0, 1, 0},
                     PointType {0, 1, 1},
                     PointType {0, 0, 1},
                     PointType {1, 0, 1});

  PolyhedronType poly = axom::primal::clip(oct, tet);
  EXPECT_EQ(0.0, poly.volume());
  EXPECT_EQ(0.0, axom::primal::intersection_volume<double>(oct, tet));
  EXPECT_EQ(0.0, axom::primal::intersection_volume<double>(tet, oct));
  EXPECT_EQ(axom::primal::clip(tet, oct).volume(), poly.volume());
}

// Tetrahedron is encapsulated by the octahedron
TEST(primal_clip, oct_tet_clip_encapsulate)
{
  using namespace Primal3D;
  constexpr double EPS = 1e-4;

  TetrahedronType tet(PointType {1, 0, 0},
                      PointType {1, 1, 0},
                      PointType {0, 1, 0},
                      PointType {1, 0, 1});
  OctahedronType oct(PointType {1, 0, 0},
                     PointType {1, 1, 0},
                     PointType {0, 1, 0},
                     PointType {0, 1, 1},
                     PointType {0, 0, 1},
                     PointType {1, 0, 1});

  PolyhedronType poly = axom::primal::clip(oct, tet);

  // Expected result should be 0.666 / 4 = 0.1666, volume of tet.
  EXPECT_NEAR(0.1666, poly.volume(), EPS);
  EXPECT_NEAR(0.1666, axom::primal::intersection_volume<double>(oct, tet), EPS);
  EXPECT_NEAR(0.1666, axom::primal::intersection_volume<double>(tet, oct), EPS);
  EXPECT_NEAR(axom::primal::clip(tet, oct).volume(), poly.volume(), EPS);
}

// Octahedron is encapsulated inside the tetrahedron
TEST(primal_clip, oct_tet_clip_encapsulate_inv)
{
  using namespace Primal3D;
  constexpr double EPS = 1e-4;

  TetrahedronType tet(PointType {0, 0, 0},
                      PointType {0, 2, 0},
                      PointType {0, 0, 2},
                      PointType {2, 0, 0});
  OctahedronType oct(PointType {1, 0, 0},
                     PointType {1, 1, 0},
                     PointType {0, 1, 0},
                     PointType {0, 1, 1},
                     PointType {0, 0, 1},
                     PointType {1, 0, 1});

  PolyhedronType poly = axom::primal::clip(oct, tet);

  // Expected result should be 0.6666, volume of oct.
  EXPECT_NEAR(0.6666, poly.volume(), EPS);
  EXPECT_NEAR(0.6666, axom::primal::intersection_volume<double>(oct, tet), EPS);
  EXPECT_NEAR(0.6666, axom::primal::intersection_volume<double>(tet, oct), EPS);
  EXPECT_NEAR(axom::primal::clip(tet, oct).volume(), poly.volume(), EPS);
}

// Half of the octahedron is clipped by the tetrahedron
TEST(primal_clip, oct_tet_clip_half)
{
  using namespace Primal3D;
  constexpr double EPS = 1e-4;

  TetrahedronType tet(PointType {0.5, 0.5, 2},
                      PointType {2, -1, 0},
                      PointType {-1, -1, 0},
                      PointType {-1, 2, 0});
  OctahedronType oct(PointType {1, 0, 0},
                     PointType {1, 1, 0},
                     PointType {0, 1, 0},
                     PointType {0, 1, 1},
                     PointType {0, 0, 1},
                     PointType {1, 0, 1});

  PolyhedronType poly = axom::primal::clip(oct, tet);

  // Expected result should be 0.3333, half the volume of oct.
  EXPECT_NEAR(0.3333, poly.volume(), EPS);
  EXPECT_NEAR(0.3333, axom::primal::intersection_volume<double>(oct, tet), EPS);
  EXPECT_NEAR(0.3333, axom::primal::intersection_volume<double>(tet, oct), EPS);
  EXPECT_NEAR(axom::primal::clip(tet, oct).volume(), poly.volume(), EPS);
}

// Octahedron is adjacent to tetrahedron
TEST(primal_clip, oct_tet_clip_adjacent)
{
  using namespace Primal3D;

  TetrahedronType tet(PointType {0, -1, 0},
                      PointType {0, 0, 1},
                      PointType {1, 0, 1},
                      PointType {1, 0, 0});
  OctahedronType oct(PointType {1, 0, 0},
                     PointType {1, 1, 0},
                     PointType {0, 1, 0},
                     PointType {0, 1, 1},
                     PointType {0, 0, 1},
                     PointType {1, 0, 1});

  PolyhedronType poly = axom::primal::clip(oct, tet);

  EXPECT_EQ(0.0, poly.volume());
  EXPECT_EQ(0.0, axom::primal::intersection_volume<double>(oct, tet));
  EXPECT_EQ(0.0, axom::primal::intersection_volume<double>(tet, oct));
  EXPECT_EQ(axom::primal::clip(tet, oct).volume(), poly.volume());
}

// Tetrahedron clips octahedron at a single vertex
TEST(primal_clip, oct_tet_clip_point)
{
  using namespace Primal3D;

  TetrahedronType tet(PointType {-1, -1, 0},
                      PointType {-0.5, 0.5, 0},
                      PointType {0, 0, 2},
                      PointType {0.5, -0.5, 0});
  OctahedronType oct(PointType {1, 0, 0},
                     PointType {1, 1, 0},
                     PointType {0, 1, 0},
                     PointType {0, 1, 1},
                     PointType {0, 0, 1},
                     PointType {1, 0, 1});

  PolyhedronType poly = axom::primal::clip(oct, tet);

  EXPECT_EQ(0.0, poly.volume());
  EXPECT_EQ(0.0, axom::primal::intersection_volume<double>(oct, tet));
  EXPECT_EQ(0.0, axom::primal::intersection_volume<double>(tet, oct));
  EXPECT_EQ(axom::primal::clip(tet, oct).volume(), poly.volume());
}

TEST(primal_clip, oct_tet_clip_special_case_1)
{
  using namespace Primal3D;
  constexpr double EPS = 1e-4;
  constexpr bool CHECK_SIGN = true;

  TetrahedronType tet(PointType {0.5, 0.5, 0.5},
                      PointType {1, 1, 0},
                      PointType {1, 0, 0},
                      PointType {0.5, 0.5, 0});

  OctahedronType oct(PointType {0.5, 0.853553, 0.146447},
                     PointType {0.853553, 0.853553, 0.5},
                     PointType {0.853553, 0.5, 0.146447},
                     PointType {1, 0.5, 0.5},
                     PointType {0.5, 0.5, 0},
                     PointType {0.5, 1, 0.5});

  // NOTE: Order of vertices 1,2 and 4,5 are flipped due
  // to winding being opposite of what volume() expects
  PolyhedronType octPoly;
  octPoly.addVertex(PointType {0.5, 0.853553, 0.146447});
  octPoly.addVertex(PointType {0.853553, 0.5, 0.146447});
  octPoly.addVertex(PointType {0.853553, 0.853553, 0.5});
  octPoly.addVertex(PointType {1, 0.5, 0.5});
  octPoly.addVertex(PointType {0.5, 1, 0.5});
  octPoly.addVertex(PointType {0.5, 0.5, 0});

  octPoly.addNeighbors(0, {1, 5, 4, 2});
  octPoly.addNeighbors(1, {0, 2, 3, 5});
  octPoly.addNeighbors(2, {0, 4, 3, 1});
  octPoly.addNeighbors(3, {1, 2, 4, 5});
  octPoly.addNeighbors(4, {0, 5, 3, 2});
  octPoly.addNeighbors(5, {0, 1, 3, 4});

  EXPECT_NEAR(0.0251, octPoly.volume(), EPS);

  PolyhedronType poly = axom::primal::clip(oct, tet, EPS, CHECK_SIGN);

  EXPECT_NEAR(0.0041, poly.volume(), EPS);
  EXPECT_NEAR(
    0.0041,
    axom::primal::intersection_volume<double>(oct, tet, EPS, CHECK_SIGN),
    EPS);
  EXPECT_NEAR(
    0.0041,
    axom::primal::intersection_volume<double>(tet, oct, EPS, CHECK_SIGN),
    EPS);
  EXPECT_NEAR(axom::primal::clip(tet, oct, EPS, CHECK_SIGN).volume(),
              poly.volume(),
              EPS);
}

TEST(primal_clip, oct_tet_clip_special_case_2)
{
  using namespace Primal3D;
  constexpr double EPS = 1e-4;
  constexpr bool CHECK_SIGN = true;

  TetrahedronType tet(PointType {0.5, 0.5, 0.5},
                      PointType {0, 1, 0},
                      PointType {1, 1, 0},
                      PointType {0.5, 0.5, 0});

  OctahedronType oct(PointType {0.5, 0.853553, 0.146447},
                     PointType {0.853553, 0.853553, 0.5},
                     PointType {0.853553, 0.5, 0.146447},
                     PointType {1, 0.5, 0.5},
                     PointType {0.5, 0.5, 0},
                     PointType {0.5, 1, 0.5});

  // NOTE: Order of vertices 1,2 and 4,5 are flipped due
  // to winding being opposite of what volume() expects
  PolyhedronType octPoly;
  octPoly.addVertex(PointType {0.5, 0.853553, 0.146447});
  octPoly.addVertex(PointType {0.853553, 0.5, 0.146447});
  octPoly.addVertex(PointType {0.853553, 0.853553, 0.5});
  octPoly.addVertex(PointType {1, 0.5, 0.5});
  octPoly.addVertex(PointType {0.5, 1, 0.5});
  octPoly.addVertex(PointType {0.5, 0.5, 0});

  octPoly.addNeighbors(0, {1, 5, 4, 2});
  octPoly.addNeighbors(1, {0, 2, 3, 5});
  octPoly.addNeighbors(2, {0, 4, 3, 1});
  octPoly.addNeighbors(3, {1, 2, 4, 5});
  octPoly.addNeighbors(4, {0, 5, 3, 2});
  octPoly.addNeighbors(5, {0, 1, 3, 4});

  EXPECT_NEAR(0.0251, octPoly.volume(), EPS);

  PolyhedronType poly = axom::primal::clip(oct, tet, EPS, CHECK_SIGN);

  EXPECT_NEAR(0.0041, poly.volume(), EPS);
  EXPECT_NEAR(
    0.0041,
    axom::primal::intersection_volume<double>(oct, tet, EPS, CHECK_SIGN),
    EPS);
  EXPECT_NEAR(
    0.0041,
    axom::primal::intersection_volume<double>(tet, oct, EPS, CHECK_SIGN),
    EPS);
  EXPECT_NEAR(axom::primal::clip(tet, oct, EPS, CHECK_SIGN).volume(),
              poly.volume(),
              EPS);
}

// Tetrahedron does not clip tetrahedron.
TEST(primal_clip, tet_tet_clip_nonintersect)
{
  using namespace Primal3D;

  {
    TetrahedronType tet1(PointType {-1, -1, -1},
                         PointType {-1, 0, 0},
                         PointType {-1, -1, 0},
                         PointType {0, 0, 0});
    EXPECT_TRUE(tet1.signedVolume() > 0.);

    TetrahedronType tet2(PointType {1, 0, 0},
                         PointType {1, 1, 0},
                         PointType {0, 1, 0},
                         PointType {1, 0, 1});
    EXPECT_TRUE(tet2.signedVolume() > 0.);

    PolyhedronType poly = axom::primal::clip(tet1, tet2);
    EXPECT_EQ(0.0, poly.volume());
    EXPECT_EQ(0.0, axom::primal::intersection_volume<double>(tet1, tet2));
  }

  // User-provided test case; tetrahedra do not overlap
  {
    TetrahedronType tet1(PointType {0, 1, 0},
                         PointType {-1, 0, 0},
                         PointType {0, -1, 0},
                         PointType {-0.25, 0, 0.25});
    EXPECT_TRUE(tet1.signedVolume() > 0.);

    TetrahedronType tet2(PointType {0.74899999999999922, 0, 1},
                         PointType {1.7489999999999992, 0, 0},
                         PointType {0.74899999999999922, -1, 0},
                         PointType {0.99899999999999922, 0, 0.25});
    EXPECT_TRUE(tet2.signedVolume() > 0.);

    PolyhedronType poly = axom::primal::clip(tet1, tet2);
    EXPECT_EQ(0.0, poly.volume());
    EXPECT_EQ(0.0, axom::primal::intersection_volume<double>(tet1, tet2));
  }
}

// Tetrahedron is adjacent to tetrahedron
TEST(primal_clip, tet_tet_clip_adjacent)
{
  using namespace Primal3D;

  TetrahedronType tet1(PointType {1, 0, 0},
                       PointType {1, 1, 0},
                       PointType {0, 1, 0},
                       PointType {1, 0, 1});

  TetrahedronType tet2(PointType {1, 0, 1},
                       PointType {0, 1, 0},
                       PointType {1, 0, 0},
                       PointType {0, 0, 0});

  PolyhedronType poly = axom::primal::clip(tet1, tet2);
  EXPECT_EQ(0.0, poly.volume());
  EXPECT_EQ(0.0, axom::primal::intersection_volume<double>(tet1, tet2));
}

// Tetrahedron clips tetrahedron at a single vertex
TEST(primal_clip, tet_tet_clip_point)
{
  using namespace Primal3D;

  TetrahedronType tet1(PointType {1, 0, 0},
                       PointType {1, 1, 0},
                       PointType {0, 1, 0},
                       PointType {1, 0, 1});

  TetrahedronType tet2(PointType {0, 1, 0},
                       PointType {0, 0, 0},
                       PointType {-1, 0, 0},
                       PointType {0, 0, 1});

  PolyhedronType poly = axom::primal::clip(tet1, tet2);
  EXPECT_EQ(0.0, poly.volume());
  EXPECT_EQ(0.0, axom::primal::intersection_volume<double>(tet1, tet2));
}

// Tetrahedrons are the same
TEST(primal_clip, tet_tet_equal)
{
  using namespace Primal3D;
  const double EPS = 1e-4;

  TetrahedronType tet(PointType {1, 0, 0},
                      PointType {1, 1, 0},
                      PointType {0, 1, 0},
                      PointType {1, 0, 1});

  PolyhedronType poly = axom::primal::clip(tet, tet);

  // Expected result should be 0.666 / 4 = 0.1666, volume of tet.
  EXPECT_NEAR(0.1666, poly.volume(), EPS);
  EXPECT_NEAR(0.1666, axom::primal::intersection_volume<double>(tet, tet), EPS);
}

// Tetrahedron is encapsulated inside the other tetrahedron
TEST(primal_clip, tet_tet_encapsulate)
{
  using namespace Primal3D;
  constexpr double EPS = 1e-4;

  TetrahedronType tet1(PointType {1, 0, 0},
                       PointType {1, 1, 0},
                       PointType {0, 1, 0},
                       PointType {1, 0, 1});

  TetrahedronType tet2(PointType {3, 0, 0},
                       PointType {0, 3, 0},
                       PointType {-3, 0, 0},
                       PointType {0, 0, 3});

  PolyhedronType poly = axom::primal::clip(tet1, tet2);

  // Expected result should be 0.666 / 4 = 0.1666, volume of tet.
  EXPECT_NEAR(0.1666, poly.volume(), EPS);
  EXPECT_NEAR(0.1666, axom::primal::intersection_volume<double>(tet1, tet2), EPS);
}

// Half of the tetrahedron is clipped by the other tetrahedron
TEST(primal_clip, tet_tet_half)
{
  using namespace Primal3D;
  constexpr double EPS = 1e-4;

  TetrahedronType tet1(PointType {1, 0, 0},
                       PointType {1, 1, 0},
                       PointType {0, 1, 0},
                       PointType {1, 1, 1});

  TetrahedronType tet2(PointType {0, 0, 0},
                       PointType {1, 0, 0},
                       PointType {1, 1, 0},
                       PointType {1, 1, 1});

  PolyhedronType poly = axom::primal::clip(tet1, tet2);

  // Expected result should be 0.666 / 4 / 2 = 0.0833, volume of tet.
  EXPECT_NEAR(0.0833, poly.volume(), EPS);
  EXPECT_NEAR(0.0833, axom::primal::intersection_volume<double>(tet1, tet2), EPS);
}

// Half of the octahedron is clipped by the tetrahedron.
// Then, we split the octahedron into 8 tetrahedrons,
// and verify the total volume from clipping the 8 split tetrahedrons
// by the starting tetrahedron are the same as with the octahedron.
TEST(primal_clip, tet_tet_clip_split)
{
  using namespace Primal3D;
  constexpr double EPS = 1e-4;
  constexpr bool CHECK_SIGN = true;

  TetrahedronType tet(PointType {0.5, 0.5, 2},
                      PointType {2, -1, 0},
                      PointType {-1, -1, 0},
                      PointType {-1, 2, 0});
  OctahedronType oct(PointType {1, 0, 0},
                     PointType {1, 1, 0},
                     PointType {0, 1, 0},
                     PointType {0, 1, 1},
                     PointType {0, 0, 1},
                     PointType {1, 0, 1});

  PolyhedronType poly = axom::primal::clip(oct, tet);

  // Expected result should be 0.3333, half the volume of oct.
  EXPECT_NEAR(0.3333, poly.volume(), EPS);
  EXPECT_NEAR(0.3333, axom::primal::intersection_volume<double>(oct, tet), EPS);
  EXPECT_NEAR(0.3333, axom::primal::intersection_volume<double>(tet, oct), EPS);

  // Split the octahedron into 8 tetrahedrons
  double tet_volumes = 0.0;

  axom::Array<TetrahedronType> split_tets(8);

  axom::primal::split(oct, split_tets);

  for(int i = 0; i < split_tets.size(); i++)
  {
    tet_volumes +=
      (axom::primal::clip(split_tets[i], tet, EPS, CHECK_SIGN)).volume();
  }

  // Expected result should still be 0.3333
  EXPECT_NEAR(0.3333, tet_volumes, EPS);
}

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  int result = RUN_ALL_TESTS();

  return result;
}

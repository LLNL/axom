// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
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

#include <limits>

namespace Primal3D
{
using PointType = axom::primal::Point<double, 3>;
using VectorType = axom::primal::Vector<double, 3>;
using BoundingBoxType = axom::primal::BoundingBox<double, 3>;
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
void check_oct_tet_clip(double EPS)
{
  using namespace Primal3D;

  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();

  // Save current/default allocator
  const int current_allocator = axom::getDefaultAllocatorID();

  // Determine new allocator (for CUDA policy, set to Unified)
  umpire::Allocator allocator =
    rm.getAllocator(axom::execution_space<ExecPolicy>::allocatorID());

  // Set new default to device
  axom::setDefaultAllocator(allocator.getId());

  // Allocate memory for shapes
  TetrahedronType* tet = axom::allocate<TetrahedronType>(1);
  OctahedronType* oct = axom::allocate<OctahedronType>(1);

  PolyhedronType* res = (axom::execution_space<ExecPolicy>::onDevice()
                           ? axom::allocate<PolyhedronType>(
                               1,
                               rm.getAllocator(umpire::resource::Unified).getId())
                           : axom::allocate<PolyhedronType>(1));

  tet[0] = TetrahedronType(PointType({1, 0, 0}),
                           PointType({1, 1, 0}),
                           PointType({0, 1, 0}),
                           PointType({1, 0, 1}));

  oct[0] = OctahedronType(PointType({1, 0, 0}),
                          PointType({1, 1, 0}),
                          PointType({0, 1, 0}),
                          PointType({0, 1, 1}),
                          PointType({0, 0, 1}),
                          PointType({1, 0, 1}));

  axom::for_all<ExecPolicy>(
    1,
    AXOM_LAMBDA(int i) { res[i] = axom::primal::clip(oct[i], tet[i]); });

  EXPECT_NEAR(0.1666, res[0].volume(), EPS);

  axom::deallocate(tet);
  axom::deallocate(oct);
  axom::deallocate(res);

  axom::setDefaultAllocator(current_allocator);
}

TEST(primal_clip, unit_poly_clip_vertices_sequential)
{
  unit_check_poly_clip<axom::SEQ_EXEC>();
}

TEST(primal_clip, clip_oct_tet_sequential)
{
  const double EPS = 1e-4;
  check_oct_tet_clip<axom::SEQ_EXEC>(EPS);
}

  #ifdef AXOM_USE_OPENMP
TEST(primal_clip, unit_poly_clip_vertices_omp)
{
  unit_check_poly_clip<axom::OMP_EXEC>();
}

TEST(primal_clip, clip_oct_tet_omp)
{
  const double EPS = 1e-4;
  check_oct_tet_clip<axom::OMP_EXEC>(EPS);
}
  #endif /* AXOM_USE_OPENMP */

  #if defined(AXOM_USE_CUDA)
AXOM_CUDA_TEST(primal_clip, unit_poly_clip_vertices_gpu)
{
  unit_check_poly_clip<axom::CUDA_EXEC<256>>();
}

AXOM_CUDA_TEST(primal_clip, clip_oct_tet_cuda)
{
  const double EPS = 1e-4;
  check_oct_tet_clip<axom::CUDA_EXEC<256>>(EPS);
}
  #endif /* AXOM_USE_CUDA */
#endif   /* AXOM_USE_RAJA && AXOM_USE_UMPIRE */

// Tetrahedron does not clip octahedron.
TEST(primal_clip, oct_tet_clip_nonintersect)
{
  using namespace Primal3D;

  TetrahedronType tet(PointType({-1, -1, -1}),
                      PointType({-1, 0, 0}),
                      PointType({-1, -1, 0}),
                      PointType({0, 0, 0}));
  OctahedronType oct(PointType({1, 0, 0}),
                     PointType({1, 1, 0}),
                     PointType({0, 1, 0}),
                     PointType({0, 1, 1}),
                     PointType({0, 0, 1}),
                     PointType({1, 0, 1}));

  PolyhedronType poly = axom::primal::clip(oct, tet);
  EXPECT_EQ(0.0, poly.volume());
}

// Tetrahedron is encapsulated by the octahedron
TEST(primal_clip, oct_tet_clip_encapsulate)
{
  using namespace Primal3D;
  const double EPS = 1e-4;

  TetrahedronType tet(PointType({1, 0, 0}),
                      PointType({1, 1, 0}),
                      PointType({0, 1, 0}),
                      PointType({1, 0, 1}));
  OctahedronType oct(PointType({1, 0, 0}),
                     PointType({1, 1, 0}),
                     PointType({0, 1, 0}),
                     PointType({0, 1, 1}),
                     PointType({0, 0, 1}),
                     PointType({1, 0, 1}));

  PolyhedronType poly = axom::primal::clip(oct, tet);

  // Expected result should be 0.666 / 4 = 0.1666, volume of tet.
  EXPECT_NEAR(0.1666, poly.volume(), EPS);
}

// Octahedron is encapsulated inside the tetrahedron
TEST(primal_clip, oct_tet_clip_encapsulate_inv)
{
  using namespace Primal3D;
  const double EPS = 1e-4;

  TetrahedronType tet(PointType({0, 0, 0}),
                      PointType({0, 2, 0}),
                      PointType({0, 0, 2}),
                      PointType({2, 0, 0}));
  OctahedronType oct(PointType({1, 0, 0}),
                     PointType({1, 1, 0}),
                     PointType({0, 1, 0}),
                     PointType({0, 1, 1}),
                     PointType({0, 0, 1}),
                     PointType({1, 0, 1}));

  PolyhedronType poly = axom::primal::clip(oct, tet);

  // Expected result should be 0.6666, volume of oct.
  EXPECT_NEAR(0.6666, poly.volume(), EPS);
}

// Half of the octahedron is clipped by the tetrahedron
TEST(primal_clip, oct_tet_clip_half)
{
  using namespace Primal3D;
  const double EPS = 1e-4;

  TetrahedronType tet(PointType({0.5, 0.5, 2}),
                      PointType({2, -1, 0}),
                      PointType({-1, -1, 0}),
                      PointType({-1, 2, 0}));
  OctahedronType oct(PointType({1, 0, 0}),
                     PointType({1, 1, 0}),
                     PointType({0, 1, 0}),
                     PointType({0, 1, 1}),
                     PointType({0, 0, 1}),
                     PointType({1, 0, 1}));

  PolyhedronType poly = axom::primal::clip(oct, tet);

  // Expected result should be 0.3333, half the volume of oct.
  EXPECT_NEAR(0.3333, poly.volume(), EPS);
}

// Octahedron is adjacent to tetrahedron
TEST(primal_clip, oct_tet_clip_adjacent)
{
  using namespace Primal3D;

  TetrahedronType tet(PointType({0, -1, 0}),
                      PointType({0, 0, 1}),
                      PointType({1, 0, 1}),
                      PointType({1, 0, 0}));
  OctahedronType oct(PointType({1, 0, 0}),
                     PointType({1, 1, 0}),
                     PointType({0, 1, 0}),
                     PointType({0, 1, 1}),
                     PointType({0, 0, 1}),
                     PointType({1, 0, 1}));

  PolyhedronType poly = axom::primal::clip(oct, tet);

  EXPECT_EQ(0.0, poly.volume());
}

// Tetrahedron clips octahedron at a single vertex
TEST(primal_clip, oct_tet_clip_point)
{
  using namespace Primal3D;

  TetrahedronType tet(PointType({-1, -1, 0}),
                      PointType({-0.5, 0.5, 0}),
                      PointType({0, 0, 2}),
                      PointType({0.5, -0.5, 0}));
  OctahedronType oct(PointType({1, 0, 0}),
                     PointType({1, 1, 0}),
                     PointType({0, 1, 0}),
                     PointType({0, 1, 1}),
                     PointType({0, 0, 1}),
                     PointType({1, 0, 1}));

  PolyhedronType poly = axom::primal::clip(oct, tet);

  EXPECT_EQ(0.0, poly.volume());
}

TEST(primal_clip, oct_tet_clip_special_case_1)
{
  using namespace Primal3D;
  const double EPS = 1e-4;

  TetrahedronType tet(PointType({0.5, 0.5, 0.5}),
                      PointType({1, 1, 0}),
                      PointType({1, 0, 0}),
                      PointType({0.5, 0.5, 0}));

  OctahedronType oct(PointType({0.5, 0.853553, 0.146447}),
                     PointType({0.853553, 0.853553, 0.5}),
                     PointType({0.853553, 0.5, 0.146447}),
                     PointType({1, 0.5, 0.5}),
                     PointType({0.5, 0.5, 0}),
                     PointType({0.5, 1, 0.5}));

  // NOTE: Order of vertices 1,2 and 4,5 are flipped due
  // to winding being opposite of what volume() expects
  PolyhedronType octPoly;
  octPoly.addVertex(PointType({0.5, 0.853553, 0.146447}));
  octPoly.addVertex(PointType({0.853553, 0.5, 0.146447}));
  octPoly.addVertex(PointType({0.853553, 0.853553, 0.5}));
  octPoly.addVertex(PointType({1, 0.5, 0.5}));
  octPoly.addVertex(PointType({0.5, 1, 0.5}));
  octPoly.addVertex(PointType({0.5, 0.5, 0}));

  octPoly.addNeighbors(0, {1, 5, 4, 2});
  octPoly.addNeighbors(1, {0, 2, 3, 5});
  octPoly.addNeighbors(2, {0, 4, 3, 1});
  octPoly.addNeighbors(3, {1, 2, 4, 5});
  octPoly.addNeighbors(4, {0, 5, 3, 2});
  octPoly.addNeighbors(5, {0, 1, 3, 4});

  EXPECT_NEAR(0.0251, octPoly.volume(), EPS);

  PolyhedronType poly = axom::primal::clip(oct, tet);

  EXPECT_NEAR(0.0041, poly.volume(), EPS);
}

TEST(primal_clip, oct_tet_clip_special_case_2)
{
  using namespace Primal3D;
  const double EPS = 1e-4;

  TetrahedronType tet(PointType({0.5, 0.5, 0.5}),
                      PointType({0, 1, 0}),
                      PointType({1, 1, 0}),
                      PointType({0.5, 0.5, 0}));

  OctahedronType oct(PointType({0.5, 0.853553, 0.146447}),
                     PointType({0.853553, 0.853553, 0.5}),
                     PointType({0.853553, 0.5, 0.146447}),
                     PointType({1, 0.5, 0.5}),
                     PointType({0.5, 0.5, 0}),
                     PointType({0.5, 1, 0.5}));

  // NOTE: Order of vertices 1,2 and 4,5 are flipped due
  // to winding being opposite of what volume() expects
  PolyhedronType octPoly;
  octPoly.addVertex(PointType({0.5, 0.853553, 0.146447}));
  octPoly.addVertex(PointType({0.853553, 0.5, 0.146447}));
  octPoly.addVertex(PointType({0.853553, 0.853553, 0.5}));
  octPoly.addVertex(PointType({1, 0.5, 0.5}));
  octPoly.addVertex(PointType({0.5, 1, 0.5}));
  octPoly.addVertex(PointType({0.5, 0.5, 0}));

  octPoly.addNeighbors(0, {1, 5, 4, 2});
  octPoly.addNeighbors(1, {0, 2, 3, 5});
  octPoly.addNeighbors(2, {0, 4, 3, 1});
  octPoly.addNeighbors(3, {1, 2, 4, 5});
  octPoly.addNeighbors(4, {0, 5, 3, 2});
  octPoly.addNeighbors(5, {0, 1, 3, 4});

  EXPECT_NEAR(0.0251, octPoly.volume(), EPS);

  PolyhedronType poly = axom::primal::clip(oct, tet);

  EXPECT_NEAR(0.0041, poly.volume(), EPS);
}

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;  // create & initialize test logger,

  int result = RUN_ALL_TESTS();

  return result;
}

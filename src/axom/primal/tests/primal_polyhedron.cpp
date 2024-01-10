// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core.hpp"
#include "axom/slic.hpp"

#include "axom/primal/geometry/Hexahedron.hpp"
#include "axom/primal/geometry/Octahedron.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Tetrahedron.hpp"
#include "axom/primal/geometry/Polygon.hpp"
#include "axom/primal/geometry/Polyhedron.hpp"

#include "axom/primal/operators/in_polyhedron.hpp"
#include "axom/primal/operators/winding_number.hpp"

#include <math.h>
#include "gtest/gtest.h"

namespace primal = axom::primal;

//------------------------------------------------------------------------------
TEST(primal_polyhedron, polyhedron_empty)
{
  using PolyhedronType = primal::Polyhedron<double, 3>;

  PolyhedronType poly;
  EXPECT_FALSE(poly.isValid());
  EXPECT_FALSE(poly.hasNeighbors());

  EXPECT_NEAR(0., poly.volume(), 1e-12);
}

//------------------------------------------------------------------------------
TEST(primal_polyhedron, polyhedron_unit_cube)
{
  using PolyhedronType = primal::Polyhedron<double, 3>;
  using PointType = primal::Point<double, 3>;

  PolyhedronType poly;
  poly.addVertex({0, 0, 0});
  poly.addVertex({1, 0, 0});
  poly.addVertex({1, 1, 0});
  poly.addVertex({0, 1, 0});
  poly.addVertex({0, 0, 1});
  poly.addVertex({1, 0, 1});
  poly.addVertex({1, 1, 1});
  poly.addVertex({0, 1, 1});

  poly.addNeighbors(poly[0], {1, 4, 3});
  poly.addNeighbors(poly[1], {5, 0, 2});
  poly.addNeighbors(poly[2], {3, 6, 1});
  poly.addNeighbors(poly[3], {7, 2, 0});
  poly.addNeighbors(poly[4], {5, 7, 0});
  poly.addNeighbors(poly[5], {1, 6, 4});
  poly.addNeighbors(poly[6], {2, 7, 5});
  poly.addNeighbors(poly[7], {4, 6, 3});

  EXPECT_EQ(1, poly.volume());

  // Example usage of experimental getFaces() function
  // (input parameters and/or output may change in the future)

  // Euler characteristic to verify provided sizes are adequate
  // (Characteristic = 2 for convex polyhedron)
  constexpr int EULER_CHAR = 2;
  constexpr int NUM_VERTICES = 8;
  constexpr int NUM_EDGES = 12;
  constexpr int NUM_FACES = 6;

  EXPECT_EQ(EULER_CHAR, NUM_VERTICES - NUM_EDGES + NUM_FACES);

  // Cube - 6 faces, 4 vertex indices for each face.
  int faces[NUM_FACES * 4];
  int face_size[NUM_FACES];
  int face_offset[NUM_FACES];
  int face_count;
  poly.getFaces(faces, face_size, face_offset, face_count);

  // Verify we got the expected number of faces
  EXPECT_EQ(face_count, NUM_FACES);

  // Verify face vertex indices
  // clang-format off
  int faces_expect [NUM_FACES * 4] = {0, 1, 5, 4,
                                      0, 4, 7, 3,
                                      0, 3, 2, 1,
                                      1, 2, 6, 5,
                                      2, 3, 7, 6,
                                      4, 5, 6, 7};
  // clang-format on

  for(int f = 0; f < face_count; ++f)
  {
    for(int f_index = face_offset[f]; f_index < face_offset[f] + face_size[f];
        f_index++)
    {
      EXPECT_EQ(faces[f_index], faces_expect[f_index]);
    }
  }

  // Test containment using winding numbers
  bool includeBoundary = true;
  bool useNonzeroRule = true;

  // Use zero numerical tolerances
  bool edge_tol = 0.0;
  bool EPS = 0.0;

  for(double x = -0.5; x <= 1.5; x += 0.1)
  {
    for(double y = -0.5; y <= 1.5; y += 0.1)
    {
      for(double z = -0.5; z <= 1.5; z += 0.1)
      {
        if((x >= 0.0 && x <= 1.0) && (y >= 0.0 && y <= 1.0) &&
           (z >= 0.0 && z <= 1.0))
        {
          EXPECT_TRUE(in_polyhedron(PointType({x, y, z}),
                                    poly,
                                    includeBoundary,
                                    useNonzeroRule,
                                    edge_tol,
                                    EPS));
        }
        else
        {
          EXPECT_FALSE(in_polyhedron(PointType({x, y, z}),
                                     poly,
                                     includeBoundary,
                                     useNonzeroRule,
                                     edge_tol,
                                     EPS));
        }
      }
    }
  }

  // Verify includeBoundary behavior
  EXPECT_TRUE(
    in_polyhedron(poly[0], poly, includeBoundary, useNonzeroRule, edge_tol, EPS));
  EXPECT_FALSE(
    in_polyhedron(poly[0], poly, !includeBoundary, useNonzeroRule, edge_tol, EPS));

  EXPECT_EQ(winding_number(poly[0], poly, includeBoundary, edge_tol, EPS), 1);
  EXPECT_EQ(winding_number(poly[0], poly, !includeBoundary, edge_tol, EPS), 0);
}

//------------------------------------------------------------------------------
TEST(primal_polyhedron, polyhedron_tetrahedron)
{
  using PointType = primal::Point<double, 3>;
  using TetrahedronType = primal::Tetrahedron<double, 3>;
  using PolyhedronType = primal::Polyhedron<double, 3>;

  constexpr double EPS = 1e-4;
  {
    TetrahedronType tet {PointType {1, 1, 1},
                         PointType {-1, 1, -1},
                         PointType {1, -1, -1},
                         PointType {-1, -1, 1}};

    EXPECT_TRUE(tet.signedVolume() > 0.);

    PolyhedronType poly;
    poly.addVertex(tet[0]);
    poly.addVertex(tet[1]);
    poly.addVertex(tet[2]);
    poly.addVertex(tet[3]);

    poly.addNeighbors(0, {1, 3, 2});
    poly.addNeighbors(1, {0, 2, 3});
    poly.addNeighbors(2, {0, 3, 1});
    poly.addNeighbors(3, {0, 1, 2});

    EXPECT_NEAR(tet.volume(), poly.volume(), EPS);
    EXPECT_NEAR(2.6666, poly.volume(), EPS);

    // Test containment using winding numbers
    const bool includeBoundary = true;
    const bool useNonzeroRule = true;

    // Use zero numerical tolerances
    const bool edge_tol = 0.0;
    const bool EPS = 0.0;

    for(double z = -1.5; z < 1.5; z += 0.1)
    {
      if((z >= -1.0) && (z <= 1.0))
      {
        EXPECT_TRUE(in_polyhedron(PointType({0.0, 0.0, z}),
                                  poly,
                                  includeBoundary,
                                  useNonzeroRule,
                                  edge_tol,
                                  EPS));
      }
      else
      {
        EXPECT_FALSE(in_polyhedron(PointType({0.0, 0.0, z}),
                                   poly,
                                   includeBoundary,
                                   useNonzeroRule,
                                   edge_tol,
                                   EPS));
      }
    }
  }

  {
    TetrahedronType tet {PointType {1, 0, 0},
                         PointType {1, 1, 0},
                         PointType {0, 1, 0},
                         PointType {1, 0, 1}};

    EXPECT_TRUE(tet.signedVolume() > 0.);

    PolyhedronType poly;
    poly.addVertex(tet[0]);
    poly.addVertex(tet[1]);
    poly.addVertex(tet[2]);
    poly.addVertex(tet[3]);

    poly.addNeighbors(0, {1, 3, 2});
    poly.addNeighbors(1, {0, 2, 3});
    poly.addNeighbors(2, {0, 3, 1});
    poly.addNeighbors(3, {0, 1, 2});

    EXPECT_NEAR(tet.volume(), poly.volume(), EPS);
    EXPECT_NEAR(0.1666, poly.volume(), EPS);
  }
}

//------------------------------------------------------------------------------
TEST(primal_polyhedron, polyhedron_octahedron)
{
  using PolyhedronType = primal::Polyhedron<double, 3>;

  constexpr double EPS = 1e-4;
  PolyhedronType octA;
  octA.addVertex({0, 0, -1});
  octA.addVertex({1, 0, 0});
  octA.addVertex({0, 1, 0});
  octA.addVertex({0, 0, 1});
  octA.addVertex({-1, 0, 0});
  octA.addVertex({0, -1, 0});

  octA.addNeighbors(octA[0], {1, 5, 4, 2});
  octA.addNeighbors(octA[1], {0, 2, 3, 5});
  octA.addNeighbors(octA[2], {0, 4, 3, 1});
  octA.addNeighbors(octA[3], {1, 2, 4, 5});
  octA.addNeighbors(octA[4], {0, 5, 3, 2});
  octA.addNeighbors(octA[5], {0, 1, 3, 4});

  EXPECT_NEAR(1.3333, octA.volume(), EPS);
  EXPECT_NEAR(1.3333, octA.signedVolume(), EPS);

  PolyhedronType octB;
  octB.addVertex({1, 0, 0});
  octB.addVertex({1, 1, 0});
  octB.addVertex({0, 1, 0});
  octB.addVertex({0, 1, 1});
  octB.addVertex({0, 0, 1});
  octB.addVertex({1, 0, 1});

  octB.addNeighbors(octB[0], {1, 5, 4, 2});
  octB.addNeighbors(octB[1], {0, 2, 3, 5});
  octB.addNeighbors(octB[2], {0, 4, 3, 1});
  octB.addNeighbors(octB[3], {1, 2, 4, 5});
  octB.addNeighbors(octB[4], {0, 5, 3, 2});
  octB.addNeighbors(octB[5], {0, 1, 3, 4});

  EXPECT_NEAR(0.6666, octB.volume(), EPS);
  EXPECT_NEAR(0.6666, octB.signedVolume(), EPS);
}

//------------------------------------------------------------------------------
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)

template <typename ExecSpace>
void check_volume()
{
  const int DIM = 3;
  using PolyhedronType = primal::Polyhedron<double, DIM>;

  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();

  // Save current/default allocator
  const int current_allocator = axom::getDefaultAllocatorID();

  // Determine new allocator (for CUDA or HIP policy, set to Device)
  umpire::Allocator allocator =
    rm.getAllocator(axom::execution_space<ExecSpace>::allocatorID());

  // Set new default to device
  axom::setDefaultAllocator(allocator.getId());

  // Initialize volume
  double* res = axom::allocate<double>(1);

  // Volume on host (for CUDA or HIP)
  double res_host;

  axom::for_all<ExecSpace>(
    1,
    AXOM_LAMBDA(int i) {
      PolyhedronType poly;
      poly.addVertex({0, 0, 0});
      poly.addVertex({1, 0, 0});
      poly.addVertex({1, 1, 0});
      poly.addVertex({0, 1, 0});
      poly.addVertex({0, 0, 1});
      poly.addVertex({1, 0, 1});
      poly.addVertex({1, 1, 1});
      poly.addVertex({0, 1, 1});

      poly.addNeighbors(0, {1, 4, 3});
      poly.addNeighbors(1, {5, 0, 2});
      poly.addNeighbors(2, {3, 6, 1});
      poly.addNeighbors(3, {7, 2, 0});
      poly.addNeighbors(4, {5, 7, 0});
      poly.addNeighbors(5, {1, 6, 4});
      poly.addNeighbors(6, {2, 7, 5});
      poly.addNeighbors(7, {4, 6, 3});

      res[i] = poly.volume();
    });

  axom::copy(&res_host, res, sizeof(double));

  EXPECT_EQ(res_host, 1);

  axom::deallocate(res);

  axom::setDefaultAllocator(current_allocator);
}

TEST(primal_polyhedron, check_volume_sequential)
{
  check_volume<axom::SEQ_EXEC>();
}

  #ifdef AXOM_USE_OPENMP
TEST(primal_polyhedron, check_volume_omp) { check_volume<axom::OMP_EXEC>(); }
  #endif /* AXOM_USE_OPENMP */

  #ifdef AXOM_USE_CUDA
AXOM_CUDA_TEST(primal_polyhedron, check_volume_cuda)
{
  constexpr int BLOCK_SIZE = 256;
  using exec = axom::CUDA_EXEC<BLOCK_SIZE>;

  check_volume<exec>();
}
  #endif /* AXOM_USE_CUDA */

  #ifdef AXOM_USE_HIP
TEST(primal_polyhedron, check_volume_hip)
{
  constexpr int BLOCK_SIZE = 256;
  using exec = axom::HIP_EXEC<BLOCK_SIZE>;

  check_volume<exec>();
}
  #endif /* AXOM_USE_HIP */

#endif /* AXOM_USE_RAJA && AXOM_USE_UMPIRE */

TEST(primal_polyhedron, polyhedron_decomposition)
{
  using PolyhedronType = primal::Polyhedron<double, 3>;
  using PointType = primal::Point<double, 3>;

  constexpr double EPS = 1e-4;

  PolyhedronType poly;
  poly.addVertex({0, 0, 0});
  poly.addVertex({1, 0, 0});
  poly.addVertex({1, 1, 0});
  poly.addVertex({0, 1, 0});
  poly.addVertex({0, 0, 1});
  poly.addVertex({1, 0, 1});
  poly.addVertex({1, 1, 1});
  poly.addVertex({0, 1, 1});

  poly.addNeighbors(poly[0], {1, 4, 3});
  poly.addNeighbors(poly[1], {5, 0, 2});
  poly.addNeighbors(poly[2], {3, 6, 1});
  poly.addNeighbors(poly[3], {7, 2, 0});
  poly.addNeighbors(poly[4], {5, 7, 0});
  poly.addNeighbors(poly[5], {1, 6, 4});
  poly.addNeighbors(poly[6], {2, 7, 5});
  poly.addNeighbors(poly[7], {4, 6, 3});

  // Hex center (hc)
  PointType hc = poly.centroid();

  //Face means (fm)
  PointType fm1 = PointType::midpoint(PointType::midpoint(poly[0], poly[1]),
                                      PointType::midpoint(poly[2], poly[3]));

  PointType fm2 = PointType::midpoint(PointType::midpoint(poly[0], poly[1]),
                                      PointType::midpoint(poly[4], poly[5]));

  PointType fm3 = PointType::midpoint(PointType::midpoint(poly[0], poly[3]),
                                      PointType::midpoint(poly[4], poly[7]));

  PointType fm4 = PointType::midpoint(PointType::midpoint(poly[1], poly[2]),
                                      PointType::midpoint(poly[5], poly[6]));

  PointType fm5 = PointType::midpoint(PointType::midpoint(poly[2], poly[3]),
                                      PointType::midpoint(poly[6], poly[7]));

  PointType fm6 = PointType::midpoint(PointType::midpoint(poly[4], poly[5]),
                                      PointType::midpoint(poly[6], poly[7]));

  PolyhedronType tets[24];

  tets[0].addVertex(hc);
  tets[0].addVertex(poly[1]);
  tets[0].addVertex(poly[0]);
  tets[0].addVertex(fm1);

  tets[1].addVertex(hc);
  tets[1].addVertex(poly[0]);
  tets[1].addVertex(poly[3]);
  tets[1].addVertex(fm1);

  tets[2].addVertex(hc);
  tets[2].addVertex(poly[3]);
  tets[2].addVertex(poly[2]);
  tets[2].addVertex(fm1);

  tets[3].addVertex(hc);
  tets[3].addVertex(poly[2]);
  tets[3].addVertex(poly[1]);
  tets[3].addVertex(fm1);

  tets[4].addVertex(hc);
  tets[4].addVertex(poly[4]);
  tets[4].addVertex(poly[0]);
  tets[4].addVertex(fm2);

  tets[5].addVertex(hc);
  tets[5].addVertex(poly[0]);
  tets[5].addVertex(poly[1]);
  tets[5].addVertex(fm2);

  tets[6].addVertex(hc);
  tets[6].addVertex(poly[1]);
  tets[6].addVertex(poly[5]);
  tets[6].addVertex(fm2);

  tets[7].addVertex(hc);
  tets[7].addVertex(poly[5]);
  tets[7].addVertex(poly[4]);
  tets[7].addVertex(fm2);

  tets[8].addVertex(hc);
  tets[8].addVertex(poly[3]);
  tets[8].addVertex(poly[0]);
  tets[8].addVertex(fm3);

  tets[9].addVertex(hc);
  tets[9].addVertex(poly[0]);
  tets[9].addVertex(poly[4]);
  tets[9].addVertex(fm3);

  tets[10].addVertex(hc);
  tets[10].addVertex(poly[4]);
  tets[10].addVertex(poly[7]);
  tets[10].addVertex(fm3);

  tets[11].addVertex(hc);
  tets[11].addVertex(poly[7]);
  tets[11].addVertex(poly[3]);
  tets[11].addVertex(fm3);

  tets[12].addVertex(hc);
  tets[12].addVertex(poly[5]);
  tets[12].addVertex(poly[1]);
  tets[12].addVertex(fm4);

  tets[13].addVertex(hc);
  tets[13].addVertex(poly[1]);
  tets[13].addVertex(poly[2]);
  tets[13].addVertex(fm4);

  tets[14].addVertex(hc);
  tets[14].addVertex(poly[2]);
  tets[14].addVertex(poly[6]);
  tets[14].addVertex(fm4);

  tets[15].addVertex(hc);
  tets[15].addVertex(poly[6]);
  tets[15].addVertex(poly[5]);
  tets[15].addVertex(fm4);

  tets[16].addVertex(hc);
  tets[16].addVertex(poly[6]);
  tets[16].addVertex(poly[2]);
  tets[16].addVertex(fm5);

  tets[17].addVertex(hc);
  tets[17].addVertex(poly[2]);
  tets[17].addVertex(poly[3]);
  tets[17].addVertex(fm5);

  tets[18].addVertex(hc);
  tets[18].addVertex(poly[3]);
  tets[18].addVertex(poly[7]);
  tets[18].addVertex(fm5);

  tets[19].addVertex(hc);
  tets[19].addVertex(poly[7]);
  tets[19].addVertex(poly[6]);
  tets[19].addVertex(fm5);

  tets[20].addVertex(hc);
  tets[20].addVertex(poly[7]);
  tets[20].addVertex(poly[4]);
  tets[20].addVertex(fm6);

  tets[21].addVertex(hc);
  tets[21].addVertex(poly[4]);
  tets[21].addVertex(poly[5]);
  tets[21].addVertex(fm6);

  tets[22].addVertex(hc);
  tets[22].addVertex(poly[5]);
  tets[22].addVertex(poly[6]);
  tets[22].addVertex(fm6);

  tets[23].addVertex(hc);
  tets[23].addVertex(poly[6]);
  tets[23].addVertex(poly[7]);
  tets[23].addVertex(fm6);

  for(int i = 0; i < 24; i++)
  {
    tets[i].addNeighbors(tets[i][0], {1, 3, 2});
    tets[i].addNeighbors(tets[i][1], {0, 2, 3});
    tets[i].addNeighbors(tets[i][2], {0, 3, 1});
    tets[i].addNeighbors(tets[i][3], {0, 1, 2});
  }

  double sum = 0;
  for(int i = 0; i < 24; i++)
  {
    EXPECT_NEAR(0.0416, tets[i].volume(), EPS);
    sum += tets[i].volume();
  }

  EXPECT_NEAR(1.000, sum, EPS);
}

//------------------------------------------------------------------------------
TEST(primal_polyhedron, polygonal_cone)
{
  namespace numerics = axom::numerics;

  using Polygon2D = primal::Polygon<double, 2>;
  using Polyhedron3D = primal::Polyhedron<double, 3>;
  using Vector3D = primal::Vector<double, 3>;
  using Point3D = primal::Point<double, 3>;
  using MatrixType = numerics::Matrix<double>;

  constexpr double EPS = 1e-8;

  // Lambda to generate a 3D rotation matrix from an angle and axis
  // Formulation from https://en.wikipedia.org/wiki/Rotation_matrix#Axis_and_angle
  auto angleAxisRotMatrix = [](double theta, const Vector3D& axis) -> MatrixType {
    const auto unitized = axis.unitVector();
    const double x = unitized[0], y = unitized[1], z = unitized[2];
    const double c = cos(theta), s = sin(theta), C = 1 - c;

    auto matx = numerics::Matrix<double>::zeros(3, 3);

    matx(0, 0) = x * x * C + c;
    matx(0, 1) = x * y * C - z * s;
    matx(0, 2) = x * z * C + y * s;

    matx(1, 0) = y * x * C + z * s;
    matx(1, 1) = y * y * C + c;
    matx(1, 2) = y * z * C - x * s;

    matx(2, 0) = z * x * C - y * s;
    matx(2, 1) = z * y * C + x * s;
    matx(2, 2) = z * z * C + c;

    return matx;
  };

  // Lambda to rotate the input point using the provided rotation matrix
  auto rotatePoint = [](const MatrixType& matx, const Point3D input) -> Point3D {
    Point3D rotated;
    numerics::matrix_vector_multiply(matx, input.data(), rotated.data());
    return rotated;
  };

  // Create a regular pentagon in the XY plane
  constexpr int N = 5;
  Polygon2D penta(N);
  for(int i = 0; i < N; ++i)
  {
    const double alpha = 2. * M_PI * i / N;
    penta.addVertex({cos(alpha), sin(alpha), 0});
  }
  SLIC_DEBUG(axom::fmt::format("Pentagon w/ signed area {}", penta.signedArea()));

  // Create several rotated cones with a polygonal base.
  // The volume of a cone is 1/3 * base_area * height, so it should equal
  // the area of the polygon when the height is 3.
  for(auto matx : {angleAxisRotMatrix(0., Vector3D {0, 0, 1}),
                   angleAxisRotMatrix(M_PI / 3., Vector3D {0, 1, 0}),
                   angleAxisRotMatrix(M_PI / 2., Vector3D {1, 1, 0}),
                   angleAxisRotMatrix(2 * M_PI / 3., Vector3D {1, 1, 1}),
                   angleAxisRotMatrix(7 * M_PI / 8., Vector3D {1, 0, 0})})
  {
    Polyhedron3D poly;

    // Add vertices for base of pentagonal cone
    for(int i = 0; i < N; ++i)
    {
      poly.addVertex(rotatePoint(matx, Point3D {penta[i][0], penta[i][1], 0}));
    }
    // Add apex of cone; z == 3 for volume to match area of polygonal base
    poly.addVertex(rotatePoint(matx, Point3D {0, 0, 3}));

    // Set up edge adjacencies
    poly.addNeighbors(0, {5, 4, 1});
    poly.addNeighbors(1, {5, 0, 2});
    poly.addNeighbors(2, {5, 1, 3});
    poly.addNeighbors(3, {5, 2, 4});
    poly.addNeighbors(4, {5, 3, 0});
    poly.addNeighbors(5, {0, 1, 2, 3, 4});

    SLIC_DEBUG(axom::fmt::format("Polyhedron {}", poly));

    EXPECT_NEAR(penta.area(), poly.volume(), EPS);

    // Example usage of experimental getFaces() function
    // (input parameters and/or output may change in the future)

    // Euler characteristic to verify provided sizes are adequate
    // (Characteristic = 2 for convex polyhedron)
    constexpr int EULER_CHAR = 2;
    constexpr int NUM_VERTICES = 6;
    constexpr int NUM_EDGES = 10;
    constexpr int NUM_FACES = 6;

    EXPECT_EQ(EULER_CHAR, NUM_VERTICES - NUM_EDGES + NUM_FACES);

    // Pentagonal Cone - 6 faces,
    // 1 face with 5 vertex indices,
    // 5 faces with 3 vertex indices
    constexpr int FACE_IDX_SIZE = 1 * 5 + 5 * 3;
    int faces[FACE_IDX_SIZE];
    int face_size[NUM_FACES];
    int face_offset[NUM_FACES];
    int face_count;
    poly.getFaces(faces, face_size, face_offset, face_count);

    // Verify we got the expected number of faces
    EXPECT_EQ(face_count, NUM_FACES);

    // Verify face vertex indices
    // clang-format off
    int faces_expect [FACE_IDX_SIZE] = {0, 5, 4,
                                        0, 4, 3, 2, 1,
                                        0, 1, 5,
                                        1, 2, 5,
                                        2, 3, 5,
                                        3, 4, 5};
    // clang-format on

    for(int f = 0; f < face_count; ++f)
    {
      for(int f_index = face_offset[f]; f_index < face_offset[f] + face_size[f];
          f_index++)
      {
        EXPECT_EQ(faces[f_index], faces_expect[f_index]);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_polyhedron, polyhedron_from_primitive)
{
  using Hexahedron3D = primal::Hexahedron<double, 3>;
  using Octahedron3D = primal::Octahedron<double, 3>;
  using Polyhedron3D = primal::Polyhedron<double, 3>;
  using Tetrahedron3D = primal::Tetrahedron<double, 3>;
  using Point3D = primal::Point<double, 3>;

  constexpr double EPS = 1e-4;
  constexpr bool CHECK_SIGN = true;

  Polyhedron3D poly;

  // Valid hexahedron
  Hexahedron3D hex(Point3D {0, 0, 1},
                   Point3D {1, 0, 1},
                   Point3D {1, 0, 0},
                   Point3D {0, 0, 0},
                   Point3D {0, 1, 1},
                   Point3D {1, 1, 1},
                   Point3D {1, 1, 0},
                   Point3D {0, 1, 0});

  poly = Polyhedron3D::from_primitive(hex, false);
  EXPECT_NEAR(1.0, poly.volume(), EPS);

  // Check signed volume
  EXPECT_NEAR(1.0, poly.signedVolume(), EPS);
  EXPECT_NEAR(hex.signedVolume(), poly.signedVolume(), EPS);

  // Negative volume
  axom::utilities::swap<Point3D>(hex[1], hex[3]);
  axom::utilities::swap<Point3D>(hex[5], hex[7]);

  poly = Polyhedron3D::from_primitive(hex, false);
  EXPECT_NEAR(-1.0, poly.signedVolume(), EPS);
  EXPECT_NEAR(-1.0, hex.signedVolume(), EPS);

  // Check sign
  poly = Polyhedron3D::from_primitive(hex, CHECK_SIGN);
  EXPECT_NEAR(1.0, poly.signedVolume(), EPS);

  // Valid octahedron
  Octahedron3D oct(Point3D {1, 0, 0},
                   Point3D {1, 1, 0},
                   Point3D {0, 1, 0},
                   Point3D {0, 1, 1},
                   Point3D {0, 0, 1},
                   Point3D {1, 0, 1});

  poly = Polyhedron3D::from_primitive(oct, false);
  EXPECT_NEAR(0.6666, poly.volume(), EPS);

  // Negative volume
  axom::utilities::swap<Point3D>(oct[1], oct[2]);
  axom::utilities::swap<Point3D>(oct[4], oct[5]);

  poly = Polyhedron3D::from_primitive(oct, false);
  EXPECT_NEAR(-0.6666, poly.signedVolume(), EPS);

  // Check sign
  poly = Polyhedron3D::from_primitive(oct, CHECK_SIGN);
  EXPECT_NEAR(0.6666, poly.signedVolume(), EPS);

  // Valid tetrahedron
  Tetrahedron3D tet(Point3D {1, 1, 1},
                    Point3D {-1, 1, -1},
                    Point3D {1, -1, -1},
                    Point3D {-1, -1, 1});

  poly = Polyhedron3D::from_primitive(tet, false);
  EXPECT_NEAR(2.6666, poly.volume(), EPS);

  // Check signed volume
  EXPECT_NEAR(2.6666, poly.signedVolume(), EPS);
  EXPECT_NEAR(tet.signedVolume(), poly.signedVolume(), EPS);

  // Negative volume
  axom::utilities::swap<Point3D>(tet[1], tet[2]);

  poly = Polyhedron3D::from_primitive(tet, false);
  EXPECT_NEAR(-2.6666, poly.signedVolume(), EPS);
  EXPECT_NEAR(tet.signedVolume(), poly.signedVolume(), EPS);

  // Check sign
  poly = Polyhedron3D::from_primitive(tet, CHECK_SIGN);
  EXPECT_NEAR(2.6666, poly.signedVolume(), EPS);
}

//------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger(axom::slic::message::Info);

  result = RUN_ALL_TESTS();

  return result;
}

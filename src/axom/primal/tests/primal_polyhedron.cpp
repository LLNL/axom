// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Tetrahedron.hpp"
#include "axom/primal/geometry/Polyhedron.hpp"

#include "axom/core.hpp"
#include "axom/slic.hpp"

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

  // Determine new allocator (for CUDA or HIP policy, set to Unified)
  umpire::Allocator allocator =
    rm.getAllocator(axom::execution_space<ExecSpace>::allocatorID());

  // Set new default to device
  axom::setDefaultAllocator(allocator.getId());

  // Initialize polyhedron on device,
  // volume results in unified memory to check results on host.
  PolyhedronType* polys = axom::allocate<PolyhedronType>(1);
  bool* res =
    (axom::execution_space<ExecSpace>::onDevice()
       ? axom::allocate<bool>(1,
                              rm.getAllocator(umpire::resource::Unified).getId())
       : axom::allocate<bool>(1));

  polys[0] = PolyhedronType();
  polys[0].addVertex({0, 0, 0});
  polys[0].addVertex({1, 0, 0});
  polys[0].addVertex({1, 1, 0});
  polys[0].addVertex({0, 1, 0});
  polys[0].addVertex({0, 0, 1});
  polys[0].addVertex({1, 0, 1});
  polys[0].addVertex({1, 1, 1});
  polys[0].addVertex({0, 1, 1});

  polys[0].addNeighbors(0, {1, 4, 3});
  polys[0].addNeighbors(1, {5, 0, 2});
  polys[0].addNeighbors(2, {3, 6, 1});
  polys[0].addNeighbors(3, {7, 2, 0});
  polys[0].addNeighbors(4, {5, 7, 0});
  polys[0].addNeighbors(5, {1, 6, 4});
  polys[0].addNeighbors(6, {2, 7, 5});
  polys[0].addNeighbors(7, {4, 6, 3});

  axom::for_all<ExecSpace>(
    1,
    AXOM_LAMBDA(int i) { res[i] = polys[i].volume(); });

  EXPECT_EQ(res[0], 1);

  axom::deallocate(polys);
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

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger(axom::slic::message::Info);

  result = RUN_ALL_TESTS();

  return result;
}

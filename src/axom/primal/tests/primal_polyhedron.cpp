// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/primal/geometry/Polyhedron.hpp"

#include "axom/core/execution/for_all.hpp"
#include "axom/core/execution/execution_space.hpp"
#include "axom/core/memory_management.hpp"

#include "fmt/fmt.hpp"
#include "axom/slic/interface/slic.hpp"

using namespace axom;

//------------------------------------------------------------------------------
TEST(primal_polyhedron, polyhedron_empty)
{
  using PolyhedronType = primal::Polyhedron<double, 3>;

  PolyhedronType poly;
  EXPECT_FALSE(poly.isValid());
  EXPECT_FALSE(poly.hasNeighbors());
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
}

//------------------------------------------------------------------------------
TEST(primal_polyhedron, polyhedron_tetrahedron)
{
  using PolyhedronType = primal::Polyhedron<double, 3>;

  static const double EPS = 1e-4;
  PolyhedronType poly;
  poly.addVertex({1, 1, 1});
  poly.addVertex({-1, 1, -1});
  poly.addVertex({1, -1, -1});
  poly.addVertex({-1, -1, 1});

  poly.addNeighbors(poly[0], {1, 3, 2});
  poly.addNeighbors(poly[1], {0, 2, 3});
  poly.addNeighbors(poly[2], {0, 3, 1});
  poly.addNeighbors(poly[3], {0, 1, 2});

  EXPECT_NEAR(2.6666, poly.volume(), EPS);

  PolyhedronType polyB;
  polyB.addVertex({1, 0, 0});
  polyB.addVertex({1, 1, 0});
  polyB.addVertex({0, 1, 0});
  polyB.addVertex({1, 0, 1});

  polyB.addNeighbors(polyB[0], {1, 3, 2});
  polyB.addNeighbors(polyB[1], {0, 2, 3});
  polyB.addNeighbors(polyB[2], {0, 3, 1});
  polyB.addNeighbors(polyB[3], {0, 1, 2});

  EXPECT_NEAR(0.1666, polyB.volume(), EPS);
}

//------------------------------------------------------------------------------
TEST(primal_polyhedron, polyhedron_octahedron)
{
  using PolyhedronType = primal::Polyhedron<double, 3>;

  static const double EPS = 1e-4;
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
  using PointType = primal::Point<double, DIM>;
  using PolyhedronType = primal::Polyhedron<double, DIM>;

  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();

  // Save current/default allocator
  const int current_allocator = axom::getDefaultAllocatorID();

  // Determine new allocator (for CUDA policy, set to Unified)
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
    AXOM_LAMBDA(int i) { res[0] = polys[0].volume(); });

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

#endif /* AXOM_USE_RAJA && AXOM_USE_UMPIRE */

//------------------------------------------------------------------------------
#include "axom/slic/core/SimpleLogger.hpp"
using axom::slic::SimpleLogger;

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  SimpleLogger logger;  // create & initialize test logger,
  axom::slic::setLoggingMsgLevel(axom::slic::message::Info);

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}

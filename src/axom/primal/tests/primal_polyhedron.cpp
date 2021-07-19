// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/primal/geometry/Polyhedron.hpp"

#include "fmt/fmt.hpp"
#include "axom/slic/interface/slic.hpp"

using namespace axom;

typedef primal::Polyhedron<double, 3> PolyhedronType;

//------------------------------------------------------------------------------
TEST(primal_polyhedron, polyhedron_empty)
{
  PolyhedronType poly;
  EXPECT_FALSE(poly.isValid());
  EXPECT_FALSE(poly.hasNeighbors());
}

//------------------------------------------------------------------------------
TEST(primal_polyhedron, polyhedron_unit_cube)
{
  PolyhedronType poly;
  poly.addVertex({0, 0, 0});
  poly.addVertex({1, 0, 0});
  poly.addVertex({1, 1, 0});
  poly.addVertex({0, 1, 0});
  poly.addVertex({0, 0, 1});
  poly.addVertex({1, 0, 1});
  poly.addVertex({1, 1, 1});
  poly.addVertex({0, 1, 1});

  poly.addNeighbor({1, 4, 3});
  poly.addNeighbor({5, 0, 2});
  poly.addNeighbor({3, 6, 1});
  poly.addNeighbor({7, 2, 0});
  poly.addNeighbor({5, 7, 0});
  poly.addNeighbor({1, 6, 4});
  poly.addNeighbor({2, 7, 5});
  poly.addNeighbor({4, 6, 3});

  EXPECT_EQ(1, poly.volume());
}

//------------------------------------------------------------------------------
TEST(primal_polyhedron, polyhedron_tetrahedron)
{
  static const double EPS = 1e-4;
  PolyhedronType poly;
  poly.addVertex({1, 1, 1});
  poly.addVertex({-1, 1, -1});
  poly.addVertex({1, -1, -1});
  poly.addVertex({-1, -1, 1});

  poly.addNeighbor({1, 3, 2});
  poly.addNeighbor({0, 2, 3});
  poly.addNeighbor({0, 3, 1});
  poly.addNeighbor({0, 1, 2});

  EXPECT_NEAR(2.6666, poly.volume(), EPS);

  PolyhedronType polyB;
  polyB.addVertex({1, 0, 0});
  polyB.addVertex({1, 1, 0});
  polyB.addVertex({0, 1, 0});
  polyB.addVertex({1, 0, 1});

  polyB.addNeighbor({1, 3, 2});
  polyB.addNeighbor({0, 2, 3});
  polyB.addNeighbor({0, 3, 1});
  polyB.addNeighbor({0, 1, 2});

  EXPECT_NEAR(0.1666, polyB.volume(), EPS);
}

//------------------------------------------------------------------------------
TEST(primal_polyhedron, polyhedron_octahedron)
{
  static const double EPS = 1e-4;
  PolyhedronType octA;
  octA.addVertex({0, 0, -1});
  octA.addVertex({1, 0, 0});
  octA.addVertex({0, 1, 0});
  octA.addVertex({0, 0, 1});
  octA.addVertex({-1, 0, 0});
  octA.addVertex({0, -1, 0});

  octA.addNeighbor({1, 5, 4, 2});
  octA.addNeighbor({0, 2, 3, 5});
  octA.addNeighbor({0, 4, 3, 1});
  octA.addNeighbor({1, 2, 4, 5});
  octA.addNeighbor({0, 5, 3, 2});
  octA.addNeighbor({0, 1, 3, 4});

  EXPECT_NEAR(1.3333, octA.volume(), EPS);

  PolyhedronType octB;
  octB.addVertex({1, 0, 0});
  octB.addVertex({1, 1, 0});
  octB.addVertex({0, 1, 0});
  octB.addVertex({0, 1, 1});
  octB.addVertex({0, 0, 1});
  octB.addVertex({1, 0, 1});

  octB.addNeighbor({1, 5, 4, 2});
  octB.addNeighbor({0, 2, 3, 5});
  octB.addNeighbor({0, 4, 3, 1});
  octB.addNeighbor({1, 2, 4, 5});
  octB.addNeighbor({0, 5, 3, 2});
  octB.addNeighbor({0, 1, 3, 4});

  EXPECT_NEAR(0.6666, octB.volume(), EPS);
}

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

// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"
#include "axom/mint.hpp"

#include "axom/quest/InOutOctree.hpp"

#include "quest_test_utilities.hpp"

#include <cstdlib>
#include <limits>

namespace
{
const int DIM = 2;
using Octree2D = axom::quest::InOutOctree<DIM>;
using GeometricBoundingBox = Octree2D::GeometricBoundingBox;
using SpacePt = Octree2D::SpacePt;

using SpaceVector = Octree2D::SpaceVector;
using GridPt = Octree2D::GridPt;
using BlockIndex = Octree2D::BlockIndex;
}  // namespace

/// Returns a SpacePt corresponding to the given vertex id \a vIdx  in \a mesh
SpacePt getVertex(axom::mint::Mesh*& mesh, int vIdx)
{
  SpacePt pt;
  mesh->getNode(vIdx, pt.data());

  return pt;
}

GeometricBoundingBox computeBoundingBox(axom::mint::Mesh*& mesh)
{
  GeometricBoundingBox bbox;
  for(int i = 0; i < mesh->getNumberOfNodes(); ++i)
  {
    bbox.addPoint(getVertex(mesh, i));
  }

  return bbox;
}

TEST(quest_inout_quadtree, triangle_boundary_mesh)
{
  SLIC_INFO("*** Exercises InOutOctree queries for several thresholds.\n");

  namespace mint = axom::mint;
  namespace quest = axom::quest;

  std::vector<double> thresholds = {1E-9, 1E-5, 1E-1, 0.};

  for(auto thresh : thresholds)
  {
    // create a simple mesh on the boundary of an equilateral triangle
    mint::Mesh* mesh = [=]() {
      auto* mesh =
        new mint::UnstructuredMesh<mint::SINGLE_SHAPE>(DIM, mint::SEGMENT);
      mesh->appendNode(1, 1);
      mesh->appendNode(4, 1);
      mesh->appendNode(2.5, 3 * sqrt(3) / 2.);

      axom::IndexType cell[4] = {0, 1, 2, 0};
      mesh->appendCell(&cell[0]);
      mesh->appendCell(&cell[1]);
      mesh->appendCell(&cell[2]);

      return mesh;
    }();

    GeometricBoundingBox bbox = computeBoundingBox(mesh);

    Octree2D octree(bbox, mesh);
    octree.setVertexWeldThreshold(thresh);

    octree.generateIndex();

    SpacePt queryInside = quest::utilities::getCentroid(getVertex(mesh, 0),
                                                        getVertex(mesh, 1),
                                                        getVertex(mesh, 2));
    SpacePt queryOutside = SpacePt(2. * bbox.getMax().array());

    EXPECT_TRUE(octree.within(queryInside));
    EXPECT_FALSE(octree.within(queryOutside));

    delete mesh;
  }
}

TEST(quest_inout_quadtree, circle_mesh)
{
  SLIC_INFO("*** Exercises InOutOctree over the boundary of a circle.\n");

  namespace mint = axom::mint;
  namespace quest = axom::quest;

  double radius = 1.;
  int num_segments = 100;
  mint::Mesh* mesh = quest::utilities::make_circle_mesh_2d(radius, num_segments);
  mint::write_vtk(mesh, "circle_mesh.vtk");

  GeometricBoundingBox bbox = computeBoundingBox(mesh).scale(1.2);

  Octree2D octree(bbox, mesh);
  octree.generateIndex();

  SpacePt queryInside = SpacePt {0, 0};
  SpacePt queryOutside = SpacePt {2. * radius, 2. * radius};

  EXPECT_TRUE(octree.within(queryInside));
  EXPECT_FALSE(octree.within(queryOutside));

  delete mesh;
}

//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;  // create & initialize test logger,

  int result = RUN_ALL_TESTS();
  return result;
}

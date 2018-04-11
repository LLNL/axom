/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#include "gtest/gtest.h"

#include "slic/slic.hpp"

#include "primal/Point.hpp"
#include "mint/UnstructuredMesh.hpp"
#include "quest/MeshTester.hpp"
#include "quest_test_utilities.hpp"

namespace
{

static const int DIM = 3;
static const double EPS = 1e-6;

typedef axom::mint::UnstructuredMesh< MINT_TRIANGLE > TriangleMesh;
typedef axom::primal::Point<double, 3> Point3;

/*! Insert a vertex with coordinates (x,y,z) into \a mesh  */
void insertVertex(TriangleMesh* mesh, double x, double y, double z)
{
  double coords[3] = {x,y,z};
  mesh->insertNode(coords);
}

/*! Insert a triangle with vertex indices (v1,v2,v3) into \a mesh  */
void insertTriangle(TriangleMesh* mesh, int v1, int v2, int v3)
{
  int indices[3] = {v1,v2,v3};
  mesh->insertCell(indices, MINT_TRIANGLE, 3);
}

}

//------------------------------------------------------------------------------
TEST( quest_vertex_weld, emptyMesh)
{
  SLIC_INFO("*** Tests welding function on an empty triangle mesh");

  TriangleMesh* mesh = new TriangleMesh(DIM);

  EXPECT_EQ(0, mesh->getMeshNumberOfNodes());
  EXPECT_EQ(0, mesh->getMeshNumberOfCells());

  axom::quest::weldTriMeshVertices(&mesh, EPS);

  EXPECT_EQ(0, mesh->getMeshNumberOfNodes());
  EXPECT_EQ(0, mesh->getMeshNumberOfCells());

  delete mesh;
  mesh = AXOM_NULLPTR;
}

//------------------------------------------------------------------------------
TEST( quest_vertex_weld, onlyVertices)
{
  SLIC_INFO(
    "*** Tests welding function on a triangle mesh"
    << " with vertices but no triangles.");

  TriangleMesh* mesh = new TriangleMesh(DIM);
  insertVertex(mesh, 1, 0, 0);
  insertVertex(mesh, 1, 1, 0);
  insertVertex(mesh, 1, 1, 1);

  EXPECT_EQ(3, mesh->getMeshNumberOfNodes());
  EXPECT_EQ(0, mesh->getMeshNumberOfCells());

  axom::quest::weldTriMeshVertices(&mesh, EPS);

  EXPECT_EQ(3, mesh->getMeshNumberOfNodes());
  EXPECT_EQ(0, mesh->getMeshNumberOfCells());

  delete mesh;
  mesh = AXOM_NULLPTR;
}

//------------------------------------------------------------------------------
TEST( quest_vertex_weld, oneTriangle)
{
  SLIC_INFO("*** Tests welding function on a triangle mesh with one triangle.");

  TriangleMesh* mesh = new TriangleMesh(DIM);
  insertVertex(mesh, 1, 0, 0);
  insertVertex(mesh, 1, 1, 0);
  insertVertex(mesh, 1, 1, 1);

  EXPECT_EQ(3, mesh->getMeshNumberOfNodes());

  insertTriangle(mesh, 0, 1, 2);
  EXPECT_EQ(1, mesh->getMeshNumberOfCells());

  axom::quest::weldTriMeshVertices(&mesh, EPS);

  EXPECT_EQ(3, mesh->getMeshNumberOfNodes());
  EXPECT_EQ(1, mesh->getMeshNumberOfCells());

  delete mesh;
  mesh = AXOM_NULLPTR;
}

//------------------------------------------------------------------------------
TEST( quest_vertex_weld, degenerateTriangle)
{
  SLIC_INFO(
    "*** Tests welding function on a triangle mesh degenerate triangles.");

  TriangleMesh* mesh = new TriangleMesh(DIM);
  insertVertex(mesh, 1, 0, 0);
  insertVertex(mesh, 1, 1, 0);
  insertVertex(mesh, 1, 1, 1);
  insertVertex(mesh, 1.01, 1.01, 1.01); // should be welded

  EXPECT_EQ(4, mesh->getMeshNumberOfNodes());

  insertTriangle(mesh, 0, 1, 2);
  insertTriangle(mesh, 0, 1, 1); // degenerate: only two distinct vertices
  insertTriangle(mesh, 1, 1, 1); // degenerate: one one distinct vertex
  insertTriangle(mesh, 1, 2, 3); // degenerate: verts 2 and 3 should be welded
  EXPECT_EQ(4, mesh->getMeshNumberOfCells());

  const double eps = .1;
  axom::quest::weldTriMeshVertices(&mesh, eps);

  EXPECT_EQ(3, mesh->getMeshNumberOfNodes());
  EXPECT_EQ(1, mesh->getMeshNumberOfCells());

  delete mesh;
  mesh = AXOM_NULLPTR;
}

//------------------------------------------------------------------------------
TEST( quest_vertex_weld, vertexAdjacentTrianglePair)
{
  SLIC_INFO(
    "*** Tests welding function on a triangle mesh"
    << " with two triangles adjacent along a vertex.");

  TriangleMesh* mesh = new TriangleMesh(DIM);
  insertVertex(mesh, 1, 0, 0);
  insertVertex(mesh, 1, 1, 0);
  insertVertex(mesh, 1, 1, 1);

  insertVertex(mesh, 0, 1, 0);
  insertVertex(mesh, 0, 1, 1);
  insertVertex(mesh, 1, 1, 1); // Duplicate

  EXPECT_EQ(6, mesh->getMeshNumberOfNodes());

  insertTriangle(mesh, 0, 1, 2);
  insertTriangle(mesh, 3, 4, 5);
  EXPECT_EQ(2, mesh->getMeshNumberOfCells());

  axom::quest::weldTriMeshVertices(&mesh, EPS);

  // Duplicated vertices should have been removed
  EXPECT_EQ(5, mesh->getMeshNumberOfNodes());
  EXPECT_EQ(2, mesh->getMeshNumberOfCells());

  delete mesh;
  mesh = AXOM_NULLPTR;
}

//------------------------------------------------------------------------------
TEST( quest_vertex_weld, edgeAdjacentTrianglePair)
{
  SLIC_INFO(
    "*** Tests welding function on a triangle mesh"
    << " with two triangles adjacent along an edge.");

  TriangleMesh* mesh = new TriangleMesh(DIM);
  insertVertex(mesh, 1, 0, 0);
  insertVertex(mesh, 1, 1, 0);
  insertVertex(mesh, 1, 1, 1);

  insertVertex(mesh, 0, 1, 0);
  insertVertex(mesh, 1, 1, 0); // Duplicate
  insertVertex(mesh, 1, 1, 1); // Duplicate

  EXPECT_EQ(6, mesh->getMeshNumberOfNodes());

  insertTriangle(mesh, 0, 1, 2);
  insertTriangle(mesh, 3, 4, 5);
  EXPECT_EQ(2, mesh->getMeshNumberOfCells());

  axom::quest::weldTriMeshVertices(&mesh, EPS);

  // Duplicated vertices should have been removed
  EXPECT_EQ(4, mesh->getMeshNumberOfNodes());
  EXPECT_EQ(2, mesh->getMeshNumberOfCells());

  delete mesh;
  mesh = AXOM_NULLPTR;
}

//------------------------------------------------------------------------------
TEST( quest_vertex_weld, fuzzWeld)
{
  SLIC_INFO(
    "*** Tests welding function on a triangle mesh"
    << " with a pair of triangles adjacent along an edge"
    << " whose vertices are welded");

  TriangleMesh* mesh = new TriangleMesh(DIM);
  insertVertex(mesh, 1, 0, 0);
  insertVertex(mesh, 1, 1, 0);
  insertVertex(mesh, 1, 1, 1);
  insertVertex(mesh, 0.99, 0.99, -0.01); // should weld with 2nd vert
  insertVertex(mesh, 1.01, 1.01, 1.01);  // should weld with 3rd vert
  insertVertex(mesh, 0, 1, 0);

  EXPECT_EQ(6, mesh->getMeshNumberOfNodes());

  insertTriangle(mesh, 0, 1, 2);
  insertTriangle(mesh, 3, 4, 5);
  EXPECT_EQ(2, mesh->getMeshNumberOfCells());

  const double eps = .1;
  axom::quest::weldTriMeshVertices(&mesh, eps);

  EXPECT_EQ(4, mesh->getMeshNumberOfNodes());
  EXPECT_EQ(2, mesh->getMeshNumberOfCells());

  delete mesh;
  mesh = AXOM_NULLPTR;
}

//------------------------------------------------------------------------------
TEST( quest_vertex_weld, disconnectedTrianglePair)
{
  SLIC_INFO(
    "*** Tests welding function on a triangle mesh"
    << " with a pair of disconnected triangles.");

  TriangleMesh* mesh = new TriangleMesh(DIM);
  insertVertex(mesh, 1, 0, 0);    // verts for first triangle
  insertVertex(mesh, 1, 1, 0);
  insertVertex(mesh, 1, 1, 1);
  insertVertex(mesh,  0, -1,  0); // verts for second triangle
  insertVertex(mesh, -1, -1,  0);
  insertVertex(mesh, -1, -1, -1);

  EXPECT_EQ(6, mesh->getMeshNumberOfNodes());

  insertTriangle(mesh, 0, 1, 2);
  insertTriangle(mesh, 3, 4, 5);
  EXPECT_EQ(2, mesh->getMeshNumberOfCells());

  axom::quest::weldTriMeshVertices(&mesh, EPS);

  // Duplicated vertices should have been removed
  EXPECT_EQ(6, mesh->getMeshNumberOfNodes());
  EXPECT_EQ(2, mesh->getMeshNumberOfCells());

  delete mesh;
  mesh = AXOM_NULLPTR;
}

//------------------------------------------------------------------------------
TEST( quest_vertex_weld, indexedTetrahedron)
{
  SLIC_INFO("*** Tests welding function on an indexed tetrahedron");

  // Get the tetrahedron mesh.  It is already an indexed mesh.
  TriangleMesh* mesh = static_cast<TriangleMesh*>(
    axom::quest::utilities::make_tetrahedron_mesh());

  const int NV = 4;
  const int NT = 4;

  EXPECT_EQ(NV, mesh->getMeshNumberOfNodes());
  EXPECT_EQ(NT, mesh->getMeshNumberOfCells());

  axom::quest::weldTriMeshVertices(&mesh, EPS);

  // welding shouldn't change anything on this mesh
  EXPECT_EQ(NV, mesh->getMeshNumberOfNodes());
  EXPECT_EQ(NT, mesh->getMeshNumberOfCells());

  delete mesh;
  mesh = AXOM_NULLPTR;
}

//------------------------------------------------------------------------------
TEST( quest_vertex_weld, indexedOctahedron)
{
  SLIC_INFO("*** Tests welding function on an indexed octahedron");

  // Get the octahedron mesh.  It is already an indexed mesh.
  TriangleMesh* mesh = static_cast<TriangleMesh*>(
    axom::quest::utilities::make_octahedron_mesh());

  const int NV = 6;
  const int NT = 8;

  EXPECT_EQ(NV, mesh->getMeshNumberOfNodes());
  EXPECT_EQ(NT, mesh->getMeshNumberOfCells());

  axom::quest::weldTriMeshVertices(&mesh, EPS);

  // welding shouldn't change anything on this mesh
  EXPECT_EQ(NV, mesh->getMeshNumberOfNodes());
  EXPECT_EQ(NT, mesh->getMeshNumberOfCells());

  delete mesh;
  mesh = AXOM_NULLPTR;
}


//----------------------------------------------------------------------
//----------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}

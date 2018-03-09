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

// Axom includes
#include "quest/MeshTester.hpp"

#include "primal/BoundingBox.hpp"
#include "primal/intersect.hpp"
#include "primal/Point.hpp"
#include "primal/Triangle.hpp"
#include "primal/UniformGrid.hpp"

// C++ includes
#include <cmath>
#include <algorithm>
#include <map>
#include <set>

namespace axom
{
namespace quest
{

typedef mint::UnstructuredMesh< MINT_TRIANGLE > TriangleMesh;
typedef primal::Triangle<double, 3> Triangle3;

typedef primal::Point<double, 3> Point3;
typedef primal::BoundingBox<double, 3> SpatialBoundingBox;
typedef primal::UniformGrid<int, 3> UniformGrid3;
typedef primal::Vector<double, 3> Vector3;
typedef primal::Segment<double, 3> Segment3;

inline SpatialBoundingBox compute_bounds(mint::Mesh* mesh)
{
  SLIC_ASSERT( mesh != AXOM_NULLPTR );

  SpatialBoundingBox meshBB;
  Point3 pt;

  for ( int i=0 ; i < mesh->getMeshNumberOfNodes() ; ++i )
  {
    mesh->getMeshNode( i, pt.data() );
    meshBB.addPoint( pt );
  } // END for all nodes

  SLIC_ASSERT( meshBB.isValid() );

  return meshBB;
}

inline SpatialBoundingBox compute_bounds(const Triangle3 & tri)
{
  SpatialBoundingBox triBB;
  triBB.addPoint(tri[0]);
  triBB.addPoint(tri[1]);
  triBB.addPoint(tri[2]);

  SLIC_ASSERT( triBB.isValid() );

  return triBB;
}

inline Triangle3 getMeshTriangle(int i, mint::Mesh* surface_mesh)
{
  SLIC_ASSERT(surface_mesh->getMeshNumberOfCellNodes(i) == 3);
  primal::Point<int, 3> triCell;
  Triangle3 tri;
  surface_mesh->getMeshCell(i, triCell.data());

  surface_mesh->getMeshNode(triCell[0], tri[0].data());
  surface_mesh->getMeshNode(triCell[1], tri[1].data());
  surface_mesh->getMeshNode(triCell[2], tri[2].data());

  return tri;
}

/*!
 * \brief Implementation for quest::findTriMeshIntersections
 *
 * \see findTriMeshIntersections
 */
void findTriMeshIntersections(
  TriangleMesh* surface_mesh,
  std::vector<std::pair<int, int> > & intersections,
  std::vector<int> & degenerateIndices,
  int spatialIndexResolution)
{
  Triangle3 t1 = Triangle3();
  Triangle3 t2 = Triangle3();
  SLIC_INFO("Running mesh_tester with UniformGrid index");

  // Create a bounding box around mesh to find the minimum point
  SpatialBoundingBox meshBB = compute_bounds(surface_mesh);
  const Point3 & minBBPt = meshBB.getMin();
  const Point3 & maxBBPt = meshBB.getMax();

  const int ncells = surface_mesh->getMeshNumberOfCells();

  // find the specified resolution.  If we're passed a number less than one,
  // use the cube root of the number of triangles.
  if (spatialIndexResolution < 1)
  {
    spatialIndexResolution = (int)(1 + std::pow(ncells, 1/3.));
  }
  int resolutions[3]=
  {spatialIndexResolution, spatialIndexResolution, spatialIndexResolution};

  SLIC_INFO("Building UniformGrid index...");
  UniformGrid3 ugrid(minBBPt.data(), maxBBPt.data(), resolutions);
  std::vector<int> nondegenerateIndices;
  nondegenerateIndices.reserve(ncells);

  for (int i=0 ; i < ncells ; i++)
  {
    t1=getMeshTriangle(i, surface_mesh);

    if (t1.degenerate())
    {
      degenerateIndices.push_back(i);
    }
    else
    {
      nondegenerateIndices.push_back(i);

      SpatialBoundingBox triBB = compute_bounds(t1);
      ugrid.insert(triBB, i);
    }
  }


  // Iterate through triangle indices *idx.
  // Check against each other triangle with index greater than the index *idx
  // that also shares a UniformGrid bin.
  SLIC_INFO("Checking mesh with a total of " << ncells << " cells.");

  std::vector<int>::iterator
    idx = nondegenerateIndices.begin(),
    ndgend = nondegenerateIndices.end();
  for ( ; idx != ndgend ; ++idx)
  {
    // Retrieve the triangle at *idx and construct a bounding box around it
    t1 = getMeshTriangle(*idx, surface_mesh);
    SpatialBoundingBox triBB2 = compute_bounds(t1);

    // Get a list of all triangles in bins this triangle will touch,
    // whose indices are greater than this triangle's index
    std::vector<int> neighborTriangles;
    const std::vector<int> binsToCheck = ugrid.getBinsForBbox(triBB2);
    size_t checkcount = binsToCheck.size();
    for (size_t curbin = 0 ; curbin < checkcount ; ++curbin)
    {
      std::vector<int> ntlist = ugrid.getBinContents(binsToCheck[curbin]);
      std::vector<int>::iterator ntlit = ntlist.begin(), ntlend = ntlist.end();
      for ( ; ntlit != ntlend ; ++ntlit)
      {
        if (*ntlit > *idx)
        {
          neighborTriangles.push_back(*ntlit);
        }
      }
    }

    std::sort(neighborTriangles.begin(), neighborTriangles.end());
    std::vector<int>::iterator nend =
      std::unique(neighborTriangles.begin(), neighborTriangles.end());
    std::vector<int>::iterator nit = neighborTriangles.begin();

    // test any remaining neighbor tris for intersection
    while (nit != nend)
    {
      t2 = getMeshTriangle(*nit, surface_mesh);
      if (primal::intersect(t1, t2))
      {
        intersections.push_back(std::make_pair(*idx, *nit));
      }
      ++nit;
    }
  }
}

} // end namespace quest
} // end namespace axom

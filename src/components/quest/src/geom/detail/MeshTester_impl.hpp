#ifndef MESH_TESTER_IMPL_HPP_
#define MESH_TESTER_IMPL_HPP_

// Axom includes
#include "primal/BoundingBox.hpp"
#include "primal/Point.hpp"
#include "primal/Triangle.hpp"
#include "primal/intersect.hpp"
#include "primal/UniformGrid.hpp"

#include "mint/UnstructuredMesh.hpp"

// C++ includes
#include <algorithm>
#include <map>
#include <set>
#include <cmath>

namespace axom {
namespace quest {
namespace detail {

typedef mint::UnstructuredMesh< MINT_TRIANGLE > TriangleMesh;
typedef primal::Triangle<double, 3> Triangle3;

typedef primal::Point<double, 3> Point3;
typedef primal::BoundingBox<double, 3> SpatialBoundingBox;
typedef primal::UniformGrid<int, 3> UniformGrid3;
typedef primal::Vector<double, 3> Vector3;
typedef primal::Segment<double, 3> Segment3;

SpatialBoundingBox compute_bounds(mint::Mesh* mesh)
{
  SLIC_ASSERT( mesh != AXOM_NULLPTR );

  SpatialBoundingBox meshBB;
  Point3 pt;

  for ( int i=0; i < mesh->getMeshNumberOfNodes(); ++i )
  {
    mesh->getMeshNode( i, pt.data() );
    meshBB.addPoint( pt );
  } // END for all nodes

  SLIC_ASSERT( meshBB.isValid() );

  return meshBB;
}

Triangle3 getMeshTriangle(int i, mint::Mesh* surface_mesh)
{
  SLIC_ASSERT(surface_mesh->getMeshNumberOfCellNodes(i) == 3);
  primal::Point<int, 3> triCell;
  surface_mesh->getMeshCell( i, triCell.data());
  primal::Point< double,3 > A1;
  primal::Point< double,3 > B1;
  primal::Point< double,3 > C1;

  surface_mesh->getMeshNode(triCell[0], A1.data());
  surface_mesh->getMeshNode(triCell[1], B1.data());
  surface_mesh->getMeshNode(triCell[2], C1.data());
  Triangle3 t1 = Triangle3(A1,B1,C1);

  return t1;
}

bool checkTT(Triangle3& t1, Triangle3& t2)
{
  if (t2.degenerate()) return false;

  if (primal::intersect(t1, t2)) {
    return true;
  }
  return false;
}

std::vector<std::pair<int, int> > 
findTriMeshIntersections_impl(mint::Mesh* surface_mesh,
			      std::vector<int> & degenerateIndices,
			      int spatialIndexResolution)
{
  std::vector<std::pair<int, int> > retval;
  std::set<int> seen, hit;

  Triangle3 t1 = Triangle3();
  Triangle3 t2 = Triangle3();
  SLIC_INFO("Running mesh_tester with UniformGrid index");

  // Create a bounding box around mesh to find the minimum point
  SpatialBoundingBox meshBB = compute_bounds(surface_mesh);
  const Point3 minBBPt = meshBB.getMin();
  const Point3 maxBBPt = meshBB.getMax();

  const int ncells = surface_mesh->getMeshNumberOfCells();

  // find the specified resolution.  If we're passed a number less than one,
  // use the cube root of the number of triangles.
  if (spatialIndexResolution < 1) {
    spatialIndexResolution = (int)(1 + std::pow(ncells, 1/3.));
  }
  int resolutions[3]={spatialIndexResolution, spatialIndexResolution, spatialIndexResolution};

  SLIC_INFO("Building UniformGrid index...");
  UniformGrid3 ugrid(minBBPt.data(), maxBBPt.data(), resolutions);

  for (int i=0; i < ncells; i++) {
    SpatialBoundingBox triBB;
    t1=getMeshTriangle(i,  surface_mesh);
    triBB.addPoint(t1[0]);
    triBB.addPoint(t1[1]);
    triBB.addPoint(t1[2]);

    ugrid.insert(triBB, i);
  }


  // Iterate through triangle indices z from first index to last index.
  // Check against each other triangle with index greater than the index z
  // that also shares a UniformGrid bin.
  SLIC_INFO("Checking mesh with a total of " << ncells << " cells.");
  for (int z=0; z< ncells; z++) {
    seen.clear();
    hit.clear();

    // Retrieve the triangle at index z and construct a bounding box around it
    SpatialBoundingBox triBB2;
    t1 = getMeshTriangle(z,  surface_mesh);
 
    if (t1.degenerate()) { 
      degenerateIndices.push_back(z);
      continue;
    }
    triBB2.addPoint(t1[0]);
    triBB2.addPoint(t1[1]);
    triBB2.addPoint(t1[2]);

    Point3 minBBPt2,maxBBPt2;

    minBBPt2 = triBB2.getMin();
    maxBBPt2 = triBB2.getMax();

    // Get a list of all triangles in bins this triangle will touch
    std::vector<int> neighborTriangles;
    const std::vector<int> binsToCheck = ugrid.getBinsForBbox(triBB2);
    for (size_t curbin = 0; curbin < binsToCheck.size(); ++curbin) {
      std::vector<int> ntlist = ugrid.getBinContents(binsToCheck[curbin]);
      neighborTriangles.insert(neighborTriangles.end(), 
        ntlist.begin(), ntlist.end());
    }
    std::sort(neighborTriangles.begin(), neighborTriangles.end());
    std::vector<int>::iterator nend = 
      std::unique(neighborTriangles.begin(), neighborTriangles.end());
    std::vector<int>::iterator nit = neighborTriangles.begin();

    // remove triangles with indices less than or equal to this tri
    while (nit != nend && *nit <= z) {
      ++nit;
    }
    // test any remaining neighbor tris for intersection
    while (nit != nend) {
      t2 = getMeshTriangle(*nit, surface_mesh);
      if (checkTT(t1, t2)) {
        retval.push_back(std::make_pair(z, *nit));
      }
      ++nit;
    }
  }

  return retval;
}

} // end namespace detail
} // end namespace quest
} // end namespace axom

#endif  // MESH_TESTER_IMPL_HPP_

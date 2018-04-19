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

#include "axom/Types.hpp"

#include "primal/BoundingBox.hpp"
#include "primal/intersect.hpp"
#include "primal/Point.hpp"
#include "primal/Triangle.hpp"
#include "primal/UniformGrid.hpp"
#include "primal/RectangularLattice.hpp"
#include "primal/MortonIndex.hpp"

// C++ includes
#include <cmath>
#include <algorithm>
#include <vector>

#ifdef AXOM_USE_CXX11
#include <unordered_map>
#else
#include <map>
#endif

namespace axom
{
namespace quest
{

typedef mint::UnstructuredMesh< mint::Topology::SINGLE > UMesh;
typedef primal::Triangle<double, 3> Triangle3;

typedef primal::Point<double, 3> Point3;
typedef primal::BoundingBox<double, 3> SpatialBoundingBox;
typedef primal::UniformGrid<int, 3> UniformGrid3;
typedef primal::Vector<double, 3> Vector3;
typedef primal::Segment<double, 3> Segment3;

typedef axom::common::uint64 IndexType;
typedef primal::RectangularLattice<3, double, IndexType> Lattice3;
typedef primal::Mortonizer<IndexType, IndexType, 3> Morton3;

#ifdef AXOM_USE_CXX11
typedef std::unordered_map<IndexType, IndexType> MortonMap;
#else
typedef std::map<IndexType, IndexType> MortonMap;
#endif
typedef std::pair<MortonMap::iterator, bool> MortonMapResult;

inline SpatialBoundingBox compute_bounds( UMesh* mesh)
{
  SLIC_ASSERT( mesh != AXOM_NULLPTR );

  SpatialBoundingBox meshBB;
  Point3 pt;

  const double* x = mesh->getCoordinateArray( mint::X_COORDINATE );
  const double* y = mesh->getCoordinateArray( mint::Y_COORDINATE );
  const double* z = mesh->getCoordinateArray( mint::Z_COORDINATE );

  const mint::IndexType numNodes = mesh->getNumberOfNodes();
  for ( mint::IndexType i=0 ; i < numNodes ; ++i )
  {
    pt[ 0 ] = x[ i ];
    pt[ 1 ] = y[ i ];
    pt[ 2 ] = z[ i ];

    meshBB.addPoint( pt );
  } // END for all nodes

  if(numNodes > 0)
  {
    SLIC_ASSERT( meshBB.isValid() );
  }

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

inline Triangle3 getMeshTriangle(mint::IndexType i, UMesh* surface_mesh)
{
  SLIC_ASSERT( surface_mesh->getNumberOfCellNodes( i ) == 3);

  Triangle3 tri;

  const mint::IndexType* triCell = surface_mesh->getCell( i );

  const double* x = surface_mesh->getCoordinateArray( mint::X_COORDINATE );
  const double* y = surface_mesh->getCoordinateArray( mint::Y_COORDINATE );
  const double* z = surface_mesh->getCoordinateArray( mint::Z_COORDINATE );

  for ( int n=0; n < 3; ++n )
  {
    const mint::IndexType nodeIdx = triCell[ n ];
    tri[ n ][ 0 ] = x[ nodeIdx ];
    tri[ n ][ 1 ] = y[ nodeIdx ];
    tri[ n ][ 2 ] = z[ nodeIdx ];
  }

  return tri;
}

inline bool areTriangleIndicesDistinct( mint::IndexType* indices)
{
  SLIC_ASSERT(indices != AXOM_NULLPTR);

  return indices[0] != indices[1]
         && indices[1] != indices[2]
         && indices[2] != indices[0];
}

/*!
 * \brief Implementation for quest::findTriMeshIntersections
 *
 * \see findTriMeshIntersections
 */
void findTriMeshIntersections(
  UMesh* surface_mesh,
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

  const int ncells = surface_mesh->getNumberOfCells();

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


/* Weld vertices of a triangle mesh that are closer than \a eps  */
void weldTriMeshVertices(UMesh** surface_mesh,double eps)
{
  /// Implementation notes:
  ///
  /// This function welds vertices in the triangle mesh by
  /// quantizing them to an integer lattice with spacing \a eps.
  ///
  /// To avoid searching neighbor elements on the grid, we run the
  /// algorithm twice. Once on the original lattice, and once on a
  /// lattice that is shifted by half the grid spacing.
  ///
  /// Due to running this algorithm twice, it is possible in extreme cases
  /// for vertices that are at distances up to (1.5 * eps)
  /// (under the max norm) to be identified in this process.
  /// This can correspond to distances of at most 1.5 * eps * \sqrt(2)
  /// in the Euclidean norm.

  SLIC_ASSERT_MSG(
    eps > 0.,
    "Epsilon must be greater than 0. Passed in value was " << eps);
  SLIC_ASSERT_MSG(
    surface_mesh != AXOM_NULLPTR
    && *surface_mesh != AXOM_NULLPTR,
    "surface_mesh must be a valid pointer to a pointer to a triangle mesh");

  int const DIM = 3;
  UMesh* oldMesh = *surface_mesh;

  SpatialBoundingBox meshBB = compute_bounds(oldMesh).expand(eps);

  // Run the algorithm twice -- on the original grid and a translated grid
  std::vector<double> offsets;
  offsets.push_back(0);
  offsets.push_back(eps/2.);
  for(std::vector<double>::const_iterator it = offsets.begin() ;
      it != offsets.end() ; ++it)
  {
    // We will build up a new triangle mesh with the welded indices
    UMesh* newMesh = new UMesh(DIM, mint::TRIANGLE);

    // Set up the lattice for quantizing points to an integer lattice
    Point3 origin(meshBB.getMin().array() - Point3(*it).array() );
    Lattice3 lattice( origin, Point3(eps));

    // A map from Morton indices to the new vertex indices
    MortonMap vertexIndexMap;

    // First, find unique indices for the welded vertices
    const int numVerts = oldMesh->getNumberOfNodes();
    int uniqueVertCount = 0;

    std::vector<int> vertex_remap; // stores the new vertex indices
    vertex_remap.resize(numVerts); // for each old vertex

    Point3 vert;
    const double* x = oldMesh->getCoordinateArray( mint::X_COORDINATE );
    const double* y = oldMesh->getCoordinateArray( mint::Y_COORDINATE );
    const double* z = oldMesh->getCoordinateArray( mint::Z_COORDINATE );

    for(int i =0 ; i < numVerts ; ++i)
    {
      // get the vertex from the mesh
      vert[ 0 ] = x[ i ];
      vert[ 1 ] = y[ i ];
      vert[ 2 ] = z[ i ];

      // find the Morton index of the point w.r.t. the lattice
      IndexType morton = Morton3::mortonize(lattice.gridCell(vert));

      // find the new vertex index; if not present, insert vertex into new mesh
      MortonMapResult res = vertexIndexMap.insert(
        std::make_pair(morton, uniqueVertCount) );
      if(res.second == true)
      {
        uniqueVertCount++;
        newMesh->appendNodes(vert.data());
      }
      vertex_remap[i] = res.first->second;
    }

    // Next, add triangles into the new mesh using the unique vertex indices
    const int NUM_TRI_VERTS = 3;
    mint::IndexType triInds[ NUM_TRI_VERTS ];
    const mint::IndexType numTris = oldMesh->getNumberOfCells();
    for( mint::IndexType i =0 ; i < numTris ; ++i)
    {
      memcpy( triInds, oldMesh->getCell( i ),
              NUM_TRI_VERTS*sizeof( mint::IndexType ) );

      for(int d =0 ; d < NUM_TRI_VERTS ; ++d)
      {
        triInds[d] = vertex_remap [ triInds[d] ];
      }

      // Degeneracy check -- vertices need to be distinct
      if( areTriangleIndicesDistinct(triInds) )
      {
        newMesh->appendCell(triInds, mint::TRIANGLE);
      }
    }

    // Finally, delete old mesh and swap pointers
    delete oldMesh;
    oldMesh = newMesh;
  }

  // Update the original mesh pointer to the new mesh
  *surface_mesh = oldMesh;
}


} // end namespace quest
} // end namespace axom

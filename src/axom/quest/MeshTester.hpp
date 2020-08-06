// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MESH_TESTER_HPP_
#define MESH_TESTER_HPP_

// Axom includes
#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/primal.hpp"
#include "axom/spin.hpp"
#include "axom/mint.hpp"

// C++ includes
#include <cmath>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <functional> // for std::hash

// BVH includes
#if defined(AXOM_USE_RAJA)
  #include "RAJA/RAJA.hpp"
  #include "axom/spin/BVH.hpp"
  #include "axom/mint/execution/internal/structured_exec.hpp"
#endif

// C/C++ includes
#include <utility>
#include <vector>

/*!
 * \file MeshTester.hpp
 * \brief Defines functions to test Quest meshes for common defects.
 */

namespace axom
{
namespace quest
{

namespace detail
{
  using UMesh = mint::UnstructuredMesh< mint::SINGLE_SHAPE >;
  using Triangle3 = primal::Triangle<double, 3>;
  using SpatialBoundingBox = primal::BoundingBox<double, 3>;
  using UniformGrid3 = spin::UniformGrid<int, 3>;
  using Point3 = primal::Point<double, 3>;
}

/*! Enumeration indicating mesh watertightness */
enum class WatertightStatus : signed char
{
  WATERTIGHT = 0,    ///< Each edge in a surface mesh is incident in two cells
  NOT_WATERTIGHT,    ///< Each edge is incident in one or two cells
  CHECK_FAILED       ///< Calculation failed (possibly a non-manifold mesh)
};

inline detail::Triangle3 getMeshTriangle(axom::IndexType i, 
                                         detail::UMesh* surface_mesh)
{
  SLIC_ASSERT( surface_mesh->getNumberOfCellNodes( i ) == 3);

  detail::Triangle3 tri;

  const axom::IndexType* triCell = surface_mesh->getCellNodeIDs( i );

  const double* x = surface_mesh->getCoordinateArray( mint::X_COORDINATE );
  const double* y = surface_mesh->getCoordinateArray( mint::Y_COORDINATE );
  const double* z = surface_mesh->getCoordinateArray( mint::Z_COORDINATE );

  for ( int n=0 ; n < 3 ; ++n )
  {
    const axom::IndexType nodeIdx = triCell[ n ];
    tri[ n ][ 0 ] = x[ nodeIdx ];
    tri[ n ][ 1 ] = y[ nodeIdx ];
    tri[ n ][ 2 ] = z[ nodeIdx ];
  }

  return tri;
}

/// \name Mesh test and repair
/// @{

/*!
 * \brief Find self-intersections and degenerate triangles in a surface mesh
 *  utilizing a Bounding Volume Hierarchy.
 *
 * \param [in] surface_mesh A triangle mesh in three dimensions
 * \param [out] intersection Pairs of indices of intersecting mesh triangles
 * \param [out] degenerateIndices indices of degenerate mesh triangles
 * \param [in] intersectionThreshold Tolerance threshold for triangle 
 * intersection tests (default: 1E-8)
 * After running this function over a surface mesh, intersection will be filled
 * with pairs of indices of intersecting triangles and degenerateIndices will
 * be filled with the indices of the degenerate triangles in the mesh.
 * Triangles that share vertex pairs (adjacent triangles in a watertight
 * surface mesh) are not reported as intersecting.  Degenerate triangles
 * are not reported as intersecting other triangles.
 *
 */
#if defined(AXOM_USE_RAJA)
template < typename ExecSpace, typename FloatType >
void findTriMeshIntersectionsBVH(
  mint::UnstructuredMesh< mint::SINGLE_SHAPE >* surface_mesh,
  std::vector<std::pair<int, int> > & intersections,
  std::vector<int> & degenerateIndices,
  double intersectionThreshold = 1E-8)
{
  SLIC_INFO("Running BVH intersection algorithm "
            << " in execution Space: "
            << axom::execution_space< ExecSpace >::name());

  // Get allocator
  int allocatorID = axom::execution_space< ExecSpace >::allocatorID();
  axom::setDefaultAllocator( allocatorID );

  constexpr int NDIMS = 3;
  constexpr int stride = 2 * NDIMS;
  const int ncells = surface_mesh->getNumberOfCells();

  detail::Triangle3* tris = axom::allocate <detail::Triangle3> (ncells);

  FloatType* xmin = axom::allocate< FloatType >( ncells );
  FloatType* ymin = axom::allocate< FloatType >( ncells );
  FloatType* zmin = axom::allocate< FloatType >( ncells );

  FloatType* xmax = axom::allocate< FloatType >( ncells );
  FloatType* ymax = axom::allocate< FloatType >( ncells );
  FloatType* zmax = axom::allocate< FloatType >( ncells );

  // Marks each cell/triangle as degenerate (1) or not (0)
  int* degenerate = axom::allocate< int >( ncells );

  // Each access-aligned bounding box represented by 2 (x,y,z) points
  FloatType* aabbs = axom::allocate< FloatType >( ncells * stride );

  // Finding degenerate indices and bounding boxes from Triangles
  detail::Triangle3 t1 = detail::Triangle3();
  detail::Triangle3 t2 = detail::Triangle3();

  // Initialize the bounding box for each Triangle and marks
  // if the Triangle is degenerate.
  mint::for_all_cells< ExecSpace, mint::xargs::coords >(
    surface_mesh, AXOM_LAMBDA( IndexType cellIdx,
                               numerics::Matrix< double >&coords,
                               const IndexType* AXOM_NOT_USED(nodeIds) )
  {
    detail::Triangle3 tri;

    for ( IndexType inode=0 ; inode < 3 ; ++inode )
    {
      const double* node = coords.getColumn( inode );
      tri[ inode ][ 0 ] = node[ mint::X_COORDINATE ];
      tri[ inode ][ 1 ] = node[ mint::Y_COORDINATE ];
      tri[ inode ][ 2 ] = node[ mint::Z_COORDINATE ];
    } // END for all cells nodes

    if (tri.degenerate())
    {
      degenerate[cellIdx] = 1;
    }
    else
    {
      degenerate[cellIdx] = 0;
    }

    tris[cellIdx] = tri;

    detail::SpatialBoundingBox triBB = compute_bounding_box(tri);

    const IndexType offset = cellIdx * stride;

    xmin[cellIdx] = aabbs[ offset ]     = triBB.getMin()[0];
    ymin[cellIdx] = aabbs[ offset + 1 ] = triBB.getMin()[1];
    zmin[cellIdx] = aabbs[ offset + 2 ] = triBB.getMin()[2];

    xmax[cellIdx] = aabbs[ offset + 3 ] = triBB.getMax()[0];
    ymax[cellIdx] = aabbs[ offset + 4 ] = triBB.getMax()[1];
    zmax[cellIdx] = aabbs[ offset + 5 ] = triBB.getMax()[2];
  } );

  // Return degenerateIndices
  for (int i = 0 ; i < ncells; i++)
  {
    if (degenerate[i] == 1)
    {
      degenerateIndices.push_back(i);
    }
  }

  // Construct BVH
  axom::spin::BVH< NDIMS, ExecSpace, FloatType > bvh( aabbs, ncells );
  bvh.build( );

  // Run find algorithm
  IndexType* offsets    = axom::allocate< IndexType >( ncells );
  IndexType* counts     = axom::allocate< IndexType >( ncells );
  IndexType* candidates = nullptr;
  bvh.findBoundingBoxes( offsets, counts, candidates, ncells, xmin, xmax,
                         ymin, ymax, zmin, zmax );

  // Get the total number of candidates
  int totalCandidates = 0;
  for (int i = 0 ; i < ncells ; i++ )
  {
    totalCandidates += counts[i];
  }

  //Deallocate no longer needed variables
  axom::deallocate(aabbs);
  
  axom::deallocate(degenerate);

  axom::deallocate(xmin);
  axom::deallocate(ymin);
  axom::deallocate(zmin);

  axom::deallocate(xmax);
  axom::deallocate(ymax);
  axom::deallocate(zmax);

  int* intersection_pairs =
    axom::allocate <int> (totalCandidates * 2);

  using ATOMIC_POL =
          typename axom::execution_space< ExecSpace >::atomic_policy;
  int* counter = axom::allocate <int> (1);
  counter[0] = 0;

  for_all< ExecSpace >( ncells, AXOM_LAMBDA (IndexType i)
  {
    for (int j = 0 ; j < counts[i] ; j++)
    {
      int candidate_index = candidates[ offsets[i] + j];
      if (i != candidate_index && i < candidate_index)
      {
        if (primal::intersect(tris[i], tris[candidate_index],
                              false, intersectionThreshold))
        {
          auto idx = RAJA::atomicAdd<ATOMIC_POL>(counter, 2);
          intersection_pairs[idx] = i;
          intersection_pairs[idx + 1] = candidate_index;
        }
      }
    }
  } );

  // Initialize pairs of clashes
  for (int i = 0 ; i < counter[0] ; i += 2)
  {
    intersections.push_back(std::make_pair(intersection_pairs[i],
                                           intersection_pairs[i + 1]));
  }

  // Deallocate
  axom::deallocate(tris);

  axom::deallocate(offsets);
  axom::deallocate(counts);
  axom::deallocate(candidates);

  axom::deallocate(intersection_pairs);
  axom::deallocate(counter);
}
#endif

/*!
 * \brief Find self-intersections and degenerate triangles in a surface mesh
 *  utilizing a Uniform Grid.
 *
 * \param [in] surface_mesh A triangle mesh in three dimensions
 * \param [out] intersection Pairs of indices of intersecting mesh triangles
 * \param [out] degenerateIndices indices of degenerate mesh triangles
 * \param [in] spatialIndexResolution The grid resolution for the index
 * structure (default: 0)
 * \param [in] intersectionThreshold Tolerance threshold for triangle 
 * intersection tests (default: 1E-8)
 *
 * After running this function over a surface mesh, intersection will be filled
 * with pairs of indices of intersecting triangles and degenerateIndices will
 * be filled with the indices of the degenerate triangles in the mesh.
 * Triangles that share vertex pairs (adjacent triangles in a watertight
 * surface mesh) are not reported as intersecting.  Degenerate triangles
 * are not reported as intersecting other triangles.
 *
 * This function uses a quest::UniformGrid spatial index.  Input
 * spatialIndexResolution specifies the bin size for the UniformGrid.  The
 * default value of 0 causes this routine to calculate a heuristic bin size
 * based on the cube root of the number of cells in the mesh.
 */
void findTriMeshIntersections(
  mint::UnstructuredMesh< mint::SINGLE_SHAPE >* surface_mesh,
  std::vector<std::pair<int, int> > & intersections,
  std::vector<int> & degenerateIndices,
  int spatialIndexResolution = 0,
  double intersectionThreshold = 1E-8);


/*!
 * \brief Check a surface mesh for holes using its face relation.
 *
 * \param [in] surface_mesh A surface mesh in three dimensions
 *
 * \returns status If the mesh is watertight, is not watertight, or
 *    if an error occurred (possibly due to non-manifold mesh).
 *
 * \note This method marks the cells on the boundary by creating a new
 *  cell-centered field variable, called "boundary", on the given input mesh.
 *
 * \note This function computes the mesh's cell-face and face-vertex relations.
 * For large meshes, this can take a long time.  The relations are used to
 * check for holes, and remain cached with the mesh after this function
 * finishes.
 */
WatertightStatus isSurfaceMeshWatertight(
  mint::UnstructuredMesh< mint::SINGLE_SHAPE >* surface_mesh );


/*!
 * \brief Mesh repair function to weld vertices that are closer than \a eps
 *
 * \param [in,out] surface_mesh A pointer to a pointer to a triangle mesh
 * \param [in] eps Distance threshold for welding vertices (using the max norm)
 *
 * \pre \a eps must be greater than zero
 * \pre \a surface_mesh is a pointer to a pointer to a non-null triangle mesh.
 * \post The triangles of \a surface_mesh are reindexed using the welded
 * vertices and degenerate triangles are removed.  The mesh can still contain
 * vertices that are not referenced by any triangles.
 *
 * This utility function repairs an input triangle mesh (embedded in three
 * dimensional space) by 'welding' vertices that are closer than \a eps.
 * The vertices are quantized to an integer lattice with spacing \a eps
 * and vertices that fall into the same cell on this lattice are identified.
 * All identified vertices are given the coordinates of the first such vertex
 * and all incident triangles use the same index for this vertex.
 *
 * The input mesh can be a "soup of triangles", where the vertices
 * of adjacent triangles have distinct indices.  After running this function,
 * vertices that are closer than \a eps are welded, and their incident
 * triangles use the new vertex indices.  Thus, the output is an
 * "indexed triangle mesh".
 *
 * This function also removes degenerate triangles from the mesh.  These
 * are triangles without three distinct vertices after the welding.
 *
 * \note This function is destructive.  It modifies the input triangle
 * mesh in place.
 * \note The distance metric in this function uses the "max" norm (l_inf).
 */
void weldTriMeshVertices(
  mint::UnstructuredMesh< mint::SINGLE_SHAPE >** surface_mesh,
  double eps);

/// @}

} // end namespace quest
} // end namespace axom

#endif   // MESH_TESTER_HPP_

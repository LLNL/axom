// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MESH_TESTER_HPP_
#define MESH_TESTER_HPP_

// Axom includes
#include "axom/config.hpp"
#include "axom/mint/mesh/UnstructuredMesh.hpp"
#include "axom/slic/interface/slic.hpp"

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

/*! Enumeration indicating mesh watertightness */
enum class WatertightStatus : signed char
{
  WATERTIGHT = 0,    ///< Each edge in a surface mesh is incident in two cells
  NOT_WATERTIGHT,    ///< Each edge is incident in one or two cells
  CHECK_FAILED       ///< Calculation failed (possibly a non-manifold mesh)
};


/// \name Mesh test and repair
/// @{

/*!
 * \brief Find self-intersections and degenerate triangles in a surface mesh.
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

/**
 * \file
 * \brief Defines functions to test Quest meshes for common defects.
 */

#ifndef MESH_TESTER_HPP_
#define MESH_TESTER_HPP_

// Axom includes
#include "mint/Mesh.hpp"

// C/C++ includes
#include <vector>
#include <utility>

#include "slic/slic.hpp"

#include "quest/MeshTester_impl.hpp"

namespace axom {
namespace quest {

/**
 * \brief For a given Mesh, find self-intersections and degenerate triangles.
 *
 * Input surface_mesh should point to a mint::Mesh in three dimensions
 * containing 2D triangles that form a surface.  After running this function,
 * degenerateIndices will be filled with the mesh cell indices of degenerate
 * triangles.  The return value lists triangle index pairs that intersect.
 * Triangles that share vertex pairs (adjacent triangles in a watertight
 * surface mesh) are not reported as intersecting.  Degenerate triangles
 * are not reported as intersecting other triangles.
 *
 * This function uses a quest::UniformGrid spatial index.  Input
 * spatialIndexResolution specifies the bin size for the UniformGrid.  The
 * default value of 0 causes this routine to calculate a heuristic bin size
 * based on the cube root of the number of cells in the mesh.
 */
std::vector< std::pair<int, int> >
findTriMeshIntersections(mint::Mesh* surface_mesh,
			 std::vector<int> & degenerateIndices,
			 int spatialIndexResolution = 0)
{
  return detail::findTriMeshIntersections_impl(surface_mesh, 
                                               degenerateIndices, 
                                               spatialIndexResolution);
}

} // end namespace quest
} // end namespace axom

#endif   // MESH_TESTER_HPP_

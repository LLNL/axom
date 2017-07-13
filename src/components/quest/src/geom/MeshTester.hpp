

#ifndef MESH_TESTER_HPP_
#define MESH_TESTER_HPP_

// Axom includes
// #include "mint/UniformMesh.hpp"
#include "mint/Mesh.hpp"

// C/C++ includes
#include <vector>
#include <utility>

#include "slic/slic.hpp"

#include "quest/MeshTester_impl.hpp"

namespace axom {
namespace quest {

std::vector< std::pair<int, int> > findTriMeshIntersections(mint::Mesh* surface_mesh,
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

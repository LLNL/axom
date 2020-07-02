Quest {#questtop}
========

[Quest](@ref axom::quest) provides spatial queries and operations on a mesh.
- Classes to [read in](@ref axom::quest::STLReader) and functions to
  [test and repair](@ref findTriMeshIntersections()) surface meshes
  composed of 2D triangles embedded in 3D
- [Functions](@ref signed_distance_init()) to query the signed distance of a
  point to a 3D surface mesh
- [Functions](@ref inout_init()) to query whether a point is inside a 3D surface mesh
- A [class](@ref axom::quest::PointInCell) to locate a physical point in a cell
  of a high-order mesh
- A [function](@ref all_nearest_neighbors()) to find the nearest neighbor for each
  of a set of points

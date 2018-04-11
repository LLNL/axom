QUEST: Query points in meshes {#questtop}
========

[QUEST](@ref axom::quest) provides queries on meshes.
- distance() and inside() query a mesh relative to a point
- mesh_min_bounds(), mesh_max_bounds(), mesh_center_of_mass(), and compute_bounds() report characteristics of a mesh
- The [SignedDistance](@ref axom::quest::SignedDistance) and [InOutOctree](@ref axom::quest::InOutOctree) classes implement the query functionality and are available for code use independent of the function interface

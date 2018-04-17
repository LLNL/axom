QUEST: Query points in meshes {#questtop}
========

[QUEST](@ref axom::quest) provides queries on triangle surface meshes embedded in 3D.
- Functions to query whether a point is inside or outside a mesh, or the point's distance to the mesh
- Functions to retrieve a mesh's bounds and center of mass
- The [SignedDistance](@ref axom::quest::SignedDistance) and [InOutOctree](@ref axom::quest::InOutOctree) classes implement the query functionality and are available for code use independent of the function interface

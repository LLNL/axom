Spin: Spatial indexes, spatial acceleration data structures {#spintop}
========

[Spin](@ref axom::spin) provides data structures to index space and their
associated support classes.
- Tree structures
  - [BVHTree](@ref axom::spin::BVHTree) creates a bounding volume hierarchy
  - [SpatialOctree](@ref axom::spin::SpatialOctree) creates an octree
- One-level structures
  - [UniformGrid](@ref axom::spin::UniformGrid) divides a region of interest
    into equal-sized bins
  - [Mortonizer](@ref axom::spin::Mortonizer) assigns each point in a region
    of interest to a point on a 1D space-filling curve
  - [ImplicitGrid](@ref axom::spin::ImplicitGrid) creates a set of bins for
    each axis, assigning every object to the range of bins it falls into
  - [RectangularLattice](@ref axom::spin::RectangularLattice) supports
    computation of coordinates for spatial bins.

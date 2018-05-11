******************************************************
Introductory examples
******************************************************

Here is a collection of introductory examples showing Primal primitives and
operations.  

- instantiate a point, a ray, several triangles, a bounding box, a plane
- get triangles' area and normal; test for degeneracy
- make a visualization (Geomview?)

With these primitives, we can perform some geometric operations.

- clip triangle against bbox to get polygon
- find which side of the plane a triangle point lies on
- find intersection point btw ray and triangle
- find the closest point from another triangle corner to another tri, and
  the squared distance btw points

A spatial index can be used when a code compares each primitive in a collection
to every other primitive, such as when checking if a triangle mesh intersects
itself.  The following naive implementation is straightforward but runs
in $O(n^2)$ time, where $n$ is the number of triangles.

- naive implementation: compare each triangle to every other tri

A code should be able to skip the call to ``intersect`` for widely-separated
primitives.  The UniformGrid is designed for this optimization.

- UniformGrid example.  For each triangle,
  - find its bounding box
  - get the bins intersecting the bounding box
  - iterate over the contents of each bin, testing for intersection

Find examples for the MortonIndex and the BVH tree.

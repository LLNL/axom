BVHTree
^^^^^^^

The ``BVHTree`` recursively subdivides the region of interest into a "tree" of
subregions, stopping when a subregion contains less than some number of objects
or when the tree reaches a specified height.  Similar to ``UniformGrid``,
subregions are also called "buckets".

The ``BVHTree`` is well-suited for particle-mesh or ray-mesh intersection tests.
It is also well-suited to data sets where the contents are unevenly distributed,
since the buckets are subdivided based on their contents.  The figure below
shows a 2D BVH tree built up over a set of triangles.

.. figure:: figs/showUniformGrid.png
   :figwidth: 300px
   :alt: Diagram showing triangles indexed with a UniformGrid

The following code example shows how a ``BVHTree`` can be used to accelerate 
a ray-mesh intersection algorithm.  Without the spatial index, each ray must
be tested against each triangle.  The ``BVHTree`` prunes out all buckets that 
do not intersect the probe, recursing into the tree of buckets until it gets
to the leaf buckets.  Then all objects in the leaf buckets are compared to the
probe.

.. literalinclude:: ../../examples/primal_introduction.cpp
   :start-after: _ugrid_triintersect_header_start
   :end-before: _ugrid_triintersect_header_end
   :language: C++

.. literalinclude:: ../../examples/primal_introduction.cpp
   :start-after: _ugrid_triintersect_start
   :end-before: _ugrid_triintersect_end
   :language: C++

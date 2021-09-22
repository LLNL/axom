BVH
^^^

The ``BVH`` class implements a
`bounding volume hierarchy <https://en.wikipedia.org/wiki/Bounding_volume_hierarchy>`_.
This data structure recursively subdivides a rectilinear region of interest
into a "tree" of subregions.

``BVH`` is well-suited for particle-mesh or ray-mesh intersection tests.
It is also well-suited to data sets where the contents are unevenly distributed,
since the bins are subdivided based on their contents.

The following code example shows how a ``BVH`` can be used to accelerate a
point-mesh intersection algorithm.  First, we generate bounding boxes for all
triangles in the mesh, and call ``BVH::initialize()`` with the bounding boxes.

.. literalinclude:: ../../examples/spin_introduction.cpp
   :start-after: _bvh_header_start
   :end-before: _bvh_header_end
   :language: C++

.. literalinclude:: ../../examples/spin_introduction.cpp
   :start-after: _bvh_build_start
   :end-before: _bvh_build_end
   :language: C++

After the structure is built, we can use the BVH to generate a list of element IDs
that are candidate neighbors to the query points.  Call ``BVH::findPoints()``
to get the list of element IDs that the query point intersects.
The key idea of ``findPoints()`` is that testing for probe intersection with a
child node (bounding box) is cheap.  If a node intersection test fails (misses),
the child node, and all of its corresponding elements, can be skipped during the
BVH traversal.  If the probe does intersect a child node, the node's children are also
tested for probe intersection.  Without the acceleration data structure, each
probe point must be tested against each triangle.

.. literalinclude:: ../../examples/spin_introduction.cpp
   :start-after: _bvh_candidate_start
   :end-before: _bvh_candidate_end
   :language: C++

Finally, test the point against all candidate neighbor triangles.

.. literalinclude:: ../../examples/spin_introduction.cpp
   :start-after: _bvh_cand_int_start
   :end-before: _bvh_cand_int_end
   :language: C++

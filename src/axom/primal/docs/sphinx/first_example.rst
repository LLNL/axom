******************************************************
Introductory examples
******************************************************

Here is a collection of introductory examples showing Primal primitives and
operations.  We will instantiate several geometric primitives as needed, perform
geometric operations, and build the UniformGrid and BVHTree spatial index
objects around them.  These examples show representative overloads of each of
the Primal operations (see the 
`API documentation <../../../doxygen/axom_doxygen/html/primaltop.html>`_ 
for more details).

Include header files for primitives (header files for operations will be shown
next to code examples).

.. literalinclude:: ../../examples/primal_introduction.cpp
   :start-after: _prims_header_start
   :end-before: _prims_header_end
   :language: C++

.. literalinclude:: ../../examples/primal_introduction.cpp
   :start-after: _using_start
   :end-before: _using_end
   :language: C++

Operations and Primitives
-------------------------

Clip triangle against bounding box
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The clip operation clips a triangle against a bounding box, returning the
resulting polygon.  The figure shows the triangle in blue and the polygon
resulting from ``clip()`` in yellow.

.. figure:: figs/showClip.png
   :figwidth: 300px
   :alt: A polygon is produced by clipping a triangle.

.. literalinclude:: ../../examples/primal_introduction.cpp
   :start-after: _clip_header_start
   :end-before: _clip_header_end
   :language: C++

.. literalinclude:: ../../examples/primal_introduction.cpp
   :start-after: _clip_start
   :end-before: _clip_end
   :language: C++

Closest point query
^^^^^^^^^^^^^^^^^^^

The closest point operation finds the point on a triangle that is closest to a
query point.  Query point :math:`o` (shown in dark blue), at the origin, is
closest to point :math:`o'` (light blue), which lies in the triangle's interior.
Query point :math:`a` (olive) is closest to point :math:`a'` (yellow), which
lies on the triangle's edge.

.. figure:: figs/showClosestPoint.png
   :figwidth: 300px
   :alt: Diagram showing the closest point query.

.. literalinclude:: ../../examples/primal_introduction.cpp
   :start-after: _closest_point_header_start
   :end-before: _closest_point_header_end
   :language: C++

.. literalinclude:: ../../examples/primal_introduction.cpp
   :start-after: _closest_point_start
   :end-before: _closest_point_end
   :language: C++

As the code example shows, ``closest_point()`` can take a pointer to an ``int``
as an optional third parameter.  If supplied, the function writes a value into
the ``int`` that indicates which of the triangle's vertices or sides contains
the closest point (or the interior).

Compute bounding box
^^^^^^^^^^^^^^^^^^^^

Primal's bounding boxes are rectangular right prisms.  That is, they are boxes
where neighboring walls are at right angles.

The BoundingBox class represents an axis-aligned bounding box, which has two
walls perpendicular to the X-axis, two perpendicular to the Y-axis, and two
perpendicular to the Z-axis.  This is sufficient for many computations; range
and intersection operations tend to be fast.

The OrientedBoundingBox class can be oriented in any way with respect to the 
coordinate axes.  This can provide a tighter fit to the bounded data, but 
construction, intersection, and range calculation are more costly.

Here a group of points is used to create both an (axis-aligned) BoundingBox
and an OrientedBoundingBox.  The points are drawn in blue, the BoundingBox in
black, and the OrientedBoundingBox in orange.

.. figure:: figs/showBoundingBoxes.png
   :figwidth: 300px
   :alt: Diagram showing (axis-aligned) BoundingBox and OrientedBoundingBox objects bounding the same set of points.

.. literalinclude:: ../../examples/primal_introduction.cpp
   :start-after: _bbox_header_start
   :end-before: _bbox_header_end
   :language: C++

.. literalinclude:: ../../examples/primal_introduction.cpp
   :start-after: _bbox_start
   :end-before: _bbox_end
   :language: C++

Primal also provides a ``merge_boxes()`` function to produce a bounding box that
contains two input bounding boxes.  This is available for client codes to use and
also supports the operation of the BVHTree spatial index.

Intersection
^^^^^^^^^^^^

The intersection test is provided by ``intersect()``.  It takes two primitives
and returns a boolean if the primitives intersect.  Some overloads return the
point of intersection in an output argument.  The overloads for
``intersect()`` are summarized in the table below.

+-------------------+-------------------+--------------------------------------+
|Arg 1              |Arg 2              |Additional arguments and notes        |
+===================+===================+======================================+
|Triangle           |Triangle           |include boundaries [#f1]_             |
|                   |                   |(default false)                       |
+-------------------+-------------------+--------------------------------------+
|Ray                |Segment            |return intersection point. 2D only.   |
+-------------------+-------------------+--------------------------------------+
|Segment            |BoundingBox        |return intersection point             |
+-------------------+-------------------+--------------------------------------+
|Ray                |BoundingBox        |return intersection point             |
+-------------------+-------------------+--------------------------------------+
|BoundingBox        |BoundingBox        |                                      |
+-------------------+-------------------+--------------------------------------+
|Sphere             |Sphere             |specify tolerance                     |
+-------------------+-------------------+--------------------------------------+
|Triangle           |BoundingBox        |                                      |
+-------------------+-------------------+--------------------------------------+
|Triangle           |Ray                |return parameterized intersection     |
|                   |                   |point (on Ray), return barycentric    |
|                   |                   |intersection point (on Triangle)      |
+-------------------+-------------------+--------------------------------------+
|Triangle           |Segment            |return parameterized intersection     |
|                   |                   |point (on Segment), return barycentric|
|                   |                   |intersection point (on Triangle)      |
+-------------------+-------------------+--------------------------------------+
|OrientedBoundingBox|OrientedBoundingBox|specify tolerance                     |
+-------------------+-------------------+--------------------------------------+

.. [#f1] By default, the triangle intersection algorithm considers only the
         triangles' interiors, so that non-coplanar triangles that share two
         vertices are not reported as intersecting.  The caller to
         ``intersect()`` can specify an optional argument to include triangle
         boundaries in the intersection test.

The example below intersection tests between two triangles, a ray, and a
BoundingBox.

.. figure:: figs/showIntersect.png
   :figwidth: 300px
   :alt: Diagram showing intersection tests.

.. literalinclude:: ../../examples/primal_introduction.cpp
   :start-after: _intersect_header_start
   :end-before: _intersect_header_end
   :language: C++

.. literalinclude:: ../../examples/primal_introduction.cpp
   :start-after: _intersect_start
   :end-before: _intersect_end
   :language: C++

If a ray intersects a bounding box at more than one point, as in the example above,
the first point along the ray is reported as the intersection.  Thus, the ray's
entry point is shown in the diagram, rather than some point in the box's interior
or where it exits the box.

Orientation
^^^^^^^^^^^

Axom contains two overloads of ``orientation()``.  The 3D case tests a point against
a triangle and reports which side it lies on; the 2D case tests a point against a
line segment.  Here is an example use of the 3D point-triangle orientation test.

.. figure:: figs/showOrientation.png
   :figwidth: 300px
   :alt: Diagram showing 3D point-triangle orientation test.

.. literalinclude:: ../../examples/primal_introduction.cpp
   :start-after: _orient_header_start
   :end-before: _orient_header_end
   :language: C++

.. literalinclude:: ../../examples/primal_introduction.cpp
   :start-after: _orient_start
   :end-before: _orient_end
   :language: C++

The triangle is shown with its normal vector pointing out of its centroid.  The
test point on the :math:`z` axis is on the negative side of the triangle.  The
centroid lies in the triangle---"on the boundary."  The remaining test point is
on the side of the triangle pointed into by the normal vector, the positive
side.

Distance
^^^^^^^^

The various overloads of ``squared_distance()`` calculate the squared 
distance between a query point and several different primitives:

- another point,
- a BoundingBox,
- a Segment,
- a Triangle.

.. figure:: figs/showDistance.png
   :figwidth: 300px
   :alt: Diagram showing 3D point-triangle orientation test.

.. literalinclude:: ../../examples/primal_introduction.cpp
   :start-after: _sqdist_header_start
   :end-before: _sqdist_header_end
   :language: C++

.. literalinclude:: ../../examples/primal_introduction.cpp
   :start-after: _sqdist_start
   :end-before: _sqdist_end
   :language: C++

The example shows the squared distance between the query point, shown in black
in the figure, and four geometric primitives.  For clarity, the diagram also shows
the projection of the query point and its closest points to the XY plane.

Spatial Index
-------------

A spatial index is a data structure used to speed up retrieval of geometric
objects.  Primal provides two spatial indices, the cell (or Verlet) list
implemented in the ``UniformGrid`` class and the bounding volume hierarchy tree
implemented in the ``BVHTree`` class.  Both classes divide a bounding box
denoting a region of interest into bins that group objects together, avoiding
the need to process objects that do not fall into a bin of interest.

UniformGrid
^^^^^^^^^^^

The ``UniformGrid`` tiles the region of interest into non-intersecting
subregions (or "buckets") of uniform size.  Each object in the region of interest
is added an object to every bucket the object
intersects.  ``UniformGrid`` can be used when a code compares each primitive
in a collection to every other "close" primitive, such as when checking if a
triangle mesh intersects itself.  The following naive implementation is
straightforward but runs in :math:`O(n^2)` time, where :math:`n` is the number
of triangles.

.. literalinclude:: ../../examples/primal_introduction.cpp
   :start-after: _naive_triintersect_start
   :end-before: _naive_triintersect_end
   :language: C++

We want to call ``intersect()`` only for triangles that can intersect, ignoring
widely-separated triangles.  The ``UniformGrid`` enables this optimization.  In
the following figure, the ``UniformGrid`` divides the region of interest into
three by three buckets outlined in grey.  A triangle :math:`t` (shown in orange)
will be compared with neighbor triangles (shown in black) that fall into the
buckets occupied by :math:`t`.  Other triangles (shown in blue) are too far away
to intersect and are not compared with :math:`t`.

.. Unlike the figures for geometric operations, this figure is "hand-drawn," not
   generated by running the example below.

.. figure:: figs/showUniformGrid.png
   :figwidth: 300px
   :alt: Diagram showing triangles indexed with a UniformGrid

.. literalinclude:: ../../examples/primal_introduction.cpp
   :start-after: _ugrid_triintersect_header_start
   :end-before: _ugrid_triintersect_header_end
   :language: C++

.. literalinclude:: ../../examples/primal_introduction.cpp
   :start-after: _ugrid_triintersect_start
   :end-before: _ugrid_triintersect_end
   :language: C++

The ``UniformGrid`` has its best effect when objects are roughly the same size and 
evenly distributed over the region of interest, and when buckets are close to the 
characteristic size of objects in the region of interest.  

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

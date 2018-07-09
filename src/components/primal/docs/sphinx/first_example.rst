******************************************************
Introductory examples
******************************************************

Here is a collection of introductory examples showing Primal primitives and
operations.  We will instantiate several geometric primitives as needed, perform
geometric operations, and build the UniformGrid and BSPTree spatial index
objects around them.  These examples show representative overloads of each of
the Primal operations (see the `API documentation <../../foo.html>`_ for more
details).

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
resulting polygon.

.. instantiate a point, a ray, several triangles, a bounding box, a plane
   get triangles' area and normal; test for degeneracy
   make a visualization (Geomview?)

.. figure:: figs/primal_clip.png
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

The closest point operation finds the point on a triangle that is closest to 
a query point.

.. figure:: figs/primal_closest_point.png
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

As the code example shows, code calling ``closest_point()`` can supply a
pointer to an int as an optional third parameter.  If supplied, the function
writes a value into the int that indicates which of the triangle's vertices
or sides contains the closest point (or the interior).

Compute bounding box
^^^^^^^^^^^^^^^^^^^^

Primal's bounding boxes are right rectangular prisms.  That is, they are boxes
where the walls are all at right angles.  

The BoundingBox class represents an axis-aligned bounding box, which has two
walls perpendicular to the X-axis, two perpendicular to the Y-axis, and two
perpendicular to the Z-axis.  This is sufficient for many computations and range
and intersection operations tend to be fast.

The OrientedBoundingBox class can be oriented in any way with respect to the 
coordinate axes.  This can provide a tighter fit to the bounded data, but 
construction, intersection, and range calculation can be more costly.

Here a group of points is used to create both an (axis-aligned) BoundingBox
and an OrientedBoundingBox.

.. figure:: figs/primal_bbox.png
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

Orientation
^^^^^^^^^^^

Distance
^^^^^^^^




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

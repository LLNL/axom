Intersection
^^^^^^^^^^^^

The intersection test is provided by ``intersect()``.  It takes two primitives
and returns a boolean indicating if the primitives intersect.  Some overloads
return the point of intersection in an output argument.  The overloads for
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

The example below tests for intersection between two triangles, a ray, and a
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
the first intersection point along the ray (the intersection closest to the ray's
origin) is reported as the intersection.  Thus, the ray's entry point is shown in
the diagram, rather than some point in the box's interior or where it exits the box.
If a ray originates inside a bounding box, the origin will be reported as the point
of intersection.

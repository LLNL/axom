
Primal User Documentation
=========================

Primal is a component of Axom that provides efficient implementations of 
fundamental geometric ideas.  Primal provides:

- Classes to represent geometric primitives such as Point and Ray
- Functions to implement operations on Primal's classes including distance and
  intersection
- Classes implementing two spatial indexes, used to accelerate operations
  such as collision detection and containment.


Primitives
----------

Primal includes the following primitives:

- Point
- Segment, Ray, Vector
- Plane, Triangle, Polygon
- Sphere
- Tetrahedron
- BoundingBox, OrientedBoundingBox

Primal also provides the NumericAray class, which implements arithmetic
operations on numerical tuples and supports Primal's Point and Vector classes.
Classes in Primal are templated on coordinate type (double, float, etc.) and
dimension.  The primitives do not inherit from a common base class.  This was a
design choice in favor of simplicity and performance.  Geometric primitives can
be tested for equality and can be printed to strings.

Primal also includes functions to merge two BoundingBox or two OrientedBoundingBox
objects and to create new OrientedBoundingBox objects from a list of points.


Operations
----------

Primal implements geometric operations with functions.  Currently, these include:

- ``closest_point`` takes a primitive P and a query point Q and returns the point
  on P closest to Q.
- ``squared_distance`` computes the squared distance from a point to another point,
  bounding box, line segment, or triangle.
- ``orientation`` finds the side of a line segment or triangle where a query point lies.
- ``clip`` finds the polygon resulting from a bounding box clipping a triangle.
- ``compute_bounding_box`` finds the bounding box for a given primitive.
- ``intersect`` predicate tests if two primitives intersect.  Some of the combinations
  also indicate the point of intersection of a 1D primitive with another primitive.

.. note::
   Most use cases have low dimension, usually 2 or 3.  Dimensionality has been
   generalized to support other values where it does not interfere with the
   common case, but some operations such as triangle intersection do not support
   other dimensionality than 2 or 3.

Spatial Indexes
---------------

There are two spatial indexes provided by Primal, represented by the BVHTree and
UniformGrid classes.  These are supported by the RectangularLattice class and
construction functions, and several classes implementing a Morton index.

Cell list
^^^^^^^^^

The `Cell list <https://en.wikipedia.org/wiki/Cell_lists>`_ is implemented by
the UniformGrid class.  A code using the UniformGrid adds objects with their
bounding boxes to the index, then calls ``build()``.  UniformGrid divides the region
of interest into non-overlapping sub-regions ("buckets"), all of equal size.
Each bucket lists all the objects whose bounding box overlaps the bucket.

The cell list is used where each object will interact with its neighbors.
Bucket size should be chosen to be equal to or larger than the range of
interaction.  This way, when calculating interactions with a particular object,
an application can rule out all objects that do not fall into the object's
bucket or its immediate neighbors.

Bounding volume hierarchy tree
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The `bounding volume hierarchy tree <https://en.wikipedia.org/wiki/Bounding_volume_hierarchy>`_
is implemented by the BVHTree class.  A code will add objects with associated bounding
boxes to the BVHTree, then call the ``build()`` method.  This method hierarchically
subdivides a region of space until the smallest subdivisions (leaf buckets)
contain no more than a certain number of objects (default 25).  Parent buckets
entirely contain their child buckets, and child buckets may overlap.

This spatial index is often used in ray casting or collision detection, where a
probe object such as a particle or ray can potentially hit anything in the
region of interest.  The main idea is that testing for probe intersection
with a bucket (bounding box) is cheap.  If a bucket intersection test fails (misses),
the contents of the bucket are cheaply pruned out of the search.  If the probe
does intersect a bucket, the next level of buckets is tested for probe intersection.
This process continues until the lowest level of BVH tree buckets, where individual
objects are tested for probe intersection.

Space filling curves
^^^^^^^^^^^^^^^^^^^^

Two helper classes are used by the spatial indexes and are also available for
use by codes.  
The `Morton index <https://en.wikipedia.org/wiki/Z-order_curve>`_, or Z-order
curve, maps each point in a 3D bounding box to a point along a 1D space filling
curve.  The RectangularLattice maps each point in a 3D bounding box to an
integer-addressed cell.

**Contents:**

.. toctree::
   :maxdepth: 2

   first_example
   primitives
   operations
   indexes

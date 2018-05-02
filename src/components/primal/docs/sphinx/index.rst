
Primal User Documentation
=========================

The Primal component of Axom provides tools for computational geometry.  Classes
in Primal include representations of geometric primitives and spatial indexes;
functions implement geometric queries such as intersection, closest-point, and
squared distance.  Primal's goal is to provide an efficient implementation of 
fundamental geometric ideas.  Primal aims to be useful within Axom as well as 
for application codes that use shapes and geometric operations.

Development in Primal has focused on addressing demand from application codes rather
than completeness, and the shapes, operations, and spatial indexes available in
Primal reflect this priority.  If you need a new one, we will build it.

Primitives
----------

Primal includes the following primitives:
- Point
- Segment, Ray, Vector
- Plane, Triangle, Polygon
- Sphere
- BoundingBox, OrientedBoundingBox

The classes are implemented using a NumericArray class.  NumericArray represents a
tuple of numbers and provides arithmetic and output operators.  Classes in Primal
are templated on coordinate type (double, float, etc.) and dimension.  The
primitives do not inherit from a common base class.  This was a design choice
in favor of simplicity and performance.  Geometric primitives each have equality
and inequality operators and an output left-shift operator.  


Operations
----------

Primal implements geometric operations with functions, listed below:
- ``closest_point`` takes a primitive P and a query point Q, and returns the point P
  closest to Q.
- ``intersect`` predicate tests if two primitives intersect.  Some of the combinations
  also indicate the point of intersection of a 1D and a 2D primitive.
- ``clip`` finds the polygon resulting from a bounding box clipping a triangle.
- ``orientation`` finds the side of a line segment or triangle where a query point lies.
- ``squared_distance`` computes the squared distance from a point to another point,
  bounding box, line segment, or triangle.
- ``compute_bounding_box`` finds the bounding box for a given primitive.

File elsewhere:
- Morton functions (list with Morton spatial index)
- Rectangular lattice functions (list below with class)
- merge_boxes (list with bounding box and OBB objects)
- compute_oriented_bounding_box (list with OBB object)

Most use cases have low dimension, usually 2 or 3.  Dimensionality has been
generalized to support other values where it does not interfere with the common case, 
but some operations such as triangle intersection do not support other dimensionality
than 2 or 3.  

Spatial Indexes
---------------

There are two spatial indexes provided by Primal, represented by the BVHTree and
UniformGrid classes.  These are supported by the RectangularLattice class and
several classes implementing a Morton index.

The `bounding volume hierarchy tree <http://foo.boo.com>`_ is implemented by BVHTree.
A code will add objects with associated bounding boxes to the BVHTree, then call
the ``build()`` method.  This method hierarchically subdivides a region of space
until the smallest subdivisions (leaf buckets) intersect no more than 25 objects.

The `Verlet list <https://en.wikipedia.org/wiki/Cell_lists>`_ or cell list is
implemented by UniformGrid.  As with BVHTree, a code adds objects with bounding
boxes to the UniformGrid, then calls ``build()``.  Unlike BVHTree, UniformGrid
divides the region space into one layer of non-overlapping bins, all of equal size.

Both BVHTree and UniformGrid allow a code to get the list of bins that intersect a query
bounding box, then perform an action on each candidate object from the list of bins.

Two helper classes are used by the spatial indexes and are also available for use by codes.
The `Morton index <https://en.wikipedia.org/wiki/Z-order_curve>`_, or Z-order curve, is
implemented by several classes and functions within Primal.  It maps each point in a 3D
bounding box to a point along a 1D space filling curve.  The RectangularLattice maps
each point in a 3D bounding box to an integer-addressed cell.

**Contents:**

.. toctree::
   :maxdepth: 2

   first_example
   primitives
   operations
   indexes

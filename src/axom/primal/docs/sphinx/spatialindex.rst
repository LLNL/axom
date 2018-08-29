Spatial Query Acceleration Data Structures
------------------------------------------

Primal provides two data structures for accelerating spatial queries, sometimes
called "spatial indexes."  The cell list is
implemented in the UniformGrid class and the bounding volume hierarchy tree
implemented in the BVHTree class.  Both classes divide an axis-aligned 
bounding box denoting a
region of interest into bins that group objects together, avoiding the need to
process objects that do not fall into a bin of interest.  The UniformGrid and
BVHTree classes are supported by the RectangularLattice class and construction
functions, and several classes implementing a Morton index.

.. toctree::
   :maxdepth: 2

   idxuniformgrid
   idxbvhtree


Space filling curves
--------------------

Two helper classes are used by the UniformGrid and BVHTree classes
and are also available for use by codes.
The `Morton index <https://en.wikipedia.org/wiki/Z-order_curve>`_, or Z-order
curve, maps each point in a 3D bounding box to a point along a 1D space filling
curve.  The RectangularLattice maps each point in a 3D bounding box to an
integer-addressed cell.

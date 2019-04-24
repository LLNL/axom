SpatialOctree
^^^^^^^^^^^^^

Axom provides an implementation of the octree spatial index.  The
``SpatialOctree`` recursively divides a bounding box into a hierarchy of
non-intersecting bounding boxes.  Each level of subdivision divides the bounding
box of interest along each of its dimensions, so 2D ``SpatialOctree`` objects
contain four child bounding boxes at each level, while 3D objects contain eight
children at each level.  The process stops when each a bounding box contains
less than a specified number of objects or when the tree reaches a specified
height.

The ``Octree`` class hierarchy is useful for building custom spatial acceleration
data structures, such as ``quest::InOutOctree``.

The figure below shows the construction of several levels of a 2D octree.

.. figure:: figs/showOctree0.png
   :figwidth: 300px
   :alt: Diagram showing points in a bounding box

.. figure:: figs/showOctree1.png
   :figwidth: 300px
   :alt: Diagram showing first division of a SpatialOctree

.. figure:: figs/showOctree2.png
   :figwidth: 300px
   :alt: Diagram showing second division of a SpatialOctree

.. figure:: figs/showOctree3.png
   :figwidth: 300px
   :alt: Diagram showing third division of a SpatialOctree

The following code example shows use of the ``SpatialOctree``.  (Placeholder
showing BVHTree; SpatialOctree needs to be written.)

.. literalinclude:: ../../examples/spin_introduction.cpp
   :start-after: _bvhtree_candidate_start
   :end-before: _bvhtree_candidate_end
   :language: C++

Some ancillary classes used in the implementation of ``SpatialOctree`` include
``BlockData``, which ties data to a block, ``Brood``, used to organize sibling
blocks, ``OctreeBase``, implementing non-geometric operations such as refinement
and identification of parent or child nodes, and ``SparseOctreeLevel`` and
``DenseOctreeLevel``, which hold the blocks at any one level of the
``SpatialOctree``.

The Morton index (link), implemented by ``Mortonizer`` and related classes,
translates an N-D point into a point on a space-filling 1D curve.  The
``RectangularLattice`` class (link) establishes a regular grid with
a block spacing for each dimension, and provides the integer coordinates of
any point in this regular grid.

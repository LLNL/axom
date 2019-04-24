RectangularLattice
^^^^^^^^^^^^^^^^^^

The ``RectangularLattice`` is a helper class that maps all of N-D space into
a regular, rectangular grid of cells identified by integer coordinates.  The
grid is defined by an origin point and a spacing for each dimension.

The following example shows the use of the ``RectangularLattice`` (placeholder,
replace when ready):

.. literalinclude:: ../../examples/spin_introduction.cpp
   :start-after: _bvhtree_candidate_start
   :end-before: _bvhtree_candidate_end
   :language: C++

Mortonizer
^^^^^^^^^^

The ``Mortonizer``, along with its associated ``MortonBase`` and ``PointHash``
classes, implements the Morton index, an operation that associates each point
in N-D space with a point on a space-filling curve.

The Z-index curve is a fractal space-filling curve.  The first four iterations
of the 2D curve's construction are shown below, and a 3D analogue exists as well.

(insert figures)

When properly templated, the ``Mortonizer`` maps an N-D point onto the Z-index
curve, and can also reverse the mapping.  This is shown in the following
example:

(insert example)

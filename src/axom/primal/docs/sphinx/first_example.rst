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

.. toctree::
   :maxdepth: 2

   opclip
   opclosestpoint
   opbbox
   opintersect
   oporientation
   opsqdist

Spatial Index
-------------

A spatial index is a data structure used to speed up retrieval of geometric
objects.  Primal provides two spatial indices, the cell (or Verlet) list
implemented in the ``UniformGrid`` class and the bounding volume hierarchy tree
implemented in the ``BVHTree`` class.  Both classes divide a bounding box
denoting a region of interest into bins that group objects together, avoiding
the need to process objects that do not fall into a bin of interest.

.. toctree::
   :maxdepth: 2

   idxuniformgrid
   idxbvhtree



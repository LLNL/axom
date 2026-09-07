.. ## Copyright (c) Lawrence Livermore National Security, LLC and other
.. ## Axom Project Contributors. See top-level LICENSE and COPYRIGHT
.. ## files for dates and other details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

Slam user guide
===============

Slam is Axom's Set-theoretic Lightweight API for Meshes. It provides building
blocks for developers who implement mesh data structures in simulation codes.
Applications choose the mesh topology, storage and parallel decomposition.
Slam supplies the types used to describe and access their entities (sets),
connectivity (relations) and fields (maps).

Start with a mesh
-----------------

The :doc:`introductory example <first_example>` builds a quadrilateral mesh, traverses its connectivity
and computes fields. It uses three kinds of objects:

* A **set** identifies mesh entities, such as vertices or cells.
* A **relation** records connections between two sets, such as the vertices of each cell.
* A **map** attaches values to set elements, such as a temperature at each vertex.

The :doc:`core_concepts` page explains their indexing rules.

Choose a representation
-----------------------

A quadrilateral has four vertices, so the cell-to-vertex incidence relation
for a quad mesh can use a compile-time cardinality. The reverse relation
usually needs a different number of cells at each vertex.
Slam's policies express these choices, along with storage, offsets and strides.

Slam provides several common :doc:`type aliases <aliases>` and construction helpers.
For custom storage or a less common configuration, you can customize
type configurations via policies as described in :doc:`implementation_details`. 
Ownership of the underlying storage depends on the containing type:
a default map owns its value buffer, while a static relation borrows its
connectivity buffers and sets.

Details about using Slam on GPUs and other execution spaces 
is described in :doc:`portability`.

API documentation
-----------------

The `Doxygen API documentation <../../../../doxygen/html/slamtop.html>`_
contains class and function reference material. :doc:`examples` points to
larger mesh examples in the source tree.

.. toctree::
   :caption: Contents
   :maxdepth: 2

   first_example
   core_concepts
   aliases
   implementation_details
   portability
   examples
   more

.. ## Copyright (c) Lawrence Livermore National Security, LLC and other
.. ## Axom Project Contributors. See top-level LICENSE and COPYRIGHT
.. ## files for dates and other details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _srm-label:

=============
Core concepts
=============

This section describes Slam's concepts: what they mean and how they are used.

.. figure:: figs/set_relation_map.png
   :figwidth: 400px
   :alt: Sets, relations and maps in slam
   :align: center

   A **relation** (blue lines) between two **sets** (ovals with red and green dots, as elements)
   and a **map** of scalar values (brown) on the second set.

.. _set-concept-label:

Set
===

Sets model mesh entities such as vertices, zones, materials or refinement levels.
Ordered sets associate each entity with a position so Slam code can iterate
and index efficiently.

Use ``RangeSet`` for contiguous ranges. Use ``ArraySet`` or ``ArrayViewSet``
when set elements are stored in Axom buffers.

.. Future
   Discuss different indexing schemes for ProductSets


.. _relation-concept-label:

Relation
========

A relation connects elements from a `from-set` to those of a `to-set`
and can be used to encode mesh incidence and adjacency relations, such 
as those from cells to vertices, from vertices to cells, from zones to materials
or from elements to neighboring elements.

Slam classifies relations along a few independent axes:

* cardinality: constant per from-set entity or variable per from-set entity
* mutability: static after construction or dynamically editable
* storage: implicit, such as a product set, or explicit, such as index buffers

Use ``ConstantRelation`` / ``ConstantRelationView`` for static fixed-cardinality
relations and ``VariableRelation`` / ``VariableRelationView`` for static
CSR-shaped relations. Use ``DynamicConstantRelation`` or
``DynamicVariableRelation`` when the connectivity needs to be edited.


.. _map-concept-label:

Map
===

Maps attach values to the members of a set. In mesh terms, a map is the Slam
abstraction for scalar, vector or tensor data associated with vertices, cells,
materials or other entity sets.

A ``Map`` stores its values in an ``axom::Array`` by default, allocating and freeing
that buffer as part of the map's lifetime. To point a map at a buffer whose lifetime is
managed elsewhere, e.g. storage owned by an application, or data to be captured in a device kernel,
construct it with an ``axom::ArrayView`` indirection instead.
``BivariateMap`` attaches values to a bivariate set, such as a product set or relation set,
with the same choice of storage.

See :ref:`aliases-label` for how the map storage default interacts with
those of other Slam containers, and for the relation aliases.

Set positions and elements
--------------------------

The set whose elements receive values is the map's mathematical *domain*.
For ``Map``, this is the set returned by ``set()``. For ``BivariateMap``, it is
the bivariate set, whose elements are pairs of positions in its first and
second sets.

Map access uses positions. If a set contains elements ``{10, 20, 30, 40}``,
``map(1)`` accesses the value associated with element ``20``. ``map.index(1)``
returns that element. With several components per element, ``map(1, c)``
accesses component ``c`` of the same entry. ``map[1 * map.numComp() + c]``
accesses that component by its flat storage position.

``BivariateMap::index(flatPosition)`` returns the coordinate pair stored at
that flat position in the bivariate set. Its existing two-argument overload,
``index(firstPosition, secondPosition)``, searches for the position within
the selected row.

Submaps
-------

A ``SubMap`` selects entries from a parent map and accesses their values in
the parent. It stores the parent pointer and a set of selected parent positions.
Its component count and shape come from the parent.

``submap.set()`` returns those selected positions. ``submap.index(i)`` follows
the selection to the set element associated with the value. Suppose the parent
set contains ``{10, 20, 30, 40}`` and the submap selects positions ``{3, 1}``:

.. list-table:: Positions and elements in a submap
   :header-rows: 1

   * - Submap position ``i``
     - Parent position ``submap.set()->at(i)``
     - Selected element ``submap.index(i)``
   * - 0
     - 3
     - 40
   * - 1
     - 1
     - 20

A nested submap selecting position ``1`` of this submap therefore accesses
the parent's value for element ``20``. Each level uses the same rule:

.. code-block:: cpp

   submap.index(i) == parent.index(submap.set()->at(i))

For a submap of a bivariate map, the result is a coordinate pair. This rule
holds at every nesting depth.

A const submap preserves the parent's value access. A submap of a mutable
owning map permits writes even through a const submap wrapper. A submap of a
const owning map returns const references. View-backed parents retain their
own const behavior.

The parent map and any buffers referenced by the selected-position set must
outlive the submap and its iterators. SubMap copies the selected-position set
itself, so a temporary range of positions is safe. If the parent is reassigned,
the selected positions must still be valid. Recreate iterators after changing
the parent's storage or component shape.


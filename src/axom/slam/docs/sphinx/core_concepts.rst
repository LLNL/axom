.. ## Copyright (c) Lawrence Livermore National Security, LLC and other
.. ## Axom Project Contributors. See top-level LICENSE and COPYRIGHT
.. ## files for dates and other details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _srm-label:

=============
Core concepts
=============

A mesh can use a **set** for its cells, a **relation** for the vertices of each cell,
and a **map** for a temperature field. The same models apply to materials,
refinement levels and other entities. This page explains what each model
means, how its positions work, and what its C++ concept requires.

.. figure:: figs/set_relation_map.png
   :figwidth: 400px
   :alt: Sets, relations and maps in slam
   :align: center

   The blue lines form a relation between the two sets of dots.
   The brown values form a map on the second set.

.. _set-concept-label:

Set
===

Sets model mesh entities such as vertices, zones, materials or refinement levels.
Ordered sets associate each entity with a position for indexing and traversal.
A position is an integer used to access the set. The element at that position
may be an integer identifier, an application-specific handle or a coordinate pair.

Use ``PositionSet`` for integers starting at zero, ``RangeSet`` for a contiguous
range with an offset, and ``ArrayIndirectionSet`` or
``ArrayViewIndirectionSet`` for (non-contiguous) an ordered collection 
of indexes using an indirection array.

``ProductSet`` represents ordered pairs of positions from two sets.
For example, a product of cell and material sets includes every cell-material
pair. ``at(flatPosition)`` accesses a pair by a single position, while
``getElements(firstPosition)`` returns a subset of second-set positions.
Referenced sets must outlive the product and retain their sizes while it is used.
See :ref:`adapter-contracts-label` for its representation limits.


.. _relation-concept-label:

Relation
========

A relation connects a *from-set* to a *to-set*. For each from-set entry,
it provides the related to-set entries. A cell-to-vertex boundary relation,
for example, gives the vertices associated with each cell. Reversing it gives a
vertex-to-cell co-boundary relation identifying the cells indident in each vertex.

We choose a relation by its cardinality and how its connections change:

* Constant cardinality gives every from-set element the same number of related entries.
  Variable cardinality allows that number to differ between elements.
* Static relations bind existing connectivity buffers. Dynamic relations
  provide operations to edit connectivity.

Use ``ConstantRelation`` / ``ConstantRelationView`` for static fixed-cardinality
relations and ``VariableRelation`` / ``VariableRelationView`` for static
variable-cardinality relations. Use ``DynamicConstantRelation`` or
``DynamicVariableRelation`` when the connectivity needs to be edited.

``RelationSet`` presents a static relation as a bivariate set of coordinate pairs.
``ProductSet`` represents all pairs without an explicit connectivity buffer.
Both model the ``BivariateSetLike`` concept, which is a different C++ contract
from the per-element connectivity required by ``RelationLike``.


.. _map-concept-label:

Map
===

Maps attach values to set elements. For example, a temperature field might have
one value per vertex, while a velocity field has several components, 
and a tensor field has a component shape. A component can also hold
an application-defined value type, such as the point class
in our :doc:`introductory example <first_example>`.

A ``Map`` stores its values in an ``axom::Array`` by default, allocating and freeing
that buffer as part of the map's lifetime. To access a buffer managed elsewhere,
construct the map with ``policies::ArrayViewIndirection`` instead. 
``BivariateMap`` attaches values to a bivariate set, such as a product set
or relation set, with the same choice of storage.

See :ref:`aliases-label` for storage and ownership choices.

Sizes and component shapes
--------------------------

A bound map has a positive number of components per set element, even when
the set is empty. Every dimension of a tensor component shape must be positive.
Custom stride policies must obey the same contract. The scalar count must fit
both the map's position type and ``axom::IndexType``. A supplied ``ArrayView``
must have exactly that many entries. For raw-pointer helpers,
the caller remains responsible for providing sufficient storage.

These are map restrictions, not requirements on every use of a stride policy.
For example, ordered sets may use a negative stride to traverse a range in reverse.

The referenced set's size must remain consistent with the map's value buffer.

Value constness follows the backing storage. A const owning map returns const references.
A const map backed by ``ArrayView<T>`` returns ``T&``, while one
backed by ``ArrayView<const T>`` returns ``const T&``. Direct access, scalar
iterator access, and submaps preserve these reference types.

Set positions and elements
--------------------------

The set whose elements receive values is the map's mathematical *domain*.
For ``Map``, this is the set returned by its ``set()`` method. 
For ``BivariateMap``, it is the bivariate set, whose elements are pairs
of positions in its first and second sets.

Map access uses positions. If a set contains elements ``{10, 20, 30, 40}``,
``map(1)`` accesses the value associated with element ``20``. ``map.index(1)``
returns that element. With several components per element, ``map(1, c)``
accesses component ``c`` of the same entry. ``map[1 * map.numComp() + c]``
accesses that component by its flat storage position.

``BivariateMap::index(flatPosition)`` returns the coordinate pair stored at
that flat position in the bivariate set. Its two-argument overload,
``index(firstPosition, secondPosition)``, searches for the position within the selected subset.

Submaps
-------

A ``SubMap`` selects entries from a parent map and accesses their values in the parent.
It stores the parent pointer and a set of selected parent positions.
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

For a submap of a bivariate map, the result is a coordinate pair.

The parent map and any buffers referenced by the selected-position set must
outlive the submap and its iterators. ``SubMap`` copies the selected-position
set itself, so a temporary ``RangeSet`` of positions is safe. Copying a set
that borrows a buffer does not extend the buffer's lifetime. If the parent is
reassigned, the selected positions must still be valid. Recreate iterators after changing
the parent's storage or component shape.

C++ concept contracts
=====================

``Concepts.hpp`` defines the operations generic code can use.

.. list-table:: Public object contracts
   :header-rows: 1
   :widths: 25 75

   * - Concept
     - Required operations and meaning
   * - ``SetLike<S>``
     - A signed integral ``PositionType``, an ``ElementType``, and const
       ``size()``, ``empty()``, and ``at(position)``. ``size()`` returns
       ``PositionType``, ``empty()`` converts to ``bool``, and ``at()`` converts
       to ``ElementType``. The element type need not be integral.
   * - ``IterableSetLike<S>``
     - ``SetLike`` plus const ``begin()`` and ``end()``. Traversal supports
       dereference, increment, and comparison with the end, and visits elements
       in positional order. Iterator aliases and standard-range conformance
       are not required.
   * - ``BivariateSetLike<S>``
     - ``SetLike`` with first and second sets that also model ``SetLike``.
       Coordinate members have their respective sets' position types.
       ``getFirstSet()`` and ``getSecondSet()`` return those sets, and
       ``getElements(first)`` returns a sized, iterable subset of second-set
       positions.
   * - ``RelationLike<R>``
     - From and To sets that model ``SetLike``, returned by ``fromSet()`` and
       ``toSet()``. ``relation[from]`` returns a sized, iterable collection of to-set
       positions. Neither set is restricted to scalar elements.
   * - ``MapLike<M>``
     - Signed integral ``PositionType``, component ``DataType``, ``size()``,
       ``numComp()``, ``index(position)``, and mutable and const
       ``flatValue(position, component)``. Value access returns lvalue
       references to ``DataType``, possibly adding constness.
   * - ``MapOver<M, S>``
     - ``MapLike`` with an explicit whole-set binding to exactly ``S``.
       ``MappedSetType`` names ``S``, ``set()`` returns a pointer to const ``S``,
       and the map and set use the same position type. ``index(position)``
       converts to the set's element type.
   * - ``Validatable<T>``
     - Const ``isValid(false)`` returns a result convertible to ``bool``.
       This is an optional capability, independent of the object concepts above.

Concept semantics
-----------------

A concept checks expressions and types, not runtime consistency.
A set's ``size()`` must be nonnegative and agree with ``empty()``.
Access requires valid positions. Dynamic containers can retain invalid entries
within their extent and their validity API determines which positions may be accessed.

Bivariate coordinates must identify valid positions in both component sets.
Traversing the related subsets in first-set order must agree with ``at(flatPosition)``.
The flat cardinality equals the sum of the subset sizes. A relation must contain
valid to-set positions, and each related collection's size must agree with its traversal.
Search methods and flat-position projections are optional conveniences,
not requirements of these public concepts.

Map access and set binding
--------------------------

``flatValue(entry, component)`` uses a local scalar-component offset in ``[0, numComp())``,
even for tensor maps. For a ``Map`` with shape ``{2, 3}``, ``flatValue(entry, 5)``
and ``value(entry, 1, 2)`` access the same component. In a bivariate map,
``entry`` is a single flat position rather than a pair of first-set and second-set
positions. It is not the global component-storage position used by ``operator[]``.

``MapLike`` does not require a ``set()`` binding. ``index(entry)`` identifies
the element receiving the values, and the component count is positive for a
bound map. A default, unbound ``SubMap`` may be empty with zero components.
Writable access is a separate requirement for algorithms that modify values.
A const owning map is read-only, while a const view can retain mutable access
to its backing storage.

``Map``, ``DynamicMap``, and ``BivariateMap`` provide ``MappedSetType`` to name
the set returned by ``set()``. For ``BivariateMap``, this is the bivariate set,
not the internal flat ``SetType``.
For a valid bound ``MapOver<M, S>``, ``size() == set()->size()`` and
``index(i)`` identifies the same element as ``set()->at(i)``.

``SubMap`` models ``MapLike`` but not ``MapOver``. Its ``set()`` selects positions
in its parent. Those positions are not the mapped elements returned by ``index()``.
Generic algorithms that only read or write mapped values should require ``MapLike``.
Require ``MapOver`` when the algorithm needs the whole-set binding.

Contract boundaries
-------------------

The public concepts describe operations that generic algorithms can require.
They do not promise that every matching type works with every Slam class.
For example, ``RelationSet`` also needs a relation's flat storage. Those
additional requirements belong to the adapter, not to ``RelationLike``.
See :ref:`adapter-contracts-label` for the operations consumed by
``RelationSet``, ``BivariateMap`` and ``SubMap``.

Size, stride, offset, subset, ordered-set indirection, and map-storage policy protocols
remain separate extension contracts. Their owner-specific combinations do not
enter the public object concepts. See :ref:`policy-contracts-label` for the
operations each owner uses.

``PositionLike`` means a signed integral type. An application-specific handle
can be a set element, but it cannot replace the integer position type through
an opt-in declaration. Slam's indexing arithmetic requires signed integers.

``TriviallyCopyableRepresentation`` checks only the C++ representation.
It does not certify device-callable operations, device-accessible allocations,
or the lifetime of referenced objects. In particular, a trivially copyable
``SubMap`` can retain a pointer to a host-only parent. See :doc:`portability`
before using a type in a kernel.

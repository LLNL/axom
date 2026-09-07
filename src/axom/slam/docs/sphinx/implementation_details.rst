.. ## Copyright (c) Lawrence Livermore National Security, LLC and other
.. ## Axom Project Contributors. See top-level LICENSE and COPYRIGHT
.. ## files for dates and other details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

======================
Implementation details
======================

.. _policy-label:

Policy-based design
===================

Slam combines policies to select a type's storage and indexing behavior.
For example, an ordered set can store its elements in an array or compute
them from an offset and stride.

Start with the common :doc:`alias types <aliases>` and use policies directly
when you need custom storage or a configuration the aliases do not cover.

Some examples:

* Size, stride and offset policies select compile-time or runtime values.
* Indirection policies choose how set elements, relation indices and map
  values are reached through backing storage:

  * ``NoIndirection`` computes an element from the position after applying
    the set's offset and stride. No element buffer is needed.
  * ``ArrayIndirection``, backed by an ``axom::Array``. Maps own their array.
    Ordered sets and static relations bind an array managed elsewhere.
  * ``ArrayViewIndirection``, backed by an ``axom::ArrayView`` of a buffer managed elsewhere.
    Device use also requires device-callable operations and accessible backing storage.
  * ``CArrayIndirection`` and ``STLVectorIndirection`` for interoperation
    with raw-pointer and ``std::vector`` storage.
  * Custom policies can adapt other buffer types to the required operations.
* Subsetting policies select no parent, a virtual parent set, or a concrete parent set.

``Map`` and ``BivariateMap`` default to ``ArrayIndirection``. Buffer ownership
comes from the containing type and its buffer representation, not from a
separate ownership policy.

The following feature diagram of ``slam::OrderedSet`` policies shows how these policies
combine in the subscript operator:

.. figure:: figs/orderedset_feature_diagram.png
   :figwidth: 100%
   :alt: Feature diagram for slam's ordered set


.. _policy-contracts-label:

Policy extension contracts
==========================

Policy concepts describe operations, rather than membership in Slam's
built-in policy families. Custom policies do not need to inherit ``ValuePolicy``
or provide ``TagType``, ``IntType``, or device-accessibility flags.

``SizePolicy``, ``StridePolicy``, ``OffsetPolicy``, and ``SubsetPolicy`` describe
const queries and ignore top-level cv/ref qualification.

Ordered sets
------------

``OrderedSet`` constructs its scalar policies from position values, copies them
into the set, and default-constructs and assigns them in its builder.
Its iterators also copy the policies. Each scalar policy supplies
``DEFAULT_VALUE`` of exactly the set's position type.

.. list-table:: Operations consumed by ordered-set policies
   :header-rows: 1
   :widths: 20 80

   * - Policy
     - Const operations
   * - Size
     - ``size()``, ``empty()``, ``isValid(verbose)``
   * - Offset
     - ``offset()``, ``isValid(verbose)``
   * - Stride
     - ``stride()``, ``isValid(verbose)``
   * - Indirection
     - ``indirection(pos)``, ``isValid(size, offset, stride, verbose)``
   * - Subset
     - ``isSubset()``, ``parentSet()``, ``isValid(begin, end, verbose)``

Size, offset, and stride queries return the set's position type. A scalar
ordered-set stride does not need map shape metadata. ``DynamicSet`` additionally
requires mutable ``size()`` to return a position reference, since insertion
and reset update the stored size.

``SubsetPolicy`` checks ``ParentSetType`` and the parent queries. ``OrderedSet``
also requires default construction, copying, and construction from a parent
pointer. Its validation checks ``isValid(begin, end, verbose)`` with the actual
const iterator type. The parent and any bound indirection storage must outlive the set and its
iterators. Relation policies can require additional operations, such as flat
buffer access, which ordered-set indirection does not promise.

Maps
----

A map stride policy supplies ``IndexType``, ``ShapeType``, a positive compile-time
``NumDims``, ``DefaultSize()``, and construction from a shape. Const ``stride()``
returns the scalar component count and ``shape()`` returns the shape. A
one-dimensional shape converts to ``IndexType``, while a multidimensional shape
supports subscripting, and ``strides()`` returns its per-dimension strides.

``Map`` does not need the built-in shape container, shape iterators, a default
stride constructor, or a policy ``isValid()`` method. Custom policies must
provide positive dimensions, a representable product equal to ``stride()``,
and row-major per-dimension strides. ``Map`` checks the scalar storage size.
The policy must perform its own shape arithmetic safely.

.. _adapter-contracts-label:

Additional type requirements
============================

Some Slam classes need more operations than the public object concepts
require. These checks live in ``detail``. They constrain the particular
consumer without adding requirements to every set, relation or map.

ProductSet
----------

``ProductSet`` requires a signed integral flat position type that can represent
the position types of both constituent sets. Construction checks that their
sizes are nonnegative and that the product fits the flat position type.
Either set may be empty.

BivariateMap
------------

``BivariateMap`` needs its bivariate set to support search and conversion of
flat positions for a selected subset to its concrete ``RangeSet``. It derives the subset type from
``getElements()`` and iterator coordinates from ``at()``. The set does not
need an iterator or a ``SubsetType`` alias.


RelationSet
-----------

``RelationSet`` obtains related subsets from const ``relation[from]``. It searches the
flat storage using ``offset(from)`` and the subset's size, and projects flat
positions through ``firstIndex(flat)`` and ``relationData()[flat]``. These
operations must agree with subset traversal. The subset needs no subscript or
offset operation, and the relation needs no ``RelationSubset`` alias.

A concrete ``RelationSet`` returns the relation's subset directly. A virtual
``RelationSet`` requires conversion to the fixed subset type of ``BivariateSet``.
If that conversion is unavailable, the concrete adapter's ``VirtualSet``
and ``OtherSet`` aliases are ``void``.

SubMap
------

``SubMap`` construction requires scalar access, size, component count, shape
and element lookup from its parent. Traversal also requires a parent iterator
that can be constructed at a selected position, dereferenced, copied and
queried for its original flat position. The parent needs no ``set_end()``.

The selected-position set is copied into the submap. It supplies positional
access and does not need its own iteration or subscript operations. Shaped
component access is checked only when requested. See :ref:`srm-label` for
the selection rule and parent lifetime requirements.

.. _setup-label:

Simplifying mesh setup
======================

Use ``make_map`` and the ``make_*_relation`` helpers when the sets and buffers
are already available. The helpers deduce template arguments and construct
the corresponding Slam types. :doc:`first_example` uses both constant- and
variable-cardinality relation helpers.

Builders let you name individual construction parameters. For example,
this range contains the elements ``10``, ``11`` and ``12``:

.. code-block:: cpp

   using IDs = axom::slam::RangeSet<int>;
   IDs ids = IDs::SetBuilder().size(3).offset(10);

A builder configures the chosen type. It does not change its policies or take
ownership of borrowed storage. The sets and buffers passed to construction
helpers must still satisfy the resulting object's lifetime requirements.

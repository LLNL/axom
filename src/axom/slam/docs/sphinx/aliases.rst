.. ## Copyright (c) Lawrence Livermore National Security, LLC and other
.. ## Axom Project Contributors. See top-level LICENSE and COPYRIGHT
.. ## files for dates and other details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _aliases-label:

===============================
Choosing policies and aliases
===============================

A Slam set, relation or map is assembled from policy template parameters including:
cardinality, stride, indirection, offset, size, subsetting and interface.
Choosing those policies is how you describe the data structure, i.e. how its
connectivity is shaped, where its data lives, and what is fixed at compile time.
This page explains the choices that come up most often, and a small set of
aliases in ``axom/slam/Aliases.hpp`` that name the most common relation configurations.

The aliases are a convenience for these common types and are preferred when applicable.
Use the policies directly for configurations not covered by an alias, such as indirection buffers
backed by ``std::vector`` or C-arrays, or specialized cardinalities.

Buffer binding and ownership
============================

An indirection policy selects how values are accessed. The containing Slam
type determines whether the buffer is owned or borrowed.

``policies::ArrayIndirection``
  An ordered set or static relation stores a pointer to an external
  ``axom::Array``. That array object must outlive the binding.
  A ``Map`` instead stores the policy's buffer type by value, so a map using
  this policy owns its array. ``BivariateMap`` uses such a map internally.

``policies::ArrayViewIndirection``
  Sets, static relations, and maps store an ``axom::ArrayView`` by value.
  They do not own its allocation. The allocation must remain valid while
  the view is used, including during asynchronous execution.

For example, this map owns its values but borrows its set:

.. code-block:: C++

   using Cells = slam::RangeSet<>;
   Cells cells(numCells);
   slam::Map<double, Cells> density(&cells);

An ``ArrayView`` binding does not make the entire object device-usable.
The interface and operations must support device execution, and every
referenced set and buffer must be accessible there.

``policies::STLVectorIndirection`` and ``policies::CArrayIndirection`` remain
available for adapting ``std::vector`` and raw-pointer storage. The vector
policy is host-only. As with ``ArrayIndirection``, sets and static relations
borrow the vector object; a vector-backed map holds its vector by value.


Fixing parameters at compile time
==================================

Where a quantity is known at compile time, encoding it in a policy lets the
compiler specialize the generated code and data structures.

The stride of a relation or map is the common case. A quad mesh's
element-to-vertex relation has exactly four vertices per element,
so its stride can be a compile-time constant:

.. code-block:: C++

   using ElemVertRelation =
     slam::ConstantRelation<ElemSet, VertSet, /* stride */ 4>;   // CompileTimeStride

When the same count is only known at runtime, use the runtime form
(``RuntimeConstantRelation``, or ``policies::RuntimeStride``) directly.
The same distinction applies to a set's size (``policies::CompileTimeSize`` vs. a runtime size) and offset.

The set and relation aliases
============================

A handful of relation configurations recur across mesh code and require spelling out all
six ``StaticRelation`` policy parameters, so a named shorthand can be helpful.
The aliases below cover them, together with the two indirection set types.
Each relation alias has a ``*View`` form that stores an ``axom::ArrayView``
instead of a pointer to an external ``axom::Array``. Both forms borrow storage
and their from/to sets. The suffix names the binding type rather than a conversion operation.

Sets:

``RangeSet<P,E>``
  A contiguous range of positions, with no separate storage. Use it for dense
  ranges of mesh entities such as cells, nodes, materials or levels.

``ArrayIndirectionSet<P,E>`` / ``ArrayViewIndirectionSet<P,E>``
  Sets defined in ``axom/slam/IndirectionSet.hpp``. The first binds an external
  ``axom::Array`` by pointer while the second stores an ``axom::ArrayView`` by value.
  Neither owns its elements.

Relations:

``ConstantRelation<FromSet, ToSet, N>``
  A static relation with exactly ``N`` to-set entities per from-set entity,
  with ``N`` fixed at compile time.

``RuntimeConstantRelation<FromSet, ToSet>``
  As above, but with the (constant) cardinality supplied at runtime.

``VariableRelation<FromSet, ToSet>``
  A static relation whose cardinality varies per from-set entity.
  It binds begin offsets and to-set positions in external buffers.

Every relation alias fixes its entry type to ``ToSet::PositionType``.
Entries identify positions in the to-set, not its element values.
The final optional template argument, ``FlatPosType``, controls flat-storage positions
and begin offsets. It defaults to a signed type that can represent both sets' positions.
It must represent from-set positions and the total number of stored entries;
it need not represent to-set positions, which have their own type.

For example, ``VariableRelation<FromSet, ToSet, std::int64_t>`` uses 64-bit flat positions
regardless of the to-set's position type. The compile-time constant form uses the same choice as
``ConstantRelation<FromSet, ToSet, N, std::int64_t>``.
Use ``make_variable_relation`` or the constant-relation helpers when the
buffers already exist and deduction is more convenient than naming the type.

The dynamic relation classes, ``DynamicConstantRelation`` and
``DynamicVariableRelation``, keep their own names because their connectivity can
be edited after construction.

There are intentionally no map aliases. 
A ``Map`` or ``BivariateMap`` already defaults to an ``axom::Array`` indirection with stride one, 
so the common cases read clearly on their own:

.. code-block:: C++

   slam::Map<double, Cells> density(&cells);            // one value per cell
   slam::BivariateMap<double, CellMatSet> volfrac(&cm); // one value per (cell, material)

The cases that vary a map, e.g. a view into a buffer managed elsewhere, or a runtime stride,
tend to be the cases where a fixed alias name would not capture the configuration cleanly:

.. code-block:: C++

   using CellPos = Cells::PositionType;

   // A field that views a buffer managed elsewhere, with a runtime component count:
   using ViewedField =
     slam::Map<double,
               Cells,
               slam::policies::ArrayViewIndirection<CellPos, double>,
               slam::policies::RuntimeStride<CellPos>>;

   // An std::vector-backed field, for interoperation with existing storage:
   using VectorField =
     slam::Map<double,
               Cells,
               slam::policies::STLVectorIndirection<CellPos, double>>;

FieldRegistry
=============

``FieldRegistry`` is a host-side convenience that keeps a set of named fields and buffers.
Its field type (``Registry::MapType``) is a ``slam::Map`` with the default ``axom::Array`` indirection,
and its buffer type (``Registry::BufferType``) is an ``axom::Array``:

.. code-block:: C++

   using CellSet = slam::RangeSet<>;
   using Registry = slam::FieldRegistry<CellSet, double>;

   CellSet cells(numCells);
   Registry fields;
   auto& density = fields.addField("density", &cells);    // Registry::MapType
   auto& buffer  = fields.addBuffer("tmp", cells.size()); // Registry::BufferType (axom::Array<double>)

Because the registry's buffers are now ``axom::Array`` rather than
``std::vector``, code that consumes them should:

* use ``auto&`` or ``Registry::BufferType&`` rather than ``std::vector<T>&``;
* use ``buffer.data()``, ``buffer.size()`` and ``buffer.resize(n)`` for the
  common vector-like operations;
* call ``buffer.view()`` when an ``axom::ArrayView`` is the right interface;
* register externally-managed storage with ``addFieldView()`` or ``addBufferView()``.

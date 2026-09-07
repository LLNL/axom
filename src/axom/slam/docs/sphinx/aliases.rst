.. ## Copyright (c) Lawrence Livermore National Security, LLC and other
.. ## Axom Project Contributors. See top-level LICENSE and COPYRIGHT
.. ## files for dates and other details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _aliases-label:

=============================
Choosing policies and aliases
=============================

A relation with four vertices per cell needs different cardinality information
from one with a variable number of neighbors. A map that owns its values needs
a different buffer representation from one that borrows application storage.
Policies express these choices. Aliases name common combinations so you do
not have to spell out every template argument.

Use the common aliases where they fit and defer to the policies directly for other
configurations, such as ``std::vector`` or raw-pointer indirection buffers.
``axom/slam/Aliases.hpp`` provides the static-relation aliases listed below.

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
referenced set and buffer must be accessible there. See :doc:`portability`.

Slam also provides ``policies::STLVectorIndirection`` and ``policies::CArrayIndirection``
for adapting ``std::vector`` and raw-pointer storage. The std::vector
policy is host-only. As with ``ArrayIndirection``, sets and static relations
borrow the vector object. A vector-backed map holds its vector by value.


Fixing parameters at compile time
=================================

Where a quantity is known at compile time, encoding it in a policy lets the
compiler specialize the generated code and data structures.

The stride of a relation or map is the common case. A quad mesh's
element-to-vertex relation has exactly four vertices per element,
so its stride can be a compile-time constant:

.. code-block:: C++

   using ElemVertRelation =
     slam::ConstantRelation<ElemSet, VertSet, 4>;

When the relation's common count is only known at runtime, use
``RuntimeConstantRelation``. The same distinction applies to other policies.
For example, a set can fix its size with ``policies::CompileTimeSize`` or
store it with ``policies::RuntimeSize``.

The set and relation aliases
============================

The relation aliases select the six ``StaticRelation`` template arguments
for common configurations. The set aliases below cover contiguous ranges
and externally stored elements.
Each relation alias has a ``*View`` form that stores an ``axom::ArrayView``
instead of a pointer to an external ``axom::Array``. Both forms borrow storage
and their from/to sets. The suffix names the binding type rather than a conversion operation.

Sets:

``RangeSet<P,E>``
  A contiguous range of element values, computed from an offset without a
  separate buffer. ``PositionSet<P,E>`` fixes that offset at zero.
  Use these types for dense ranges of cells, nodes, materials or levels.

``ArrayIndirectionSet<P,E>`` / ``ArrayViewIndirectionSet<P,E>``
  Sets defined in ``axom/slam/IndirectionSet.hpp``. The first binds an external
  ``axom::Array`` by pointer while the second stores an ``axom::ArrayView`` by value.
  Neither owns its elements.

Relations:

``ConstantRelation<FromSet, ToSet, N>``
  A static relation with exactly ``N`` to-set entities per from-set entity,
  with ``N`` fixed at compile time.

``RuntimeConstantRelation<FromSet, ToSet>``
  A static relation with the same cardinality for every from-set element, supplied at runtime.

Both constant-cardinality forms require a positive count, even when the
from-set is empty.

``VariableRelation<FromSet, ToSet>``
  A static relation whose cardinality varies per from-set entity.
  It binds begin offsets and to-set positions in external buffers.

Every relation alias fixes its entry type to ``ToSet::PositionType``.
Entries identify positions in the to-set, not its element values.
The final optional template argument, ``FlatPosType``, controls flat-storage positions
and begin offsets. It defaults to a signed type that can represent both sets' positions.
It must represent from-set positions and the total number of stored entries.
It need not represent to-set positions, which have their own type.

For example, ``VariableRelation<FromSet, ToSet, std::int64_t>`` uses 64-bit flat positions
regardless of the to-set's position type. The compile-time constant form uses the same choice as
``ConstantRelation<FromSet, ToSet, N, std::int64_t>``.
Use ``make_variable_relation`` or the constant-relation helpers when the
buffers already exist and deduction is more convenient than naming the type.

Map configurations
==================

``Map`` and ``BivariateMap`` default to ``axom::Array`` indirection
with one component per entry:

.. code-block:: C++

   slam::Map<double, Cells> density(&cells);            // one value per cell
   slam::BivariateMap<double, CellMatSet> volfrac(&cm); // one value per (cell, material)

Name the policies explicitly when the storage or component count differs:

.. code-block:: C++

   using CellPos = Cells::PositionType;

   // Borrowed storage with a runtime component count
   using ViewedField =
     slam::Map<double,
               Cells,
               slam::policies::ArrayViewIndirection<CellPos, double>,
               slam::policies::RuntimeStride<CellPos>>;

   // A field that stores its values in an std::vector
   using VectorField =
     slam::Map<double,
               Cells,
               slam::policies::STLVectorIndirection<CellPos, double>>;

FieldRegistry
=============

``FieldRegistry`` keeps named fields and buffers on the host.
``Registry::MapType`` is a ``slam::Map`` with the default ``axom::Array``
indirection, and ``Registry::BufferType`` is an ``axom::Array``:

.. code-block:: C++

   using CellSet = slam::RangeSet<>;
   using Registry = slam::FieldRegistry<CellSet, double>;

   CellSet cells(numCells);
   Registry fields;
   auto& density = fields.addField("density", &cells);    // Registry::MapType
   auto& buffer  = fields.addBuffer("tmp", cells.size()); // Registry::BufferType (axom::Array<double>)

Use ``auto&`` or ``Registry::BufferType&`` to refer to a registry buffer.
It supports ``data()``, ``size()`` and ``resize(n)``.
Call ``buffer.view()`` when you need an ``axom::ArrayView``.

Use ``addFieldView()`` or ``addBufferView()`` to register externally managed
storage. The registry retains a view, not ownership of the allocation.
The storage must outlive every use of that view, and a field's set must also
remain valid while the field is used.

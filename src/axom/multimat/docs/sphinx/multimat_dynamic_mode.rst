.. ## Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

******************************************************
Dynamic Mode
******************************************************

The distribution of materials in MultiMat is controlled by the
Cell-Material Relation (CMR). For many MultiMat use cases, this is set
once (static mode) and then fields are defined on the MultiMat object. MultiMat
also supports a dynamic mode that permits materials to move around in the mesh.

When creating a MultiMat object using the default constructor, it will default
to static mode using ``CELL_DOM`` data layout with sparse data. The data layout
argument later impacts the ``addEntry()`` and ``removeEntry()`` methods that modify
the MultiMat object's CMR in dynamic mode. For those methods, when MultiMat is
created with a ``CELL_DOM`` data layout, it means that the first argument to
``addEntry()`` will be a cell number and the second will be a material number.

To convert the MultiMat object to dynamic mode, call the ``convertToDynamic()``
method. This method changes some internal representations (including field 
organization) to better support dynamic modifications of the CMR. For example,
when changing to dynamic mode, ``SPARSE`` fields are converted to ``DENSE``
so further changes do not require field data to be reallocated/reorganized
again. The CMR is modified using calls to the ``addEntry()`` and ``removeEntry()``
methods.

.. code-block:: cpp

    constexpr int nmats = 3;
    constexpr int ncells = 9;

    axom::multimat::MultiMat mm;

    // Multimat initialization omitted

    // Switch to dynamic mode
    mm.convertToDynamic();

    // Add material 2 in zone 3 that was not there before.
    mm.addEntry(3, 2);

    // Remove material 1 in zone 5
    mm.removeEntry(5, 1);

    // Volume fraction updates omitted (iterate Volfrac field, set new values)


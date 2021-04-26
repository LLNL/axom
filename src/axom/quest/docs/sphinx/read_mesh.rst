.. ## Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _reading-mesh:

*****************
Reading in a mesh
*****************

Applications commonly need to read a mesh file from disk.  Quest provides the
``STLReader`` class, which can read binary or ASCII `STL`_ files, as well as the
``PSTLReader`` class for use in parallel codes.  STL (stereolithography)
is a common file format for triangle surface meshes.  The STL reader classes
will read the file from disk and build a ``mint::Mesh`` object.

.. _STL: https://en.wikipedia.org/wiki/STL_(file_format)

The code examples are excerpts from the file ``<axom>/src/tools/mesh_tester.cpp``.

We include the STL reader header

.. literalinclude:: ../../../../tools/mesh_tester.cpp
   :start-after: _read_stl_include1_start
   :end-before: _read_stl_include1_end
   :language: C++

and also the mint Mesh and UnstructuredMesh headers.

.. literalinclude:: ../../../../tools/mesh_tester.cpp
   :start-after: _read_stl_include2_start
   :end-before: _read_stl_include2_end
   :language: C++

For convenience, we use typedefs in the axom namespace.

.. literalinclude:: ../../../../tools/mesh_tester.cpp
   :start-after: _read_stl_typedefs_start
   :end-before: _read_stl_typedefs_end
   :language: C++

The following example shows usage of the STLReader class:

.. literalinclude:: ../../../../tools/mesh_tester.cpp
   :start-after: _read_stl_file_start
   :end-before: _read_stl_file_end
   :language: C++

After reading the STL file, the ``STLReader::getMesh`` method gives access to the
underlying mesh data.  The reader may then be deleted.


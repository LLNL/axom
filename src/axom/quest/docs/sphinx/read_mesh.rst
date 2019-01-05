.. ##
.. ## Copyright (c) 2017-18, Lawrence Livermore National Security, LLC.
.. ##
.. ## Produced at the Lawrence Livermore National Laboratory
.. ##
.. ## LLNL-CODE-741217
.. ##
.. ## All rights reserved.
.. ##
.. ## This file is part of Axom.
.. ##
.. ## For details about use and distribution, please read axom/LICENSE.
.. ##

*****************
Reading in a mesh
*****************

Applications commonly need to read a mesh file from disk.  Quest provides the
STLReader class, which can read binary or ASCII STL files, as well as the
PSTLReader class for use in parallel codes.  STL (stereolithography) is a common
file format for triangle surface meshes.  The STL reader classes will read the
file from disk and build a mint Mesh object.

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


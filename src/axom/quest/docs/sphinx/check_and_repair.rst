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

***********************
Check and repair a mesh
***********************

The STL file format specifies triangles without de-duplicating their
vertices.  Vertex welding is needed for several mesh algorithms that
require a watertight manifold.  Additionally, mesh files often contain
errors and require some kind of cleanup.

The following code examples are excerpts from the file
axom/src/tools/mesh_tester.cpp.

Quest provides a function to weld vertices within a distance of some
specified epsilon.  This function takes arguments
``mint::UnstructuredMesh< mint::SINGLE_SHAPE > **surface_mesh`` and
``double epsilon``, and modifies ``surface_mesh``.  In addition to
the mint Mesh and UnstructuredMesh headers (see previous page), we
include the headers declaring the functions for checking and repairing
surface meshes.

.. literalinclude:: ../../../../tools/mesh_tester.cpp
   :start-after: _check_repair_include_start
   :end-before: _check_repair_include_end
   :language: C++

The function call itself:

.. literalinclude:: ../../../../tools/mesh_tester.cpp
   :start-after: _check_repair_weld_start
   :end-before: _check_repair_weld_end
   :language: C++

One problem that can occur in a surface mesh is self-intersection.  A
well-formed mesh will have each triangle touching the edge of each of its
neighbors.  Intersecting or degenerate triangles can cause problems for some
spatial algorithms.  To detect such problems using Quest, we first make some
containers for the list of defects that might be found.

.. literalinclude:: ../../../../tools/mesh_tester.cpp
   :start-after: _check_repair_intersections_containers_start
   :end-before: _check_repair_intersections_containers_end
   :language: C++

Then, we call the function to detect self-intersections and degenerate
triangles.

.. literalinclude:: ../../../../tools/mesh_tester.cpp
   :start-after: _check_repair_intersections_start
   :end-before: _check_repair_intersections_end
   :language: C++

After calling ``findTriMeshIntersections``, ``collisions`` will record
the indexes of each pair of intersecting triangles and ``degenerate`` will
contain the index of each degenerate triangle.  The user code can then address or
report any triangles found.  Mesh repair beyond welding close vertices is
beyond the scope of the Quest component.

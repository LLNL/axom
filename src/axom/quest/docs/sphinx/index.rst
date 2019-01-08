.. ##
.. ## Copyright (c) 2017-19, Lawrence Livermore National Security, LLC.
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

Quest User Documentation
========================

The Quest component of Axom provides several spatial queries and operations on
a ``mint::Mesh``.

  - Read a surface mesh from an STL file
  - Check for some common mesh errors; deduplicate vertices

    - vertex welding: merge vertices closer than a specified distance
      "epsilon"
    - find self-intersections and degenerate triangles in a surface mesh
    - watertightness query: is a surface mesh a watertight manifold?

  - Surface mesh point queries

    - in/out query: is a point inside or outside a surface mesh?
    - signed distance query: find the minimum distance from a query point
      to a surface mesh

  - Point in cell query: for a query point, find the cell of the mesh
    that holds the point and the point's isoparametric coordinates within
    that cell
  - All nearest neighbors: given a list of point locations and regions,
    find all neighbors of each point in a different region

Quest also provides the ``OctreeBase`` class and its specializations
``SpatialOctree`` and ``InOutOctree``.  These classes support the queries
listed above and are also available for use by client codes.

.. toctree::
   :maxdepth: 2

   read_mesh
   check_and_repair
   point_mesh_query
   point_mesh_query_cpp
   point_in_cell
   all_nearest_neighbors

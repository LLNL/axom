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

Quest User Documentation
========================

The Quest component of Axom provides functions to read a triangle surface mesh
from a file (2D triangles comprising a manifold embedded in 3D) into a 
``mint::Mesh`` data structure.  

In addition to reading in a mesh, Quest provides several spatial queries and
operations:

  - vertex welding: merge vertices closer than a specified distance
    "epsilon"
  - find self-intersections and degenerate triangles in a surface mesh
  - watertightness query: is a surface mesh a watertight manifold?
  - in/out query: is a point inside a surface mesh, on the mesh, or outside?
  - distance query: find the distance from a query point to a surface mesh
  - point in cell query: for a query point, find the cell of the mesh
    that holds the point and the point's isoparametric coordinates within
    that cell
  - all nearest neighbors: given a list of point locations and regions,
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

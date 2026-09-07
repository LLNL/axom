.. ## Copyright (c) Lawrence Livermore National Security, LLC and other
.. ## Axom Project Contributors. See top-level LICENSE and COPYRIGHT
.. ## files for dates and other details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)


=================
Detailed examples
=================

The examples in ``src/axom/slam/examples`` show how sets, relations and maps
fit into mesh code. Start with ``UserDocs.cpp``, the quadrilateral mesh in
the :doc:`introductory example <first_example>`, 
then choose an example that matches your application.

Positions and application handles
=================================

``HandleMesh.cpp`` separates container positions from application-level
identifiers. Its node and zone sets use different integer position types and
return distinct handle types. The relation stores node positions, while maps
attach temperatures to entities and interpolation weights to individual
zone-node connections. This is a useful example when an entity identifier
cannot serve as an array index.

Unstructured hexahedral mesh
============================

``UnstructMeshField.cpp`` loads a hexahedral mesh from a VTK file and constructs
both zone-to-node and node-to-zone relations. The former has eight nodes per
zone. The latter has variable cardinality. The example traverses the mesh and
computes fields using both relations.

The executable is ``slam_unstructMesh_ex``. Its test uses mesh files from the
``slam`` directory of ``AXOM_DATA_DIR``.

Hydrodynamics examples
======================

``tinyHydro/`` uses Slam for a two-dimensional polygonal mesh and its fields.
``PolygonMeshXY.hpp`` defines the mesh accessors, and ``TinyHydroTypes.hpp``
collects the set, relation and field types. The accompanying tests exercise
the mesh and hydrodynamics routines.

``lulesh2.0.3/`` contains a Slam version of the LULESH proxy application.
``lulesh2.0.3_orig/`` contains the original version for comparison. These
examples are disabled in Windows builds because they use Unix APIs.

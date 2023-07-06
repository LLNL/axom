.. ## Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _isosurface-detection:

********************
Isosurface Detection
********************

Quest can generate isosurface meshes for node-centered scalar fields.
This feature takes a structured mesh with some scalar nodal field and
generates an ``UnstructureMesh`` at a user-specified isovalue.  The
isosurface mesh contains information on which elements of the field
mesh it crosses.  The output may be useful for material surface
reconstruction and visualization, among other things.

We support 2D and 3D configurations.  The isosurface mesh is a line
or a surface, respectively.

.. Note::

   Currently, only the 1987 Lorensen and Cline marching cubes
   algorithm is implemented.  Other similar or improved algorithms
   could be added in the future.

.. Note::

   If an input mesh cell contains an isosurface saddle point, the
   isocontour topology is ambiguous.  This implementation will choose
   the topology arbitrarily but consistently.

This feature is implemented in the class ``quest::MarchingCubes``.

The inputs are:

#. The input mesh containing the scalar field.  This mesh should be in
   Conduit's blueprint format.
   See https://llnl-conduit.readthedocs.io/en/latest/blueprint_mesh.html
#. The name of the blueprint coordinates data and scalar field data.
#. The contour value.

``MarchingCubes`` generates the isocontour mesh in an internal format.
Use ``populateContourMesh`` to put it in a ``mint::UnstructuredMesh``
object.  In the future, we will support outputs in blueprint format.

The following example shows usage of the ``MarchingCubes`` class.
(A complete example is provided in
``src/axom/quest/examples/quest_marching_cubes_example.cpp``.)

Relevant header files

.. sourcecode:: C++

   #include "conduit_relay_io_blueprint.hpp"
   #include "axom/quest/MarchingCubes.hpp"
   #include "axom/mint/mesh/UnstructuredMesh.hpp"

Set up the user's blueprint mesh and the ``MarchingCubes`` object.

.. sourcecode:: C++

   // Set up a blueprint mesh with coordset name "coordset"
   // and node-centered scalar field "scalarFieldName".
   // The blueprint mesh must be a structured mesh in
   // multi-domain format.  For single-domain format,
   // see the similar ``MarchingCubesSingleDomain`` class
   // in the same namespace.
   conduit::Node blueprintMesh = blueprint_mesh_from_user();

   quest::MarchingCubes mc(quest::MarchingCubesRuntimePolicy::seq,
                           blueprintMesh,
                           "coordset",
                           "scalarFieldName");

Run the algorithm.

.. sourcecode:: C++

   double contourValue = 0.5;
   mc.computeIsocontour(contourValue);

Place the isocontour in an output mesh.

.. sourcecode:: C++

   // Place output contour in a mint::UnstructuredMesh object.  For each
   // cell in the contourMesh, put the coresponding cell id from
   // blueprintMesh in field "cellIds" and the coreesponding domain
   // ids in field "domainIds".
   mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> contourMesh;
   mc.populateContourMesh(contourMesh, "cellIds", "domainIds");

After putting the isosurface in the ``UnstructuredMesh`` object,
the ``MarchingCubes`` object may be deleted.


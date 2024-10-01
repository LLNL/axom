.. ## Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

******************************************************
MIR Algorithms
******************************************************

The MIR component contains MIR algorithms that will take a Blueprint mesh as input,
perform MIR on it, and output a new Blueprint mesh with the reconstructed output.
A Blueprint mesh is contained in a ``conduit::Node`` and it follows the [https://llnl-conduit.readthedocs.io/en/latest/blueprint_mesh.html](Blueprint protocol),
which means the node contains specific items that describe the mesh coordinates, topology, fields, and materials.

#######
Inputs
#######

MIR algorithms are designed to accept a Conduit node containing various options that can
alter how the algorithm operates. the MIR algorithm copies the options node to the memory space
where it will run.

+---------------------------------+------------------------------------------------------+
| Option                          | Description                                          |
+=================================+======================================================+
| matset: name                    | A required string argument that specifies the name   |
|                                 | of the matset that will be operated on.              |
+---------------------------------+------------------------------------------------------+
| matsetName: name                | An optional string argument that specifies the name  |
|                                 | of the matset to create in the output. If the name   |
|                                 | is not given, the output matset will have the same   |
|                                 | name as the input matset.                            |
+---------------------------------+------------------------------------------------------+
| originalElementsField: name     | The name of the field in which to store the original |
|                                 | elements map.                                        |
+---------------------------------+------------------------------------------------------+
| selectedZones: [zone list]      | An optional argument that provides a list of zone ids|
|                                 | to operate on. The output mesh will only have        |
|                                 | contributions from zone numbers in this list, if it  |
|                                 | is given.                                            |
+---------------------------------+------------------------------------------------------+

###############
EquiZAlgorithm
###############

The [https://www.google.com/url?sa=t&source=web&rct=j&opi=89978449&url=https://www.osti.gov/servlets/purl/15014510&ved=2ahUKEwittMui-euIAxUzxOYEHXTWA2kQFnoECBcQAQ&usg=AOvVaw3qbX9qgwCn4qDP0iZ3Sq0J](Equi-Z algorithm) by J. Meredith 
is a useful visualization-oriented algorithm for MIR. Whereas many MIR algorithms 
produce disjointed element output, Equi-Z creates output that mostly forms continuous
surfaces and shapes. Continuity is achieved by averaging material volume fractions to
the mesh nodes for each material and then performing successive clipping for each
material, using the node-averaged volume fractions to determine where clipping occurs
along each edge. The basic algorithm is lookup-based so shape decomposition for a
clipped-zone can be easily determined. The cliping stage produces numerous zone fragments
that are marked with the appropriate material number and moved onto the next material
clipping stage. This concludes when all zones are comprised of only 1 material. From,
there points are made unique and the output mesh is created with a new coordset, topology,
fields, and matset.

Axom's implementation of Equi-Z can run on the CPU and the GPU. First, the zones of
interest are identified and they are classified as clean or mixed. The clean zones are
pulled out early into a new mesh and mixed zones are sent into the Equi-Z algorithm to
reconstruct clean zones. The two meshes are then merged together at the end in an output
Conduit node.

Axom's implementation supports 2D/3D zones from structured or unstructured topologies
made of Finite Element Zoo elements (e.g. triangles, quadrilaterals, tetrahedra, pyramids,
wedges, hexahedra, or topologically-compatible mixtures). The MIR logic for Equi-Z is
encapsulated in ``EquizAlgorithm``, which is a class that is templated on view objects.
View objects help provide an interface between the Blueprint data and the MIR algorithm.
At a minimum, an execution space and three views are required to instantiate the
``EquiZAlgorithm`` class. The execution space determines which compute backend will
be used to execute the algorithm. The Blueprint data must exist in a compatible
memory space for the execution space. The views are: _CoordsetView_, _TopologyView_, and
_MaterialView_. The CoordsetView template argument lets the algorithm access the mesh's
coordset using concrete data types and supports queries that return points. The
TopologyView provides a set of operations that can be performed on meshes, mainly a
device-aware method for retrieving individual zones that can be used in device kernels.
The MaterialView provides an interface for matsets.

Once view types have been created and views have been instantiated, the ``EquiZAlgorithm``
algorithm can be instantiated and used. The EquiZAlgorithm class provides a single
``execute()`` method that takes the input mesh, an options node, and a node to contain
the output mesh. The output mesh will exist in the same memory space as the input mesh,
which again, must be compatible with the selected execution space. The ``axom::mir::utilities::blueprint::copy()``
function can be used to copy Conduit nodes from one memory space to another.

.. literalinclude:: ../../examples/mir_concentric_circles.cpp
   :start-after: _equiz_mir_start
   :end-before: _equiz_mir_end
   :language: C++

The MIR output will contain a new field called "originalElements" that indicates which
original zone number gave rise to the reconstructed zone. This field makes it possible
to map back to the original mesh. The name of the field can be changed using options.

#####################
Example Applications
#####################

The mir_concentric_circles application generates a uniform mesh populated with circular mixed material shells and then it performs MIR on the input mesh before writing the reconstructed mesh.

+--------------------+---------------------------------------------------------------+
| Argument           | Description                                                   |
+====================+===============================================================+
| --gridsize number  | The number of zones along an axis.                            |
+--------------------+---------------------------------------------------------------+
| --numcircles number| The number of number of circles to use for material creation. |
+--------------------+---------------------------------------------------------------+
| --output filepath  | The file path for output files.                               |
+--------------------+---------------------------------------------------------------+
| --policy policy    | Set the execution policy (seq, omp, cuda, hip)                |
+--------------------+---------------------------------------------------------------+
| --caliper mode     | The caliper mode (none, report)                               |
+--------------------+---------------------------------------------------------------+

To run the example program from the Axom build directory, follow these steps:

  ./examples/mir_concentric_circles --gridsize 100 --numcircles 5 --output mir

#####################
Visualization
#####################

The [https://visit-dav.github.io/visit-website/](VisIt software) can be used to view the Blueprint output from MIR algorithms.
Blueprint data is saved in an HDF5 format and the top level file has a ".root" extension. Open the ".root" file in VisIt to
get started and then add a "FilledBoundary" plot of the material defined on the mesh topology. Plotting the mesh lines will
reveal that there is a single material per zone. If the input mesh is visualized in a similar manner, it will be evident that
there are multiple materials in some of the zones, if viewing a mixed material dataset.

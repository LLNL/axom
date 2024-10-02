.. ## Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

******************************************************
MIR Blueprint Utilities
******************************************************

The MIR component contains several useful building blocks for writing algorithms
for Blueprint meshes.

 * Structured as classes with an ``execute()`` method to permit them to divide algorithm into various methods.
 * Often templated on execution space and views

#######################
Copying Blueprint Data
#######################

If a ``conduit::Node`` containing Blueprint data is not on the desired memory space, it can be
moved using the ``axom::mir::utilities::blueprint::copy<ExecSpace>()`` function. The ``copy``
function moves the source ``conduit::Node`` to the destination ``conduit::Node``, making sure
to use the appropriate Axom allocator for non-string bulk arrays (e.g. arrays of ints, floats,
doubles, etc.). Data small enough to fit in a ``conduit::Node`` and strings are left in the host
memory space, which lets algorithms on the host side query them. For data that have been moved
to the device, their sizes and data types can still be queried using normal Conduit mechanisms
such as getting metadata via the ``conduit::Node::dtype()`` method.

  .. codeblock:: cpp

      conduit::Node hostMesh, deviceMesh, hostMesh2;
      // host->device
      axom::mir::utilities::blueprint::copy<axom::HIP_EXEC<256>>(deviceMesh, hostMesh);
      // device->host
      axom::mir::utilities::blueprint::copy<axom::SEQ_EXEC>(hostMesh2, deviceMesh);

############################
ConduitAllocateThroughAxom
############################

When writing algorithms that construct Blueprint data, it is helpful to force Conduit
to allocate its memory through Axom's allocation routines and then make an ``axom::ArrayView``
of the data in the Conduit node. This prevents data from having to be copied from an Axom
data structure into a Conduit node since it can be constructed from the start inside the
Conduit node. The size of the array must be known.

The ``axom::mir::utilities::blueprint::ConduitAllocateThroughAxom``
class is a template class that takes an execution space as a template argument and it
installs an allocation routine in Conduit that can be used to allocate data through
Axom. The Conduit allocator is set on each ``conduit::Node`` before setting data into
the object.

.. literalinclude:: ../../ClipField.cpp
   :start-after: _mir_utilities_c2a_begin
   :end-before: _mir_utilities_c2a_end
   :language: C++

##########
ClipField
##########

The ``axom::mir::clipping::ClipField`` class intersects all the zones in the input Blueprint
mesh with an implicit surface where the selected input field equals zero and produces a new
Blueprint mesh based on the selected zone fragments produced by the intersection. This can be thought
of as an isosurface algorithm but with a volumetric output mesh where the mesh is either inside or
outside of the selected isovalue. The ``ClipField`` class has multiple template arguments to
select the execution space, the type of topology view, the type of coordset view, and the
type of intersector used to determine intersections. The default intersection uses an isosurface-
based intersection method, though other intersectors could be created to perform plane
or sphere intersections.

.. literalinclude:: ../../ClipField.cpp
   :start-after: _mir_utilities_clipfield_start
   :end-before: _mir_utilities_clipfield_end
   :language: C++

################
CoordsetBlender
################

The ``axom::mir::utilities::blueprint::CoordsetBlender`` class takes a "BlendGroup" and makes
a new explicit coordset where each new point corresponds to one blend group. A "BlendGroup" is
an object that groups several array views that describe a set of blend groups. Each blend group
is formed from a list of node ids and weight values and the new coordinate is formed by looking
up the points in the blend group in the source coordset and multiplying them by their weights
and summing them together to produce the new point in the output coordset. Classes such as
``ClipField`` use ``CoordsetBlender`` to make new coordsets that contain points that were a
combination of multiple points in the input coordset.

.. literalinclude:: ../../ClipField.cpp
   :start-after: _mir_utilities_coordsetblender_start
   :end-before: _mir_utilities_coordsetblender_end
   :language: C++

################
CoordsetSlicer
################

The ``axom::mir::utilities::blueprint::CoordsetSlicer`` class takes ``SliceData`` and makes a
new explicit coordset where each point corresponds to a single index from the node indices
stored in SliceData. This class can be used to select a subset of a coordset, reorder nodes
in a coordset, or repeat nodes in a coordset.

.. literalinclude:: ../../ExtractZones.cpp
   :start-after: _mir_utilities_coordsetslicer_begin
   :end-before: _mir_utilities_coordsetslicer_end
   :language: C++

##################
ExtractZones
##################

The ``axom::mir::utilities::ExtractZones`` class takes a list of selected zone ids and extracts
a new mesh from a source mesh that includes only the selected zones. There is a derived class
``ExtractZonesAndMatset`` that also extracts a matset, if present.

.. literalinclude:: ../../ExtractZones.cpp
   :start-after: _mir_utilities_extractzones_begin
   :end-before: _mir_utilities_extractzones_end
   :language: C++

#############
FieldBlender
#############

The ``axom::mir::utilities::blueprint::FieldBlender`` class is similar to the ``CoordsetBlender``
class, except that it operates on a field instead of coordsets. The class is used to create a
new field that includes values derived from multiple weighted source values.

############
FieldSlicer
############

The ``axom::mir::utilities::blueprint::FieldSlicer`` class selects specific indices from a
field and makes a new field.

.. literalinclude:: ../../tests/mir_blueprint_utilities.cpp
   :start-after: _mir_utilities_fieldslicer_begin
   :end-before: _mir_utilities_fieldslicer_end
   :language: C++

##################
MakeUnstructured
##################

The ``axom::mir::utilities::blueprint::MakeUnstructured`` class takes a structured topology
and creates a new unstructured topology. This class does not need views to wrap the input
structured topology.

.. literalinclude:: ../../tests/mir_blueprint_utilities.cpp
   :start-after: _mir_utilities_makeunstructured_begin
   :end-before: _mir_utilities_makeunstructured_end
   :language: C++

##################
MatsetSlicer
##################

The ``axom::mir::utilities::blueprint::MatsetSlicer`` class is similar to the ``FieldSlicer``
class except it slices matsets instead of fields. The same ``SliceData`` can be passed to
MatsetSlicer to pull out and assemble a new matset data for a specific list of zones.

.. literalinclude:: ../../ExtractZones.cpp
   :start-after: _mir_utilities_matsetslicer_begin
   :end-before: _mir_utilities_matsetslicer_end
   :language: C++

##################
MergeMeshes
##################

The ``axom::mir::utilities::blueprint::MergeMeshes`` class merges data for coordsets,
topology, and fields from multiple input meshes into a new combined mesh. The class also
supports renaming nodes using a map that converts a local mesh's node ids to the final
output node numbering, enabling meshes to be merged such that some nodes get combined.
A derived class can also merge matsets.

.. literalinclude:: ../../EquiZAlgorithm.hpp
   :start-after: _mir_utilities_mergemeshes_begin
   :end-before: _mir_utilities_mergemeshes_end
   :language: C++

###########################
NodeToZoneRelationBuilder
###########################

The ``axom::mir::utilities::blueprint::NodeToZoneRelationBuilder`` class creates a Blueprint
O2M (one to many) relation that relates node numbers to the zones that contain them. This mapping
is akin to inverting the normal mesh connectivity which is a map of zones to node ids. The O2M
relation is useful for recentering data from the zones to the nodes.

.. literalinclude:: ../../tests/mir_blueprint_utilities.hpp
   :start-after: _mir_utilities_n2zrel_begin
   :end-before: _mir_utilities_n2zrel_end
   :language: C++

###############
RecenterFields
############### 

The ``axom::mir::utilities::blueprint::RecenterFields`` class uses an O2M relation to average
field data from multiple values to an averaged value. In Axom, this is used to convert a field
associated with the elements to a new field associated with the nodes.

.. literalinclude:: ../../tests/mir_blueprint_utilities.hpp
   :start-after: _mir_utilities_recenterfield_begin
   :end-before: _mir_utilities_recenterfield_end
   :language: C++

##################
Unique
##################

The ``axom::mir::utilities::Unique`` class can take an unsorted list of values and produce a
sorted list of unique outputs, along with a list of offsets into the original values to identify
one representative value in the original list for each unique value. This class is used to help
merge points.

.. literalinclude:: ../../tests/mir_clipfield.hpp
   :start-after: _mir_utilities_unique_begin
   :end-before: _mir_utilities_unique_end
   :language: C++

##################
ZoneListBuilder
##################

The ``axom::mir::utilities::blueprint::ZoneListBuilder`` class takes a matset view and a list
of selected zone ids and makes two output lists of zone ids that correspond to clean zones and
mixed zones (more than 1 material in the zone).

.. literalinclude:: ../../EquiZAlgorithm.hpp
   :start-after: _mir_utilities_zlb_begin
   :end-before: _mir_utilities_zlb_end
   :language: C++

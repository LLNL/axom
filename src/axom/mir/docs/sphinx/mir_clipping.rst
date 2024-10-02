.. ## Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

*************
Clipping
*************

The MIR component provides a clipping algorithm that can perform isosurface-based
clipping and return volumetric output for zones and partial zones that are "inside"
or "outside" the clip boundary. The clipping algorithm is implemented in the
``axom::mir::clipping::ClipField`` class. The class can be instantiated with several
template arguments that govern where it will execute, which coordset and topology
types it supports, and how it performs intersection. The input to the algorithm is
a Blueprint mesh and when instantiated with coordset and topology views appropriate
for the input data, the algorithm can operate on a wide variety of mesh types. This
includes 2D/3D structured and unstructured topologies, that can be represented using
finite elements.

By default, the algorithm will clip using a field but other intersection routines
can be substituted via a template argument to facilitate creation of clipping using
planes, spheres, surfaces of revolution, etc. The Equi-Z algorithm uses ClipField
with an intersector that uses material volume fractions to determine the clipped geometry.

#######
Inputs
#######

Like the MIR algorithms, the clipping algorithm is designed to accept a Conduit node
containing various options that can influence how the algorithm operates. The clipping
algorithm copies the options node to the memory space where it will be used.

+---------------------------------+------------------------------------------------------+
| Option                          | Description                                          |
+=================================+======================================================+
| clipField: name                 | A required string argument that specifies the name   |
|                                 | of the field that is used for clipping. At present,  |
|                                 | the field must be a vertex-associated field.         |
+---------------------------------+------------------------------------------------------+
| clipValue: value                | An optional numeric argument that specifies the      |
|                                 | value in the field at which the clip boundary is     |
|                                 | defined. The default is 0.                           |
+---------------------------------+------------------------------------------------------+
| colorField: name                | If inside=1 and outside=1 then a color field is      |
|                                 | generated so it is possible to tell apart regions of |
|                                 | the clip output that were inside or outside the clip |
|                                 | boundary. This field permits the user to change the  |
|                                 | name of the color field, which is called "color" by  |
|                                 | default.                                             |
+---------------------------------+------------------------------------------------------+
| inside: number                  | Indicates to the clipping algorithm that it should   |
|                                 | preserve zone fragments that were "inside" the clip  |
|                                 | boundary. Set to 1 to enable, 0 to disable. The      |
|                                 | algorithm will generate these fragments by default.  |
+---------------------------------+------------------------------------------------------+
| originalElementsField: name     | The name of the field in which to store the original |
|                                 | elements map.                                        |
+---------------------------------+------------------------------------------------------+
| outside: number                 | Indicates to the clipping algorithm that it should   |
|                                 | preserve zone fragments "outside" the clip boundary. |
|                                 | Set to 1 to enable, 0 to disable. These fragments are|
|                                 | not on by default.                                   |
+---------------------------------+------------------------------------------------------+
| selectedZones: [zone list]      | An optional argument that provides a list of zone ids|
|                                 | on which to operate. The output mesh will only have  |
|                                 | contributions from zone numbers in this list, if it  |
|                                 | is given.                                            |
+---------------------------------+------------------------------------------------------+

##########
ClipField
##########

To use the ``ClipField`` class, one must have Blueprint data with at least one vertex-associated
field. Views for the coordset and topology are created and used to instantiate the ``ClipField``
class, which then takes an input node for the Blueprint mesh, an input node that contains
the options, and a 3rd output node that will contain the clip output. The input mesh node
needs to contain data arrays for coordinates, mesh topology, and fields that are in the
memory space of the targeted device. Other Conduit nodes that contain strings or single numbers
that can fit within a node are safe remaining in host memory. If the mesh is not in the
desired memory space, it can be moved using ``axom::mir::utilities::blueprint::copy()``.

.. literalinclude:: ../../ClipField.cpp
   :start-after: _mir_utilities_clipfield_start
   :end-before: _mir_utilities_clipfield_end
   :language: C++

#############
Intersectors
#############

An intersector is a class that is passed as a template argument to ``ClipField``. The intersector
determines how the ``ClipField`` algorithm will generate intersection cases, for each zone
in the mesh. The ``ClipField`` algorithm default intersector uses a field to determine clip
cases, resulting in isosurface behavior for the geometry intersections. Alternative intersectors
can be provided to achieve other types of intersections.

An intersector needs to provide an interface like the following:

 .. codeblock{.cpp}::

    template <typename ExecSpace, typename ConnectivityType>
    class CustomIntersector
    {
    public:
      using ConnectivityView = axom::ArrayView<ConnectivityType>;

      // Internal view - runs on device.
      struct View
      {
        // Given a zone index and the node ids that comprise the zone, return
        // the appropriate clip case.
        AXOM_HOST_DEVICE
        axom::IndexType determineClipCase(axom::IndexType zoneIndex,
                                          const ConnectivityView &nodeIds) const
        {
          axom::IndexType clipcase = 0;
          for(IndexType i = 0; i < nodeIds.size(); i++)
          {
            const auto id = nodeIds[i];
            const auto value = // Compute distance from node to surface.
            clipcase |= (value > 0) ? (1 << i) : 0;
          }
          return clipcase;
        }

        // Compute the weight[0,1] of a clip value along an edge (id0, id1) using the clip field and value.
        AXOM_HOST_DEVICE
        ClipFieldType computeWeight(axom::IndexType zoneIndex,
                                    ConnectivityType id0,
                                    ConnectivityType id1) const
        {
          return 1.;
        }
      };

      // Initialize the object from options (on host).
      void initialize(const conduit::Node &n_options, const conduit::Node &n_fields)
      {
      }

      // Determine the name of the topology on which to operate.
      std::string getTopologyName(const conduit::Node &n_input,
                                  const conduit::Node &n_options) const
      {
        return "mesh";
      }

      // Return a new instance of the view.
      View view() const { return m_view; }

      View m_view;
    };

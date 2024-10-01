.. ## Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level COPYRIGHT file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

******
Views
******

----------
ArrayView
----------

Axom provides ``axom::ArrayView`` to wrap data in a non-owning data structure that can be passed to
kernels. The MIR component provides functions that help wrap arrays stored in 



----------
Coordsets
----------

----------------
Topology Views
----------------

Topology views provide a layer on top of the Blueprint mesh topology that enables Axom algorithms
to be written while not needing to care specifically about topology types and data types.
Axom provides topology views for structured meshes and unstructured meshes. 

^^^^^^^^^^^^^^^^^^^^^^
Structured Mesh Views
^^^^^^^^^^^^^^^^^^^^^^

The structured mesh topology view, ``StructuredTopologyView``, pertains to any of the Blueprint
topology types. The ``StructuredTopologyView`` class is a template that takes an indexing policy
as a template argument. The indexing policy computes zone indices and converts to/from
logical/global indices. The ``StridedStructuredIndexingPolicy`` class supports indexing for
strided-structured Blueprint meshes, which are structured meshes that exist over a sub-window
of the overall mesh. There are helper functions for creating structured topology views from
a Conduit node.

  .. codeblock{.cpp}::

     conduit::Node &n_topo1 = n_mesh["topologies/mesh2d"];
     conduit::Node &n_topo2 = n_mesh["topologies/mesh3d"];
     conduit::Node &n_topo3 = n_mesh["topologies/mesh2dss"];
     // Make a 2D structured mesh view from the topology.
     auto topologyView1 = axom::mir::views::make_structured<2>::view(n_topo1);
     // Make a 3D structured mesh view from the topology.
     auto topologyView2 = axom::mir::views::make_structured<2>::view(n_topo2);
     // Make a 2D strided-structured mesh view from the topology.
     auto topologyView3 = axom::mir::views::make_strided_structured<2>::view(n_topo3);

^^^^^^^^^^^^^^^^^^^^^^^^
Unstructured Mesh Views
^^^^^^^^^^^^^^^^^^^^^^^^

There are 3 unstructured mesh views. The ``UnstructuredTopologySingleShapeView`` class wraps a
Blueprint topology that contains a single zone/shape type. The zone type is a template argument
that determines the type of zone that is held within the topology.

  .. codeblock{.cpp}::

     // Make a topology view for a tetrahedral mesh with int connectivity.
     namespace bputils = axom::mir::utilities::blueprint;
     const conduit::Node &n_topo = n_mesh["topologies/mesh"];
     const auto connView = bputils::make_array_view<int>(n_topo["elements/connectivity"]);
     axom::mir::views::UnstructuredTopologySingleShapeView<axom::mir::views::TetShape<int>> view(connView);

There are multiple shape types defined in ``axom/mir/views/Shapes.hpp`` that can be used with
the ``UnstructuredTopologySingleShapeView``: TriShape, QuadShape, TetShape, PyramidShape,
WedgeShape, and HexShape.

Blueprint supports "mixed" topologies that contain multiple shape types. These topologies are
handled using the ``axom::mir::views::UnstructuredTopologyMixedShapeView``. Additional array
views are needed to supply the sizes, offsets, and shapes arrays.

  .. codeblock{.cpp}::

     // A shape map helps map values from the values used in the Blueprint topology to
     // the shape ids used in Axom.
     const conduit::Node &n_topo = n_mesh["topologies/mesh"];
     axom::Array<int> ids, values;
     auto shapeMap = axom::mir::views::buildShapeMap(n_topo);

     namespace bputils = axom::mir::utilities::blueprint;
     axom::mir::views::UnstructuredTopologyMixedShapeView<int> view(
       bputils::make_array_view(n_topo["elements/connectivity"),
       bputils::make_array_view(n_topo["elements/sizes"),
       bputils::make_array_view(n_topo["elements/offsets"),
       bputils::make_array_view(n_topo["elements/shapes"),
       shapeMap);

  .. codeblock{.cpp}::

     topologyView = ...
     axom::for_all<ExecSpace>(topologyView.numberOfZones(), AXOM_LAMBDA(axom::IndexType zoneIndex)
     {
       // Get the current zone.
       const auto zone = topologyView.zone(zoneIndex);

       // Iterate over this zone's nodes.
       for(const auto &nodeId : zone.getIds())
       {
         // Do something.
       }
     });

----------
Matsets
----------

----------
Dispatch
----------

There are several helper functions for converting a Conduit node to a specific view type
and passing the view to a lambda for further processing. These dispatch functions take
care of wrapping a Conduit node containing Blueprint data in various view types before
passing the views to a user-supplied lambda function. The lambda function will typically
be instantiated multiple times to handle cases when there are multiple data types and
object types (e.g. coordsets). Generic lambdas can be used to process multiple view
and it is possible to nest multiple dispatch functions.

 .. codeblock{.cpp}::

    const conduit::Node &n_coordset = n_mesh["coordsets/coords"];
    axom::mir::views::dispatch_coordset(n_coordset, [&](auto coordsetView) {
      // Get the C++ type of the coordset.
      using CoordsetView = decltype(CoordsetView);

      // Implement algorithm using coordsetView.
    });

Dispatch functions for topologies enable creation of algorithms that can operate on multiple
topology types through a topology view. These dispatch functions can be called for specific
topology types such as unstructured topologies or they can be called to implement algorithms
that can operate on any topology.

 .. codeblock{.cpp}::

    const conduit::Node &n_topo = n_mesh["topologies/mesh"];
    // Handle rectilinear topology type.
    axom::mir::views::dispatch_rectilinear_topology(n_topo, [&](auto topologyView) {
    });
    // Handle structured topology types
    axom::mir::views::dispatch_structured_topology(n_topo, [&](auto topologyView) {
    });
    // Handle unstructured topology types
    axom::mir::views::dispatch_unstructured_topology(n_topo, [&](auto topologyView) {
    });
    // Handle any topology type.
    axom::mir::views::dispatch_topologies(n_topo, [&](auto topologyView) {
    });

Nesting dispatch functions permits the calling code to handle both coordset views and
topology views using a single lambda function for the algorithm. For portability, the
algorithm should be placed in a function or class member method when instantiated from
the anonymous lambda function from the dispatch functions.

  .. codeblock{.cpp}::

    struct Algorithm
    {
      void execute(const conduit::Node &n_mesh)
      {
        // Handle product of coordset types and topology types.
        axom::mir::views::dispatch_coordset(n_mesh["coordsets/coords"], [&](auto coordsetView)
        {
          axom::mir::views::dispatch_topologies(n_mesh["topologies/mesh"], [&](auto topologyView)
          {
            implementation(coordsetView, topologyView);
          });
        });
      }

      template <typename CoordsetView, typename TopologyView>
      void implementation(CoordsetView coordsetView, TopologyView topologyView) const
      {
        // Do algorithm that involves coordsetView and topologyView.
      }
    };


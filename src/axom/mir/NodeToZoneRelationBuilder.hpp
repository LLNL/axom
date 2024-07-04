// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_NODE_TO_ZONE_RELATION_BUILDER_HPP_
#define AXOM_MIR_NODE_TO_ZONE_RELATION_BUILDER_HPP_

#include "axom/core/memory_management.hpp"
#include "axom/core/execution/execution_space.hpp"
#include "axom/core/execution/for_all.hpp"
#include "axom/core/Array.hpp"
#include "axom/core/ArrayView.hpp"
#include "axom/mir/utilities.hpp"
#include "axom/mir/blueprint_utilities.hpp"
#include "axom/mir/views/dispatch_unstructured_topology.hpp"

#include <conduit/conduit.hpp>
#include <conduit/conduit_blueprint.hpp>
#include <conduit/conduit_blueprint_mesh_utils.hpp>

#include <RAJA/RAJA.hpp>

namespace axom
{
namespace mir
{
namespace utilities
{

/**
 * \brief Build an o2m relation that lets us look up the zones for a node.
 */
template <typename ExecSpace>
class NodeToZoneRelationBuilder
{
public:

  /**
   * \brief Build a node to zone relation and store the resulting O2M relation in the \a relation conduit node.
   *
   * \param topo The topology for which we're building the O2M relation.
   * \param[out] The node that will contain the O2M relation.
   */
  void execute(const conduit::Node &topo, conduit::Node &relation);

private:
  /**
   * \brief Given views that contain the nodes and zones, sort the zones using the
   *        node numbers to produce a list of zones for each node and an offsets array
   *        that points to the start of each list of zones.
   * 
   * \param[in]    nodes_view   A view that contains the set of all of the nodes in the topology (the connectivity)
   * \param[inout[ zones_view   A view (same size as \a nodes_view) that contains the zone number of each node.
   * \param[out]   offsets_view A view that we fill with offsets so offsets_view[i] points to the start of the i'th list in \a zones_view.
   */
  template <typename ViewType>
  void buildRelation(const ViewType &nodes_view, ViewType &zones_view, ViewType &offsets_view) const;
};


template <typename ExecSpace>
template <typename ViewType>
void
NodeToZoneRelationBuilder<ExecSpace>::buildRelation(const ViewType &nodes_view, ViewType &zones_view, ViewType &offsets_view) const
{
  assert(nodes_view.size() == zones_view.size());

  using loop_policy = typename axom::execution_space<ExecSpace>::loop_policy;
  const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();

  // Make a copy of the nodes that we'll use as keys.
  const auto n = nodes_view.size();
  axom::Array<axom::IndexType> keys(n, n, allocatorID);
  auto keys_view = keys.view();
  axom::for_all<ExecSpace>(n, AXOM_LAMBDA(axom::IndexType i)
  {
    keys_view[i] = nodes_view[i];
  });

  // Sort the keys, zones in place. This sorts the zones_view which we want for output.
  RAJA::sort_pairs<loop_policy>(RAJA::make_span(keys_view, n),
                                RAJA::make_span(zones_view, n));

  // Make a mask array for where differences occur.
  axom::Array<axom::IndexType> mask(n, n, allocatorID);
  auto mask_view = mask.view();
  axom::for_all<ExecSpace>(n, AXOM_LAMBDA(axom::IndexType i)
  {
    const axom::IndexType different = (keys_view[i] != keys_view[i - 1]) ? 1 : 0;
    const axom::IndexType m = (i >= 1) ? different : 1;
    mask_view[i] = m;
  });

  // Do a scan on the mask array to build an offset array.
  axom::Array<axom::IndexType> dest_offsets(n, n, allocatorID);
  auto dest_offsets_view = dest_offsets.view();
  RAJA::exclusive_scan<loop_policy>(RAJA::make_span(mask_view, n),
                                    RAJA::make_span(dest_offsets_view, n),
                                    RAJA::operators::plus<axom::IndexType>{});

  // Build the offsets to each node's zone ids.
  axom::for_all<ExecSpace>(offsets_view.size(), AXOM_LAMBDA(axom::IndexType i)
  {
    offsets_view[i] = 0;
  });
  axom::for_all<ExecSpace>(n, AXOM_LAMBDA(axom::IndexType i)
  {
    if(mask_view[i])
    {
      offsets_view[dest_offsets_view[i]] = i;
    }
  });
}

template <typename ExecSpace>
void
NodeToZoneRelationBuilder<ExecSpace>::execute(const conduit::Node &topo, conduit::Node &relation)
{
  using loop_policy = typename axom::execution_space<ExecSpace>::loop_policy;
  using reduce_policy = typename axom::execution_space<ExecSpace>::reduce_policy;
  const std::string type = topo.fetch_existing("type").as_string();
  const auto allocatorID = axom::execution_space<ExecSpace>::allocatorID();

  conduit::Node &n_zones = relation["zones"];
  conduit::Node &n_sizes = relation["sizes"];
  conduit::Node &n_offsets = relation["offsets"];
  n_zones.set_allocator(allocatorID);
  n_sizes.set_allocator(allocatorID);
  n_offsets.set_allocator(allocatorID);

  const conduit::Node *coordset = conduit::blueprint::mesh::utils::find_reference_node(topo, "coordset");

  if(type == "unstructured")
  {
    conduit::blueprint::mesh::utils::ShapeType shape(topo);
    const conduit::Node &n_connectivity = topo["elements/connectivity"];
    const auto intTypeId = n_connectivity.dtype().id();
    const auto connSize = n_connectivity.dtype().number_of_elements();

    // Use the coordset to get the number of nodes. Conduit should be able to do this using only metadata.
    const auto nnodes = conduit::blueprint::mesh::utils::coordset::length(*coordset);

    if(shape.is_polyhedral())
    {
      views::dispatch_unstructured_polyhedral_topology(topo, [&](auto topoView)
      {
        const auto nzones = topoView.numberOfZones();
        axom::Array<axom::IndexType> sizes(nzones, nzones, allocatorID);
        auto sizes_view = sizes.view();

        // Run through the topology once to do a count of each zone's unique node ids.
        RAJA::ReduceSum<reduce_policy, axom::IndexType> count(0);
        topoView.template for_all_zones<ExecSpace>(AXOM_LAMBDA(auto zoneIndex, const auto &zone)
        {          
          const auto uniqueIds = zone.getUniqueIds();
          sizes_view[zoneIndex] = uniqueIds.size();
          count += uniqueIds.size();
        });
        const auto connSize = count.get();

        // Do a scan on the size array to build an offset array.
        axom::Array<axom::IndexType> offsets(nzones, nzones, allocatorID);
        auto offsets_view = offsets.view();
        RAJA::exclusive_scan<loop_policy>(RAJA::make_span(sizes_view, nzones),
                                          RAJA::make_span(offsets_view, nzones),
                                          RAJA::operators::plus<axom::IndexType>{});
        sizes.clear();

        // Allocate Conduit arrays on the device in a data type that matches the connectivity.
        conduit::Node n_conn;
        n_conn.set_allocator(allocatorID);
        n_conn.set(conduit::DataType(intTypeId, connSize));

        n_zones.set(conduit::DataType(intTypeId, connSize));
        n_sizes.set(conduit::DataType(intTypeId, nnodes));
        n_offsets.set(conduit::DataType(intTypeId, nnodes));

        views::IndexNode_to_ArrayView_same(n_conn, n_zones, n_sizes, n_offsets, [&](auto connectivityView, auto zonesView, auto sizesView, auto offsetsView)
        {         
          // Run through the data one more time to build the nodes and zones arrays.
          topoView.template for_all_zones<ExecSpace>(AXOM_LAMBDA(auto zoneIndex, const auto &zone)
          {          
            const auto uniqueIds = zone.getUniqueIds();
            auto destIdx = offsets_view[zoneIndex];
            for(axom::IndexType i = 0; i < uniqueIds.size(); i++, destIdx++)
            {
              connectivityView[destIdx] = uniqueIds[i];
              zonesView[destIdx] = zoneIndex;
            }
          });

          // Make the relation, outputting into the zonesView and offsetsView.
          using ViewType = decltype(connectivityView);
          buildRelation<ViewType>(connectivityView, zonesView, offsetsView);

          // Compute sizes from offsets.
          axom::for_all<ExecSpace>(offsetsView.size(), AXOM_LAMBDA(auto i)
          {
            sizesView[i] = (i < offsetsView.size() - 1) ? (offsetsView[i + 1]  - offsetsView[i]) : (connSize - offsetsView[i]);
          });
        });
      });
    }
    else if(shape.is_polygonal())
    {
      const conduit::Node &n_topo_sizes = topo["elements/sizes"];
      const conduit::Node &n_topo_offsets = topo["elements/offsets"];

      const auto nzones = n_topo_sizes.dtype().number_of_elements();

      // Allocate Conduit arrays on the device in a data type that matches the connectivity.
      n_zones.set(conduit::DataType(intTypeId, connSize));
      n_sizes.set(conduit::DataType(intTypeId, nnodes));
      n_offsets.set(conduit::DataType(intTypeId, nnodes));

      // Make zones for each node
      views::IndexNode_to_ArrayView_same(n_zones, n_topo_sizes, n_topo_offsets, [&](auto zonesView, auto sizesView, auto offsetsView)
      {
        using DataType = typename decltype(zonesView)::value_type;
        axom::for_all<ExecSpace>(0, nzones, AXOM_LAMBDA(axom::IndexType zoneIndex)
        {
          for(DataType i = 0; i < sizesView[zoneIndex]; i++)
            zonesView[offsetsView[zoneIndex] + i] = zoneIndex;
        });
      });

      views::IndexNode_to_ArrayView_same(n_connectivity, n_zones, n_sizes, n_offsets, [&](auto connectivityView, auto zonesView, auto sizesView, auto offsetsView)
      {
        // Make the relation, outputting into the zonesView and offsetsView.
        using ViewType = decltype(connectivityView);
        buildRelation<ViewType>(connectivityView, zonesView, offsetsView);
        
        // Compute sizes from offsets.
        axom::for_all<ExecSpace>(offsetsView.size(), AXOM_LAMBDA(auto i)
        {
          sizesView[i] = (i < offsetsView.size() - 1) ? (offsetsView[i + 1]  - offsetsView[i]) : (connSize - offsetsView[i]);
        });
      });
    }
    else
    {
      // Shapes are all the same size.
      const auto nodesPerShape = shape.indices;
      const auto nzones = connSize / nodesPerShape;

      // Allocate Conduit arrays on the device in a data type that matches the connectivity.
      n_zones.set(conduit::DataType(intTypeId, connSize));
      n_sizes.set(conduit::DataType(intTypeId, nnodes));
      n_offsets.set(conduit::DataType(intTypeId, nnodes));

      views::IndexNode_to_ArrayView_same(n_connectivity, n_zones, n_sizes, n_offsets, [&](auto connectivityView, auto zonesView, auto sizesView, auto offsetsView)
      {
        // Make zones for each node
        axom::for_all<ExecSpace>(0, connSize, AXOM_LAMBDA(axom::IndexType index)
        {
          zonesView[index] = index / nodesPerShape;
        });
        // Make the relation, outputting into the zonesView and offsetsView.
//        using ViewType = decltype(connectivityView);
        buildRelation<ExecSpace>(connectivityView, zonesView, offsetsView);
        
        // Compute sizes from offsets.
        axom::for_all<ExecSpace>(offsetsView.size(), AXOM_LAMBDA(auto i)
        {
          sizesView[i] = (i < offsetsView.size() - 1) ? (offsetsView[i + 1]  - offsetsView[i]) : (connSize - offsetsView[i]);
        });
      });
    }
  }
  else
  {
    // These are all structured topos of some sort. Make an unstructured representation and recurse.

    conduit::Node mesh;
    axom::mir::utilities::blueprint::to_unstructured<ExecSpace>(topo, *coordset, "newtopo", mesh);
        
    // Recurse using the unstructured mesh.
    execute(mesh["newtopo"], relation);
  }
}

} // end namespace utilities
} // end namespace mir
} // end namespace axom

#endif

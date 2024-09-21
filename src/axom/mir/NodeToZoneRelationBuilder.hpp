// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_NODE_TO_ZONE_RELATION_BUILDER_HPP_
#define AXOM_MIR_NODE_TO_ZONE_RELATION_BUILDER_HPP_

#include "axom/core.hpp"
#include "axom/mir/utilities.hpp"
#include "axom/mir/blueprint_utilities.hpp"
#include "axom/mir/views/dispatch_unstructured_topology.hpp"

#include <conduit/conduit.hpp>
#include <conduit/conduit_blueprint.hpp>
#include <conduit/conduit_blueprint_mesh_utils.hpp>

#include <RAJA/RAJA.hpp>

#include <map>
#include <vector>

namespace axom
{
namespace mir
{
namespace utilities
{
namespace blueprint
{
namespace details
{
/**
   * \brief Given views that contain the nodes and zones, sort the zones using the
   *        node numbers to produce a list of zones for each node and an offsets array
   *        that points to the start of each list of zones.
   * 
   * \param[in]    nodesView   A view that contains the set of all of the nodes in the topology (the connectivity)
   * \param[inout[ zonesView   A view (same size as \a nodesView) that contains the zone number of each node.
   * \param[out]   offsetsView A view that we fill with offsets so offsetsView[i] points to the start of the i'th list in \a zonesView.
   *
   * \note RAJA::sort_pairs can be slow if there are a lot of nodes (depends on ExecSpace too).
   */
template <typename ExecSpace, typename ViewType>
struct BuildRelation
{
  static void execute(const ViewType &nodesView,
                      ViewType &zonesView,
                      ViewType &sizesView,
                      ViewType &offsetsView)
  {
    AXOM_ANNOTATE_SCOPE("FillZonesAndOffsets");
    assert(nodesView.size() == zonesView.size());

    using loop_policy = typename axom::execution_space<ExecSpace>::loop_policy;
    using value_type = typename ViewType::value_type;
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    // Make a copy of the nodes that we'll use as keys.
    const auto n = nodesView.size();
    axom::Array<value_type> keys(n, n, allocatorID);
    auto keysView = keys.view();
    axom::for_all<ExecSpace>(
      n,
      AXOM_LAMBDA(axom::IndexType i) { keysView[i] = nodesView[i]; });

    // Sort the keys, zones in place. This sorts the zonesView which we want for output.
    RAJA::sort_pairs<loop_policy>(RAJA::make_span(keysView.data(), n),
                                  RAJA::make_span(zonesView.data(), n));

    // Make a mask array for where differences occur.
    axom::Array<axom::IndexType> mask(n, n, allocatorID);
    auto maskView = mask.view();
    axom::for_all<ExecSpace>(
      n,
      AXOM_LAMBDA(axom::IndexType i) {
        maskView[i] = (i >= 1) ? ((keysView[i] != keysView[i - 1]) ? 1 : 0) : 1;
      });

    // Do a scan on the mask array to build an offset array.
    axom::Array<axom::IndexType> dest_offsets(n, n, allocatorID);
    auto dest_offsetsView = dest_offsets.view();
    axom::exclusive_scan<ExecSpace>(maskView, dest_offsetsView);

    // Build the offsets to each node's zone ids.
    axom::for_all<ExecSpace>(
      offsetsView.size(),
      AXOM_LAMBDA(axom::IndexType i) { offsetsView[i] = 0; });
    axom::for_all<ExecSpace>(
      n,
      AXOM_LAMBDA(axom::IndexType i) {
        if(maskView[i])
        {
          offsetsView[dest_offsetsView[i]] = i;
        }
      });

    // Compute sizes from offsets.
    const value_type totalSize = nodesView.size();
    axom::for_all<ExecSpace>(
      offsetsView.size(),
      AXOM_LAMBDA(axom::IndexType i) {
        sizesView[i] = (i < offsetsView.size() - 1)
          ? (offsetsView[i + 1] - offsetsView[i])
          : (totalSize - offsetsView[i]);
      });
  }
};

/// Partial specialization for axom::SEQ_EXEC.
template <typename ViewType>
struct BuildRelation<axom::SEQ_EXEC, ViewType>
{
  static void execute(const ViewType &nodesView,
                      ViewType &zonesView,
                      ViewType &sizesView,
                      ViewType &offsetsView)
  {
    AXOM_ANNOTATE_SCOPE("FillZonesAndOffsets");
    assert(nodesView.size() == zonesView.size());
    using value_type = typename ViewType::value_type;
    using ExecSpace = axom::SEQ_EXEC;
    const int allocatorID = execution_space<ExecSpace>::allocatorID();

    // Count how many times a node is used.
    const auto nnodes = offsetsView.size();
    axom::for_all<ExecSpace>(
      nnodes,
      AXOM_LAMBDA(axom::IndexType index) { sizesView[index] = 0; });
    axom::for_all<ExecSpace>(
      nodesView.size(),
      AXOM_LAMBDA(axom::IndexType index) {
        sizesView[nodesView[index]]++;  // Works only because ExecSpace=SEQ_EXEC.
      });
    // Make offsets
    axom::exclusive_scan<ExecSpace>(sizesView, offsetsView);

    axom::for_all<ExecSpace>(
      nnodes,
      AXOM_LAMBDA(axom::IndexType index) { sizesView[index] = 0; });

    axom::Array<value_type> zcopy(zonesView.size(), zonesView.size(), allocatorID);
    axom::copy(zcopy.data(),
               zonesView.data(),
               zonesView.size() * sizeof(value_type));
    auto zcopyView = zcopy.view();

    // Fill in zonesView, sizesView with each node's zones.
    axom::for_all<ExecSpace>(
      nodesView.size(),
      AXOM_LAMBDA(axom::IndexType index) {
        const auto ni = nodesView[index];
        const auto destOffset = offsetsView[ni] + sizesView[ni];
        zonesView[destOffset] = zcopyView[index];
        sizesView[ni]++;  // Works only because ExecSpace=SEQ_EXEC.
      });
  }
};

}  // end namespace details

/**
 * \brief Build an o2m relation that lets us look up the zones for a node.
 *
 * \note The zone list for each point is not sorted.
 */
template <typename ExecSpace>
class NodeToZoneRelationBuilder
{
public:
  /**
   * \brief Build a node to zone relation and store the resulting O2M relation in the \a relation conduit node.
   *
   * \param topo The topology for which we're building the O2M relation.
   * \param coordset The topology's coordset.
   * \param[out] The node that will contain the O2M relation.
   */
  void execute(const conduit::Node &topo,
               const conduit::Node &coordset,
               conduit::Node &relation)
  {
    const std::string type = topo.fetch_existing("type").as_string();

    // Get the ID of a Conduit allocator that will allocate through Axom with device allocator allocatorID.
    utilities::blueprint::ConduitAllocateThroughAxom<ExecSpace> c2a;
    const int conduitAllocatorID = c2a.getConduitAllocatorID();

    conduit::Node &n_zones = relation["zones"];
    conduit::Node &n_sizes = relation["sizes"];
    conduit::Node &n_offsets = relation["offsets"];
    n_zones.set_allocator(conduitAllocatorID);
    n_sizes.set_allocator(conduitAllocatorID);
    n_offsets.set_allocator(conduitAllocatorID);

    if(type == "unstructured")
    {
      conduit::blueprint::mesh::utils::ShapeType shape(topo);
      const conduit::Node &n_connectivity = topo["elements/connectivity"];
      const std::string shapeType = topo["elements/shape"].as_string();
      const auto intTypeId = n_connectivity.dtype().id();
      const auto connSize = n_connectivity.dtype().number_of_elements();

      // Use the coordset to get the number of nodes. Conduit should be able to do this using only metadata.
      const auto nnodes =
        conduit::blueprint::mesh::utils::coordset::length(coordset);

      if(shape.is_polyhedral())
      {
        using reduce_policy =
          typename axom::execution_space<ExecSpace>::reduce_policy;
        const auto allocatorID = axom::execution_space<ExecSpace>::allocatorID();

        views::dispatch_unstructured_polyhedral_topology(
          topo,
          [&](auto AXOM_UNUSED_PARAM(shape), auto topoView) {
            const auto nzones = topoView.numberOfZones();
            axom::Array<axom::IndexType> sizes(nzones, nzones, allocatorID);
            auto sizes_view = sizes.view();

            // Run through the topology once to do a count of each zone's unique node ids.
            using ZoneType = typename decltype(topoView)::ShapeType;
            RAJA::ReduceSum<reduce_policy, axom::IndexType> count(0);
            topoView.template for_all_zones<ExecSpace>(
              AXOM_LAMBDA(axom::IndexType zoneIndex, const ZoneType &zone) {
                const auto uniqueIds = zone.getUniqueIds();
                sizes_view[zoneIndex] = uniqueIds.size();
                count += uniqueIds.size();
              });
            const auto connSize = count.get();

            // Do a scan on the size array to build an offset array.
            axom::Array<axom::IndexType> offsets(nzones, nzones, allocatorID);
            auto offsets_view = offsets.view();
            axom::exclusive_scan<ExecSpace>(sizes_view, offsets_view);
            sizes.clear();

            // Allocate Conduit arrays on the device in a data type that matches the connectivity.
            conduit::Node n_conn;
            n_conn.set_allocator(conduitAllocatorID);
            n_conn.set(conduit::DataType(intTypeId, connSize));

            n_zones.set(conduit::DataType(intTypeId, connSize));
            n_sizes.set(conduit::DataType(intTypeId, nnodes));
            n_offsets.set(conduit::DataType(intTypeId, nnodes));

            views::IndexNode_to_ArrayView_same(
              n_conn,
              n_zones,
              n_sizes,
              n_offsets,
              [&](auto connectivityView,
                  auto zonesView,
                  auto sizesView,
                  auto offsetsView) {
                // Run through the data one more time to build the nodes and zones arrays.
                topoView.template for_all_zones<ExecSpace>(
                  AXOM_LAMBDA(axom::IndexType zoneIndex, const ZoneType &zone) {
                    const auto uniqueIds = zone.getUniqueIds();
                    auto destIdx = offsets_view[zoneIndex];
                    for(axom::IndexType i = 0; i < uniqueIds.size();
                        i++, destIdx++)
                    {
                      connectivityView[destIdx] = uniqueIds[i];
                      zonesView[destIdx] = zoneIndex;
                    }
                  });

                // Make the relation.
                using ViewType = decltype(connectivityView);
                details::BuildRelation<ExecSpace, ViewType>::execute(
                  connectivityView,
                  zonesView,
                  sizesView,
                  offsetsView);
              });
          });
      }
      else if(shape.is_polygonal() || shapeType == "mixed")
      {
        const conduit::Node &n_topo_sizes = topo["elements/sizes"];
        const conduit::Node &n_topo_offsets = topo["elements/offsets"];

        const auto nzones = n_topo_sizes.dtype().number_of_elements();

        // Allocate Conduit arrays on the device in a data type that matches the connectivity.
        n_zones.set(conduit::DataType(intTypeId, connSize));
        n_sizes.set(conduit::DataType(intTypeId, nnodes));
        n_offsets.set(conduit::DataType(intTypeId, nnodes));

        // Make zones for each node
        views::IndexNode_to_ArrayView_same(
          n_zones,
          n_topo_sizes,
          n_topo_offsets,
          [&](auto zonesView, auto sizesView, auto offsetsView) {
            using DataType = typename decltype(zonesView)::value_type;
            axom::for_all<ExecSpace>(
              0,
              nzones,
              AXOM_LAMBDA(axom::IndexType zoneIndex) {
                for(DataType i = 0; i < sizesView[zoneIndex]; i++)
                  zonesView[offsetsView[zoneIndex] + i] = zoneIndex;
              });
          });

        views::IndexNode_to_ArrayView_same(
          n_connectivity,
          n_zones,
          n_sizes,
          n_offsets,
          [&](auto connectivityView, auto zonesView, auto sizesView, auto offsetsView) {
            // Make the relation.
            using ViewType = decltype(connectivityView);
            details::BuildRelation<ExecSpace, ViewType>::execute(connectivityView,
                                                                 zonesView,
                                                                 sizesView,
                                                                 offsetsView);
          });
      }
      else
      {
        // Shapes are all the same size.
        const auto nodesPerShape = shape.indices;

        // Allocate Conduit arrays on the device in a data type that matches the connectivity.
        n_zones.set(conduit::DataType(intTypeId, connSize));
        n_sizes.set(conduit::DataType(intTypeId, nnodes));
        n_offsets.set(conduit::DataType(intTypeId, nnodes));

        views::IndexNode_to_ArrayView_same(
          n_connectivity,
          n_zones,
          n_sizes,
          n_offsets,
          [&](auto connectivityView, auto zonesView, auto sizesView, auto offsetsView) {
            // Make zones for each node
            axom::for_all<ExecSpace>(
              0,
              connSize,
              AXOM_LAMBDA(axom::IndexType index) {
                zonesView[index] = index / nodesPerShape;
              });
            // Make the relation.
            using ViewType = decltype(connectivityView);
            details::BuildRelation<ExecSpace, ViewType>::execute(connectivityView,
                                                                 zonesView,
                                                                 sizesView,
                                                                 offsetsView);
          });
      }
    }
    else
    {
      // These are all structured topos of some sort. Make an unstructured representation and recurse.

      conduit::Node mesh;
      axom::mir::utilities::blueprint::to_unstructured<ExecSpace>(topo,
                                                                  coordset,
                                                                  "newtopo",
                                                                  mesh);

      // Recurse using the unstructured mesh.
      execute(mesh.fetch_existing("topologies/newtopo"), coordset, relation);
    }
  }
};

}  // end namespace blueprint
}  // end namespace utilities
}  // end namespace mir
}  // end namespace axom

#endif

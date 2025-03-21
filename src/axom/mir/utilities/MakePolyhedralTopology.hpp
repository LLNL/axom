// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_MAKE_POLYHEDRAL_TOPOLOGY_HPP_
#define AXOM_MIR_MAKE_POLYHEDRAL_TOPOLOGY_HPP_

#include "axom/core.hpp"
#include "axom/mir/utilities/utilities.hpp"
#include "axom/mir/utilities/blueprint_utilities.hpp"

#include <conduit/conduit.hpp>

// RAJA
#if defined(AXOM_USE_RAJA)
  #include "RAJA/RAJA.hpp"
#endif

namespace axom
{
namespace mir
{
namespace utilities
{
namespace blueprint
{
/**
 * \brief Make a polyhedral unstructured representation of a topology.
 *
 * \tparam ExecSpace The execution space where the algorithm runs.
 * \tparam TopologyView The topology view type.
 */
template <typename ExecSpace, typename TopologyView>
class MakePolyhedralTopology
{
public:
  using ConnectivityType = typename TopologyView::ConnectivityType;

  static_assert(!TopologyView::ShapeType::is_polyhedral(),
                "MakePolyhedralTopology does not operate on polyhedral topologies");

  /*!
   * \brief Constructor.
   *
   * \param topologyView The topology view that wraps the input topology.
   */
  MakePolyhedralTopology(const TopologyView &topologyView) : m_topologyView(topologyView)
  {
  }

  /**
   * \brief Make a polyhedral unstructured representation of a topology.
   *
   * \tparam ExecSpace The execution space where the work will be done.
   *
   * \param n_topo The input topology to be turned polyhedral.
   * \param[out] n_newTopo The node that will contain the new polyhedral topology.
   *
   */
  void execute(const conduit::Node &n_topo,
               conduit::Node &n_newTopo) const
  {
    AXOM_ANNOTATE_SCOPE("MakePolyhedralTopology");
    namespace bputils = axom::mir::utilities::blueprint;
    using reduce_policy = typename axom::execution_space<ExecSpace>::reduce_policy;
    bputils::ConduitAllocateThroughAxom<ExecSpace> c2a;

    const auto allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    const auto nzones = m_topologyView.numberOfZones();

    //--------------------------------------------------------------------------
    AXOM_ANNOTATE_BEGIN("counting");
    n_newTopo["type"] = "unstructured";
    n_newTopo["coordset"] = n_topo["coordset"].as_string();
    n_newTopo["elements/shape"] = "polyhedral";
    n_newTopo["subelements/shape"] = "polygonal";

    // This node is the number of faces in each zone.
    conduit::Node &n_elem_sizes = n_newTopo["elements/sizes"];
    n_elem_sizes.set_allocator(c2a.getConduitAllocatorID());
    n_elem_sizes.set(
      conduit::DataType(bputils::cpp2conduit<ConnectivityType>::id, nzones));
    auto elem_sizes = bputils::make_array_view<ConnectivityType>(n_elem_sizes);

    conduit::Node &n_elem_offsets = n_newTopo["elements/offsets"];
    n_elem_offsets.set_allocator(c2a.getConduitAllocatorID());
    n_elem_offsets.set(
      conduit::DataType(bputils::cpp2conduit<ConnectivityType>::id, nzones));
    auto elem_offsets = bputils::make_array_view<ConnectivityType>(n_elem_offsets);

    // Compute the total number of faces if they were all unique.
    RAJA::ReduceSum<reduce_policy, axom::IndexType> reduceTotalFaces(0), reduceTotalFaceStorage(0);
    axom::Array<axom::IndexType> zoneFaceSizes(nzones, nzones, allocatorID);
    axom::Array<axom::IndexType> zoneFaceOffsets(nzones, nzones, allocatorID);
    auto zoneFaceSizesView = zoneFaceSizes.view();
    auto zoneFaceOffsetsView = zoneFaceOffsets.view();
    const TopologyView deviceTopologyView(m_topologyView);
    axom::for_all<ExecSpace>(nzones, AXOM_LAMBDA(axom::IndexType zoneIndex)
    {
      const auto zone = deviceTopologyView.zone(zoneIndex);
      elem_sizes[zoneIndex] = zone.numberOfFaces();
      reduceTotalFaces += zone.numberOfFaces();

      axom::IndexType faceStorage = 0;
      for(axom::IndexType fi = 0; fi < zone.numberOfFaces(); fi++)
      {
        faceStorage += zone.numberOfNodesInFace(fi);
      }
      zoneFaceSizesView[zoneIndex] = faceStorage;
      reduceTotalFaceStorage += faceStorage;
    });
    axom::exclusive_scan<ExecSpace>(elem_sizes, elem_offsets);
    axom::exclusive_scan<ExecSpace>(zoneFaceSizesView, zoneFaceOffsetsView);
    zoneFaceSizes.clear();
    const axom::IndexType totalFaces = reduceTotalFaces.get();
    const axom::IndexType totalFaceStorage = reduceTotalFaceStorage.get();
    AXOM_ANNOTATE_END("counting");

    //--------------------------------------------------------------------------
    AXOM_ANNOTATE_BEGIN("elements");

    conduit::Node &n_elem_conn = n_newTopo["elements/connectivity"];
    n_elem_conn.set_allocator(c2a.getConduitAllocatorID());
    n_elem_conn.set(
      conduit::DataType(bputils::cpp2conduit<ConnectivityType>::id, totalFaces));
    auto elem_conn = bputils::make_array_view<ConnectivityType>(n_elem_conn);
    axom::for_all<ExecSpace>(totalFaces, AXOM_LAMBDA(axom::IndexType faceIndex)
    {
      elem_conn[faceIndex] = faceIndex;
    });
    AXOM_ANNOTATE_END("elements");

    //--------------------------------------------------------------------------
    AXOM_ANNOTATE_BEGIN("subelements");

    // Allocate subelement connectivity
    conduit::Node &n_se_conn = n_newTopo["subelements/connectivity"];
    n_se_conn.set_allocator(c2a.getConduitAllocatorID());
    n_se_conn.set(
      conduit::DataType(bputils::cpp2conduit<ConnectivityType>::id, totalFaceStorage));
    auto se_conn = bputils::make_array_view<ConnectivityType>(n_se_conn);

    conduit::Node &n_se_sizes = n_newTopo["subelements/sizes"];
    n_se_sizes.set_allocator(c2a.getConduitAllocatorID());
    n_se_sizes.set(
      conduit::DataType(bputils::cpp2conduit<ConnectivityType>::id, totalFaces));
    auto se_sizes = bputils::make_array_view<ConnectivityType>(n_se_sizes);

    conduit::Node &n_se_offsets = n_newTopo["subelements/offsets"];
    n_se_offsets.set_allocator(c2a.getConduitAllocatorID());
    n_se_offsets.set(
      conduit::DataType(bputils::cpp2conduit<ConnectivityType>::id, totalFaces));
    auto se_offsets = bputils::make_array_view<ConnectivityType>(n_se_offsets);

    // Populate subelement connectivity and make names for the faces.
    RAJA::ReduceSum<reduce_policy, axom::IndexType> reduceTotalStorage(0);
    axom::for_all<ExecSpace>(nzones, AXOM_LAMBDA(axom::IndexType zoneIndex)
    {
      const auto zone = deviceTopologyView.zone(zoneIndex);
      // This where the zone's faces begin in se_conn. We'll update it as we add faces.
      auto offset = zoneFaceOffsetsView[zoneIndex];

      for(axom::IndexType fi = 0; fi < zone.numberOfFaces(); fi++)
      {
        // Where this face begins in se_conn.
        auto faceIds = se_conn.data() + offset;
        // Load the face's ids into faceIds in se_conn.
        axom::IndexType numFaceIds;
        zone.getFace(fi, faceIds, numFaceIds);
        offset += numFaceIds;

        // Store the size of this face.
        const auto thisFaceIndex = elem_offsets[zoneIndex] + fi;
        se_sizes[thisFaceIndex] = numFaceIds;
      }
    });
    axom::exclusive_scan<ExecSpace>(se_sizes, se_offsets);
    AXOM_ANNOTATE_END("subelements");
  }

  TopologyView m_topologyView;
};

}  // end namespace blueprint
}  // end namespace utilities
}  // end namespace mir
}  // end namespace axom

#endif

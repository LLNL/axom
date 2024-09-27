// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_ZONELIST_BUILDER_HPP
#define AXOM_MIR_ZONELIST_BUILDER_HPP

#include <axom/core.hpp>
#include <axom/mir.hpp>

#include <conduit.hpp>

namespace axom
{
namespace mir
{
namespace utilities
{
namespace blueprint
{
/*!
 * \brief This struct builds lists of clean and mixed zones using the input topology and matset views.
 *
 * \tparam ExecSpace The execution space where the algorithm will run.
 * \tparam TopologyView The topology view type on which the algorithm will run.
 * \tparam MatsetView The matset view type on which the algorithm will run.
 */
template <typename ExecSpace, typename TopologyView, typename MatsetView>
class ZoneListBuilder
{
public:
  using SelectedZonesView = axom::ArrayView<axom::IndexType>;
  using ZoneType = typename TopologyView::ShapeType;

  /*!
   * \brief Constructor
   *
   * \param topoView The topology view to use for creating the zone lists.
   * \param matsetView The matset view to use for creating the zone lists.
   */
  ZoneListBuilder(const TopologyView &topoView, const MatsetView &matsetView)
    : m_topologyView(topoView)
    , m_matsetView(matsetView)
  { }

  /*!
   * \brief Build the list of clean and mixed zones using the number of materials
   *        per zone, maxed to the nodes.
   *
   * \param nnodes The number of nodes in the topology's coordset.
   * \param[out] cleanIndices An array that will contain the list of clean material zone ids.
   * \param[out] mixedIndices An array that will contain the list of mixed material zone ids.
   *
   * \note The clean/mixed index arrays are not strictly what could be determined by the matset alone.
   *       We figure out which nodes touch multiple materials. Then we iterate over the zones and
   *       those that touch only nodes that have 1 material are marked clean, otherwise they are
   *       considered mixed as we might have to split those zones.
   */
  void execute(axom::IndexType nnodes,
               axom::Array<axom::IndexType> &cleanIndices,
               axom::Array<axom::IndexType> &mixedIndices) const
  {
    using atomic_policy =
      typename axom::execution_space<ExecSpace>::atomic_policy;
    using reduce_policy =
      typename axom::execution_space<ExecSpace>::reduce_policy;

    AXOM_ANNOTATE_SCOPE("ZoneListBuilder");
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    AXOM_ANNOTATE_BEGIN("nMatsPerNode");
    axom::Array<int> nMatsPerNode(nnodes, nnodes, allocatorID);
    nMatsPerNode.fill(0);
    auto nMatsPerNodeView = nMatsPerNode.view();

    // Determine max number of materials a node might touch.
    MatsetView deviceMatsetView(m_matsetView);
    m_topologyView.template for_all_zones<ExecSpace>(
      AXOM_LAMBDA(axom::IndexType zoneIndex, const ZoneType &zone) {
        const int nmats = deviceMatsetView.numberOfMaterials(zoneIndex);
        const auto nnodesThisZone = zone.numberOfNodes();
        int *nodeData = nMatsPerNodeView.data();
        for(axom::IndexType i = 0; i < nnodesThisZone; i++)
        {
          const auto nodeId = zone.getId(i);
          int *nodePtr = nodeData + nodeId;
          RAJA::atomicMax<atomic_policy>(nodePtr, nmats);
        }
      });
    AXOM_ANNOTATE_END("nMatsPerNode");

    // Now, mark all zones that have 1 mat per node as clean.
    AXOM_ANNOTATE_BEGIN("mask");
    const auto nzones = m_topologyView.numberOfZones();
    axom::Array<int> mask(nzones, nzones, allocatorID);
    auto maskView = mask.view();
    RAJA::ReduceSum<reduce_policy, int> mask_reduce(0);
    m_topologyView.template for_all_zones<ExecSpace>(
      AXOM_LAMBDA(axom::IndexType zoneIndex, const ZoneType &zone) {
        bool clean = true;
        const axom::IndexType nnodesThisZone = zone.numberOfNodes();
        for(axom::IndexType i = 0; i < nnodesThisZone && clean; i++)
        {
          const auto nodeId = zone.getId(i);
          clean &= (nMatsPerNodeView[nodeId] == 1);
        }

        const int ival = clean ? 1 : 0;
        maskView[zoneIndex] = ival;
        mask_reduce += ival;
      });
    AXOM_ANNOTATE_END("mask");

    const int nClean = mask_reduce.get();
    if(nClean > 0)
    {
      AXOM_ANNOTATE_BEGIN("offsets");
      axom::Array<int> maskOffsets(nzones, nzones, allocatorID);
      auto maskOffsetsView = maskOffsets.view();
      axom::exclusive_scan<ExecSpace>(maskView, maskOffsetsView);
      AXOM_ANNOTATE_END("offsets");

      // Make the output cleanIndices array.
      AXOM_ANNOTATE_BEGIN("cleanIndices");
      cleanIndices = axom::Array<axom::IndexType>(nClean, nClean, allocatorID);
      auto cleanIndicesView = cleanIndices.view();
      axom::for_all<ExecSpace>(
        nzones,
        AXOM_LAMBDA(axom::IndexType index) {
          if(maskView[index] > 0)
          {
            cleanIndicesView[maskOffsetsView[index]] = index;
          }
        });
      AXOM_ANNOTATE_END("cleanIndices");

      // Make the mixedIndices array.
      AXOM_ANNOTATE_BEGIN("mixedIndices");
      axom::for_all<ExecSpace>(
        nzones,
        AXOM_LAMBDA(axom::IndexType index) {
          maskView[index] = (maskView[index] == 1) ? 0 : 1;
        });
      axom::exclusive_scan<ExecSpace>(maskView, maskOffsetsView);
      const int nMixed = nzones - nClean;
      mixedIndices = axom::Array<axom::IndexType>(nMixed, nMixed, allocatorID);
      auto mixedIndicesView = mixedIndices.view();
      axom::for_all<ExecSpace>(
        nzones,
        AXOM_LAMBDA(axom::IndexType index) {
          if(maskView[index] > 0)
          {
            mixedIndicesView[maskOffsetsView[index]] = index;
          }
        });
      AXOM_ANNOTATE_END("mixedIndices");
    }
    else
    {
      AXOM_ANNOTATE_SCOPE("mixedIndices");
      cleanIndices = axom::Array<axom::IndexType>();

      // There were no clean, so it must all be mixed.
      mixedIndices = axom::Array<axom::IndexType>(nzones, nzones, allocatorID);
      auto mixedIndicesView = mixedIndices.view();
      axom::for_all<ExecSpace>(
        nzones,
        AXOM_LAMBDA(axom::IndexType index) { mixedIndicesView[index] = index; });
    }
  }

  /*!
   * \brief Build the list of clean and mixed zones using the number of materials
   *        per zone, maxed to the nodes. Limit the number of zones.
   *
   * \param nnodes The number of nodes in the topology's coordset.
   * \param selectedZonesView A view containing the zone indices we're considering.
   * \param[out] cleanIndices An array that will contain the list of clean material zone ids.
   * \param[out] mixedIndices An array that will contain the list of mixed material zone ids.
   *
   * \note The clean/mixed index arrays are not strictly what could be determined by the matset alone.
   *       We figure out which nodes touch multiple materials. Then we iterate over the zones and
   *       those that touch only nodes that have 1 material are marked clean, otherwise they are
   *       considered mixed as we might have to split those zones.
   */
  void execute(axom::IndexType nnodes,
               const SelectedZonesView &selectedZonesView,
               axom::Array<axom::IndexType> &cleanIndices,
               axom::Array<axom::IndexType> &mixedIndices) const
  {
    using atomic_policy =
      typename axom::execution_space<ExecSpace>::atomic_policy;
    using reduce_policy =
      typename axom::execution_space<ExecSpace>::reduce_policy;

    AXOM_ANNOTATE_SCOPE("ZoneListBuilder");
    SLIC_ASSERT(selectedZonesView.size() > 0);

    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    AXOM_ANNOTATE_BEGIN("nMatsPerNode");
    axom::Array<int> nMatsPerNode(nnodes, nnodes, allocatorID);
    nMatsPerNode.fill(0);
    auto nMatsPerNodeView = nMatsPerNode.view();

    // Determine max number of materials a node might touch.
    MatsetView deviceMatsetView(m_matsetView);
    m_topologyView.template for_selected_zones<ExecSpace>(
      selectedZonesView,
      AXOM_LAMBDA(axom::IndexType AXOM_UNUSED_PARAM(szIndex), axom::IndexType zoneIndex, const ZoneType &zone) {
        const int nmats = deviceMatsetView.numberOfMaterials(zoneIndex);
        const auto nnodesThisZone = zone.numberOfNodes();
        int *nodeData = nMatsPerNodeView.data();
        for(axom::IndexType i = 0; i < nnodesThisZone; i++)
        {
          const auto nodeId = zone.getId(i);
          int *nodePtr = nodeData + nodeId;
          RAJA::atomicMax<atomic_policy>(nodePtr, nmats);
        }
      });
    AXOM_ANNOTATE_END("nMatsPerNode");

    // Now, mark all selected zones that have 1 mat per node as clean.
    AXOM_ANNOTATE_BEGIN("mask");
    const auto nzones = selectedZonesView.size();
    axom::Array<int> mask(nzones, nzones, allocatorID);
    auto maskView = mask.view();
    RAJA::ReduceSum<reduce_policy, int> mask_reduce(0);
    m_topologyView.template for_selected_zones<ExecSpace>(
      selectedZonesView,
      AXOM_LAMBDA(axom::IndexType szIndex, axom::IndexType AXOM_UNUSED_PARAM(zoneIndex), const ZoneType &zone) {
        bool clean = true;
        const axom::IndexType nnodesThisZone = zone.numberOfNodes();
        for(axom::IndexType i = 0; i < nnodesThisZone && clean; i++)
        {
          const auto nodeId = zone.getId(i);
          clean &= (nMatsPerNodeView[nodeId] == 1);
        }

        const int ival = clean ? 1 : 0;
        maskView[szIndex] = ival;
        mask_reduce += ival;
      });
    AXOM_ANNOTATE_END("mask");

    const int nClean = mask_reduce.get();
    if(nClean > 0)
    {
      AXOM_ANNOTATE_BEGIN("offsets");
      axom::Array<int> maskOffsets(nzones, nzones, allocatorID);
      auto maskOffsetsView = maskOffsets.view();
      axom::exclusive_scan<ExecSpace>(maskView, maskOffsetsView);
      AXOM_ANNOTATE_END("offsets");

      // Make the output cleanIndices array.
      AXOM_ANNOTATE_BEGIN("cleanIndices");
      cleanIndices = axom::Array<axom::IndexType>(nClean, nClean, allocatorID);
      auto cleanIndicesView = cleanIndices.view();
      axom::for_all<ExecSpace>(
        nzones,
        AXOM_LAMBDA(axom::IndexType index) {
          if(maskView[index] > 0)
          {
            cleanIndicesView[maskOffsetsView[index]] = selectedZonesView[index];
          }
        });
      AXOM_ANNOTATE_END("cleanIndices");

      // Make the mixedIndices array.
      AXOM_ANNOTATE_BEGIN("mixedIndices");
      axom::for_all<ExecSpace>(
        nzones,
        AXOM_LAMBDA(axom::IndexType index) {
          maskView[index] = (maskView[index] == 1) ? 0 : 1;
        });
      axom::exclusive_scan<ExecSpace>(maskView, maskOffsetsView);
      const int nMixed = nzones - nClean;
      mixedIndices = axom::Array<axom::IndexType>(nMixed, nMixed, allocatorID);
      auto mixedIndicesView = mixedIndices.view();
      axom::for_all<ExecSpace>(
        nzones,
        AXOM_LAMBDA(axom::IndexType index) {
          if(maskView[index] > 0)
          {
            mixedIndicesView[maskOffsetsView[index]] = selectedZonesView[index];
          }
        });
      AXOM_ANNOTATE_END("mixedIndices");
    }
    else
    {
      AXOM_ANNOTATE_SCOPE("mixedIndices");
      cleanIndices = axom::Array<axom::IndexType>();

      // There were no clean, so it must all be mixed.
      mixedIndices = axom::Array<axom::IndexType>(nzones, nzones, allocatorID);
      auto mixedIndicesView = mixedIndices.view();
      axom::for_all<ExecSpace>(
        nzones,
        AXOM_LAMBDA(axom::IndexType index) {
          mixedIndicesView[index] = selectedZonesView[index];
        });
    }
  }

private:
  TopologyView m_topologyView;
  MatsetView m_matsetView;
};

}  // end namespace blueprint
}  // end namespace utilities
}  // end namespace mir
}  // end namespace axom

#endif
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

/**
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
  /**
   * \brief Constructor
   *
   * \param topoView The topology view to use for creating the zone lists.
   * \param matsetView The matset view to use for creating the zone lists.
   */
  ZoneListBuilder(const TopologyView &topoView, const MatsetView &matsetView) : m_topologyView(topoView), m_matsetView(matsetView)
  {
  }

  /**
   * \brief Build the list of clean and mixed zones.
   *
   * \param n2z A Conduit node containing a node to zone relation.
   * \param[out] cleanIndices An array that will contain the list of clean material zone ids.
   * \param[out] mixedIndices An array that will contain the list of mixed material zone ids.
   *
   * \note The clean/mixed index arrays are not strictly what could be determined by the matset alone.
   *       We figure out which nodes touch multiple materials. Then we iterate over the zones and
   *       those that touch only nodes that have 1 material are marked clean, otherwise they are
   *       considered mixed as we might have to split those zones.
   */
  void execute(const conduit::Node &n2z, axom::Array<axom::IndexType> &cleanIndices, axom::Array<axom::IndexType> &mixedIndices) const
  {
    AXOM_ANNOTATE_SCOPE("ZoneListBuilder");
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    const conduit::Node &n_zones = n2z.fetch_existing("zones");
    const conduit::Node &n_indices = n2z.fetch_existing("indices");
    const conduit::Node &n_sizes = n2z.fetch_existing("sizes");
    const conduit::Node &n_offsets = n2z.fetch_existing("offsets");
    axom::mir::views::IndexNode_to_ArrayView_same(n_zones, n_indices, n_sizes, n_offsets, [&](auto zones, auto indices, auto sizes, auto offsets)
    {
      using reduce_policy = typename axom::execution_space<ExecSpace>::reduce_policy;
      const auto nnodes = sizes.size();
      axom::Array<int> nMatsPerNode(nnodes, nnodes, allocatorID);
      auto nMatsPerNodeView = nMatsPerNode.view();

      // OR to the nodes. We just take the max number of materials of adjacent
      // zones and put those on the nodes lest we need an array that stores
      // general N-bit values per element.
      MatsetView deviceMatsetView(m_matsetView);
      axom::for_all<ExecSpace>(nnodes, AXOM_LAMBDA(auto index)
      {
        int nmats = 0;
        const int size = static_cast<int>(sizes[index]);
        const auto offset = offsets[index];
        for(int i = 0; i < size; i++)
        {
          const auto zi = zones[indices[offset + i]];
          nmats = axom::utilities::max(nmats, static_cast<int>(deviceMatsetView.numberOfMaterials(zi)));
        }
        nMatsPerNodeView[index] = nmats;
      });

      // Now, mark all zones that have 1 mat per node as clean.
      const auto nzones = m_topologyView.numberOfZones();
      axom::Array<int> mask(nzones, nzones, allocatorID);
      auto maskView = mask.view();
      RAJA::ReduceSum<reduce_policy, int> mask_reduce(0);
      m_topologyView.template for_all_zones<ExecSpace>(AXOM_LAMBDA(auto zoneIndex, const auto &zone)
      {
        bool clean = true;
        const axom::IndexType nnodesThisZone = zone.numberOfNodes();
        for(axom::IndexType i = 0; i < nnodesThisZone && clean; i++)
        {
          const auto nodeId = zone.getId(i);
          clean |= nMatsPerNodeView[nodeId] == 1;
        }

        const int ival = clean ? 1 : 0;
        maskView[zoneIndex] = ival;
        mask_reduce += ival;
      });

      const int nClean = mask_reduce.get();
      if(nClean > 0)
      {
        axom::Array<int> maskOffsets(nzones, nzones, allocatorID);
        auto maskOffsetsView = maskOffsets.view();
        axom::exclusive_scan<ExecSpace>(maskView, maskOffsetsView);

        // Make the output cleanIndices array.
        cleanIndices = axom::Array<axom::IndexType>(nClean, nClean, allocatorID);
        auto cleanIndicesView = cleanIndices.view();
        axom::for_all<ExecSpace>(nzones, AXOM_LAMBDA(auto index)
        {
          if(maskView[index] > 0)
          {
            cleanIndicesView[maskOffsetsView[index]] = index;
          }
        });

        // Make the mixedIndices array.
        axom::for_all<ExecSpace>(nzones, AXOM_LAMBDA(auto index)
        {
          maskView[index] = ~maskView[index];
        });
        axom::exclusive_scan<ExecSpace>(maskView, maskOffsetsView);
        const int nMixed = nzones - nClean;
        mixedIndices = axom::Array<axom::IndexType>(nMixed, nMixed, allocatorID);
        auto mixedIndicesView = mixedIndices.view();
        axom::for_all<ExecSpace>(nzones, AXOM_LAMBDA(auto index)
        {
          if(maskView[index] > 0)
          {
            mixedIndicesView[maskOffsetsView[index]] = index;
          }
        });
      }
      else
      {
        cleanIndices = axom::Array<axom::IndexType>();

        // There were no clean, so it must all be mixed.
        mixedIndices = axom::Array<axom::IndexType>(nzones, nzones, allocatorID);
        auto mixedIndicesView = mixedIndices.view();
        axom::for_all<ExecSpace>(nzones, AXOM_LAMBDA(auto index)
        {
          mixedIndicesView[index] = index;
        });
      }
    });
  }

private:
  TopologyView m_topologyView;
  MatsetView   m_matsetView;
};

} // end namespace blueprint
} // end namespace utilities
} // end namespace mir
} // end namespace axom

#endif

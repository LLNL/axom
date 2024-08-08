// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_BLEND_GROUP_BUILDER_HPP_
#define AXOM_MIR_BLEND_GROUP_BUILDER_HPP_

#include "axom/core.hpp"
#include "axom/mir/utilities.hpp"

namespace axom
{
namespace mir
{
namespace clipping
{

/**
 * \brief This class encapsulates the logic for building blend groups.
 *
 * \tparam ExecSpace The execution space where the algorithm will run.
 * \tparam NamingPolicy The naming policy used to create names from node ids.
 *
 * \note The class only contains views to data so it can be copied into lambdas.
 */
template <typename ExecSpace, typename NamingPolicy = axom::mir::utilities::HashNaming>
class BlendGroupBuilder
{
public:
  using KeyType = typename NamingPolicy::KeyType;

  /**
   * \brief This struct holds the views that represent data for blend groups.
   */
  struct State
  {
    // clang-format off
    IndexType m_nzones;

    axom::ArrayView<IndexType> m_blendGroupsView;        // Number of blend groups in each zone.
    axom::ArrayView<IndexType> m_blendGroupsLenView;     // total size of blend group data for each zone.
    axom::ArrayView<IndexType> m_blendOffsetView;        // The offset of each zone's blend groups.
    axom::ArrayView<IndexType> m_blendGroupOffsetsView;  // Start of each zone's blend group data.

    axom::ArrayView<KeyType>   m_blendNamesView;         // Blend group names
    axom::ArrayView<IndexType> m_blendGroupSizesView;    // Size of individual blend group.
    axom::ArrayView<IndexType> m_blendGroupStartView;    // Start of individual blend group's data in m_blendIdsView/m_blendCoeffView.
    axom::ArrayView<IndexType> m_blendIdsView;           // blend group ids.
    axom::ArrayView<float>     m_blendCoeffView;         // blend group weights.

    axom::ArrayView<KeyType>   m_blendUniqueNamesView;   // Unique names of blend groups.
    axom::ArrayView<IndexType> m_blendUniqueIndicesView; // Indices of the unique names in blend group definitions.
    // clang-format on
  };

  /**
   * \brief Access the state views.
   * \return A reference to the state.
   */
  State &state() { return m_state; }
  const State &state() const { return m_state; }

  /**
   * \brief Set the number of zones.
   *
   * \param blendGroupsView The view that holds the number of blend groups for each zone.
   * \param blendGroupsLenView The view that holds the size of the blend group data for each zone.
   */
  void setBlendGroupSizes(const axom::ArrayView<IndexType> &blendGroupsView,
                          const axom::ArrayView<IndexType> &blendGroupsLenView)
  {
    m_state.m_nzones = blendGroupsView.size();
    m_state.m_blendGroupsView = blendGroupsView;
    m_state.m_blendGroupsLenView = blendGroupsLenView;
  }

  /**
   * \brief Compute the total sizes of blend group storage.
   *
   * \param[out] bgSum The total number of blend groups for all zones.
   * \param[out] bgLenSum The total size of blend group data for all zones.
   */
  void computeBlendGroupSizes(IndexType &bgSum, IndexType &bgLenSum)
  {
    using reduce_policy =
      typename axom::execution_space<ExecSpace>::reduce_policy;
    RAJA::ReduceSum<reduce_policy, IndexType> blendGroups_sum(0);
    RAJA::ReduceSum<reduce_policy, IndexType> blendGroupLen_sum(0);
    const auto localBlendGroupsView = m_state.m_blendGroupsView;
    const auto localBlendGroupsLenView = m_state.m_blendGroupsLenView;
    axom::for_all<ExecSpace>(
      m_state.m_nzones,
      AXOM_LAMBDA(auto zoneIndex) {
        blendGroups_sum += localBlendGroupsView[zoneIndex];
        blendGroupLen_sum += localBlendGroupsLenView[zoneIndex];
      });
    bgSum = blendGroups_sum.get();
    bgLenSum = blendGroupLen_sum.get();
  }

  /**
   * \brief Set the views for the blend group offsets and then fill them using a scan.
   *
   * \param blendOffsetView The offsets to each blend group for views sized: view[blendGroupSum].
   * \param blendGroupOffsetsView The offsets to each zone's blend groups data.
   */
  void setBlendGroupOffsets(const axom::ArrayView<IndexType> &blendOffsetView,
                            const axom::ArrayView<IndexType> &blendGroupOffsetsView)
  {
    m_state.m_blendOffsetView = blendOffsetView;
    m_state.m_blendGroupOffsetsView = blendGroupOffsetsView;
  }

  /**
   * brief Compute the blend group offsets that make it easier to store data.
   */
  void computeBlendGroupOffsets()
  {
    axom::exclusive_scan<ExecSpace>(m_state.m_blendGroupsLenView,
                                    m_state.m_blendOffsetView);
    axom::exclusive_scan<ExecSpace>(m_state.m_blendGroupsView,
                                    m_state.m_blendGroupOffsetsView);
  }

  /**
   * \brief Set the views that we'll use for blend groups.
   */
  void setBlendViews(const axom::ArrayView<KeyType> &blendNames,
                     const axom::ArrayView<IndexType> &blendGroupSizes,
                     const axom::ArrayView<IndexType> &blendGroupStart,
                     const axom::ArrayView<IndexType> &blendIds,
                     const axom::ArrayView<float> &blendCoeff)
  {
    m_state.m_blendNamesView = blendNames;
    m_state.m_blendGroupSizesView = blendGroupSizes;
    m_state.m_blendGroupStartView = blendGroupStart;
    m_state.m_blendIdsView = blendIds;
    m_state.m_blendCoeffView = blendCoeff;
  }

  /**
   * \brief Set the unique names and ids views. These are used in mapping blend groups to unique blend groups.
   *
   * \param uniqueNames A view containing unique, sorted blend group names.
   * \param uniqueIndices A view containing the original blend group index for each unique name.
   */
  void setUniqueNames(const axom::ArrayView<KeyType> &uniqueNames,
                      const axom::ArrayView<IndexType> &uniqueIndices)
  {
    m_state.m_blendUniqueNamesView = uniqueNames;
    m_state.m_blendUniqueIndicesView = uniqueIndices;
  }

  /**
   * \brief Get the blend names view.
   * \return The blend names view.
   */
  const axom::ArrayView<KeyType> &blendNames() const
  {
    return m_state.m_blendNamesView;
  }

  /**
   * \brief This class helps us manage blend group creation and usage for blend groups within a single zone.
   */
  class zone_blend_groups
  {
  public:
    /**
     * \brief Return the number of blend groups for this zone.
     * \return The number of blend groups for this zone.
     */
    AXOM_HOST_DEVICE
    inline IndexType size() const
    {
      return m_state->m_blendGroupsView[m_zoneIndex];
    }

    /**
     * \brief Set the number of blend groups and total size of the blend groups for a zone.
     *
     * \param zoneIndex The index of the zone we're initializing.
     * \param nBlendGroups The number of blend groups in this zone.
     * \param blendGroupSize The size of all of the blend groups in this zone.
     */
    AXOM_HOST_DEVICE
    inline void setSize(IndexType nBlendGroups, IndexType blendGroupsSize)
    {
      m_state->m_blendGroupsView[m_zoneIndex] = nBlendGroups;
      m_state->m_blendGroupsLenView[m_zoneIndex] = blendGroupsSize;
    }

    /**
     * \brief Start creating a new blend group within the allocated space for this zone.
     */
    AXOM_HOST_DEVICE
    inline void beginGroup() { m_currentDataOffset = m_startOffset; }

    /**
     * \brief Add a new blend point in the current blend group.
     *
     * \param id The node id that will be used for blending.
     * \param weight The weight that will be used for blending.
     */
    AXOM_HOST_DEVICE
    inline void add(IndexType id, float weight)
    {
      m_state->m_blendIdsView[m_currentDataOffset] = id;
      m_state->m_blendCoeffView[m_currentDataOffset] = weight;
      m_currentDataOffset++;
    }

    /**
     * \brief End the current blend group, storing its name, size, etc.
     */
    AXOM_HOST_DEVICE
    inline void endGroup()
    {
      IndexType numIds = m_currentDataOffset - m_startOffset;

      // Store where this blendGroup starts in the blendIds,blendCoeff.
      m_state->m_blendGroupStartView[m_blendGroupId] = m_startOffset;

      // Save the size for this blend group.
      m_state->m_blendGroupSizesView[m_blendGroupId] = numIds;

      // Store "name" of blend group.
      KeyType blendName =
        NamingPolicy::makeName(m_state->m_blendIdsView.data() + m_startOffset,
                               numIds);

      m_state->m_blendNamesView[m_blendGroupId] = blendName;
#if defined(AXOM_DEBUG) && !defined(AXOM_DEVICE_CODE)
      const float w = weightSum();
      constexpr float EPS = 1.e-5;
      if(w < (1. - EPS) || w > (1. + EPS))
      {
        SLIC_ERROR(fmt::format("Invalid blend group weights w={}.", w));
        print(std::cout);
      }
#endif
      m_blendGroupId++;
      m_startOffset = m_currentDataOffset;
    }

    /**
     * \brief Return the sum of the weights for the current blend group.
     * \note The blendGroup must have finished construction.
     * \return The sum of weights, which should be about 1.
     */
    AXOM_HOST_DEVICE
    inline float weightSum() const
    {
      const auto numIds = m_state->m_blendGroupSizesView[m_blendGroupId];
      const auto start = m_state->m_blendGroupStartView[m_blendGroupId];
      float w = 0.f;
      for(IndexType i = 0; i < numIds; i++)
        w += m_state->m_blendCoeffView[start + i];
      return w;
    }

    /**
     * \brief Return the name of the current blend group.
     * \return The name of the current blend group.
     */
    AXOM_HOST_DEVICE
    inline KeyType name() const
    {
      return m_state->m_blendNamesView[m_blendGroupId];
    }

    /**
     * \brief Return index of the current blend group in the unique blend groups.
     * \return The unique index of the current blend group.
     */
    AXOM_HOST_DEVICE
    inline IndexType uniqueBlendGroupIndex() const
    {
      return axom::mir::utilities::bsearch(name(),
                                           m_state->m_blendUniqueNamesView);
    }

    /**
     * \brief Advance to the next blend group.
     */
    AXOM_HOST_DEVICE
    inline void operator++() { m_blendGroupId++; }

    /**
     * \brief Advance to the next blend group.
     */
    AXOM_HOST_DEVICE
    inline void operator++(int) { m_blendGroupId++; }

#if !defined(AXOM_DEVICE_CODE)
    /**
     * \brief Print the current blend group to a stream.
     * \param os The stream to which the blend group will print.
     */
    void print(std::ostream &os) const
    {
      const auto n = m_state->m_blendGroupSizesView[m_blendGroupId];
      const auto offset = m_state->m_blendGroupStartView[m_blendGroupId];
      os << "-\n";
      os << " zoneIndex: " << m_zoneIndex << std::endl;
      os << " blendGroupId: " << m_blendGroupId << std::endl;
      os << " size: " << n << std::endl;
      os << " offset: " << offset << std::endl;

      const IndexType *ids = m_state->m_blendIdsView.data() + offset;
      os << " ids: [";
      for(int bi = 0; bi < n; bi++)
      {
        if(bi > 0) os << ", ";
        os << ids[bi];
      }
      os << "]";
      os << "\n";
      const float *weights = m_state->m_blendCoeffView.data() + offset;
      os << " weights: [";
      for(int bi = 0; bi < n; bi++)
      {
        if(bi > 0) os << ", ";
        os << weights[bi];
      }
      os << "]\n";
      os << " name: " << m_state->m_blendNamesView[m_blendGroupId];
      os << "\n";
    }
#endif

    /**
     * \brief Return the current blend group's ids.
     * \return The current blend group's ids.
     * \note This method should not be used if blend groups are still being constructed.
     */
    axom::ArrayView<IndexType> ids() const
    {
      const auto n = m_state->m_blendGroupSizesView[m_blendGroupId];
      const auto offset = m_state->m_blendGroupStartView[m_blendGroupId];
      return axom::ArrayView<IndexType>(m_state->m_blendIdsView.data() + offset,
                                        n);
    }

    /**
     * \brief Return the current blend group's ids.
     * \return The current blend group's ids.
     * \note This method should not be used if blend groups are still being constructed.
     */
    axom::ArrayView<float> weights() const
    {
      const auto n = m_state->m_blendGroupSizesView[m_blendGroupId];
      const auto offset = m_state->m_blendGroupStartView[m_blendGroupId];
      return axom::ArrayView<float>(m_state->m_blendCoeffsView.data() + offset,
                                    n);
    }

  private:
    friend class BlendGroupBuilder;

    IndexType m_zoneIndex;  // The zone that owns this set of blend groups.
    IndexType m_blendGroupId;  // The global blend group index within this current zone.
    IndexType m_startOffset;  // The data offset for the first ids/weights in this blend group.
    IndexType m_currentDataOffset;  // The current data offset.
    State *m_state;                 // Pointer to the main state.
  };

  /**
   * \brief Return a zone_blend_groups object for the current zone so we can add blend groups.
   *
   * \param zoneIndex The zone whose blend groups we want to edit.
   *
   * \note This method must be marked const because we can call it from an AXOM_LAMBDA.
   *       we pass a non-const State reference to the zone_blend_groups that we construct
   *       so we can write into the blend group data.
   */
  AXOM_HOST_DEVICE
  zone_blend_groups blendGroupsForZone(IndexType zoneIndex) const
  {
    zone_blend_groups groups;
    // The zone that owns this set of blend groups.
    groups.m_zoneIndex = zoneIndex;

    // Global blend group id for the first blend group in this zone.
    groups.m_blendGroupId = m_state.m_blendGroupOffsetsView[zoneIndex];
    // Global start
    groups.m_startOffset = groups.m_currentDataOffset =
      m_state.m_blendOffsetView[zoneIndex];

    groups.m_state = const_cast<State *>(&m_state);
    return groups;
  }

  /**
   * \brief Make a BlendData object from the views in this object.
   *
   * \return A BlendData object suitable for making new fields and coordsets.
   */
  axom::mir::utilities::blueprint::BlendData makeBlendData() const
  {
    axom::mir::utilities::blueprint::BlendData blend;

    blend.m_selectedIndicesView =
      m_state.m_blendUniqueIndicesView;  // We'll use these to select just the unique indices
    blend.m_blendGroupSizesView = m_state.m_blendGroupSizesView;
    blend.m_blendGroupStartView = m_state.m_blendGroupStartView;
    blend.m_blendIdsView = m_state.m_blendIdsView;
    blend.m_blendCoeffView = m_state.m_blendCoeffView;

    return blend;
  }

private:
  State m_state;
};

}  // end namespace clipping
}  // end namespace mir
}  // end namespace axom

#endif

// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_FIELD_BLENDER_HPP_
#define AXOM_MIR_FIELD_BLENDER_HPP_

#include "axom/core.hpp"
#include "axom/mir/views/NodeArrayView.hpp"
#include "axom/mir/utilities.hpp"
#include "axom/mir/blueprint_utilities.hpp"

#include <conduit/conduit.hpp>

namespace axom
{
namespace mir
{
namespace utilities
{
namespace blueprint
{
/**
 * \brief This class contains views of blend data.
 */
struct BlendData
{
  axom::ArrayView<IndexType> m_selectedIndicesView;  // Contains indices of the selected blend groups.

  axom::ArrayView<IndexType> m_blendGroupSizesView;  // The number of ids/weights in each blend group.
  axom::ArrayView<IndexType>
    m_blendGroupStartView;  // The starting offset for a blend group in the ids/weights.
  axom::ArrayView<IndexType> m_blendIdsView;  // Contains ids that make up the blend groups
  axom::ArrayView<float> m_blendCoeffView;  // Contains the weights that make up the blend groups.
};

/**
 * \brief This policy can be used with FieldBlender to select all blend groups.
 */
struct SelectAllPolicy
{
  AXOM_HOST_DEVICE
  static IndexType size(const BlendData &blend)
  {
    return blend.m_blendGroupSizesView.size();
  }

  AXOM_HOST_DEVICE
  static IndexType selectedIndex(const BlendData & /*blend*/, IndexType index)
  {
    return index;
  }
};

/**
 * \brief This policy can be used with FieldBlender to select a subset of blend groups, according to m_selectedIndicesView.
 */
struct SelectThroughArrayView
{
  AXOM_HOST_DEVICE
  static IndexType size(const BlendData &blend)
  {
    return blend.m_selectedIndicesView.size();
  }

  AXOM_HOST_DEVICE
  static IndexType selectedIndex(const BlendData &blend, IndexType index)
  {
    return blend.m_selectedIndicesView[index];
  }
};

/**
 * \accelerated
 * \class FieldBlender
 *
 * \brief This class uses BlendData to generate a new blended field from an input field.
 *
 * \tparam ExecSpace The execution space where the work will occur.
 * \tparam SelectionPolicy The selection policy to use.
 * \tparam IndexingPolicy A class that provides operator[] that can transform node indices.
 */
template <typename ExecSpace, typename SelectionPolicy, typename IndexingPolicy = DirectIndexing>
class FieldBlender
{
public:
  /// Constructor
  FieldBlender() : m_indexing() { }

  /**
   * \brief Constructor
   * \param indexing An object used to transform node indices.
   */
  FieldBlender(const IndexingPolicy &indexing) : m_indexing(indexing) { }

  /**
   * \brief Create a new blended field from the \a n_input field and place it in \a n_output.
   *
   * \param blend The BlendData that will be used to make the new field.
   * \param n_input The input field that we're blending.
   * \param n_output The output node that will contain the new field.
   */
  void execute(const BlendData &blend,
               const conduit::Node &n_input,
               conduit::Node &n_output) const
  {
    n_output.reset();
    n_output["association"] = n_input["association"];
    n_output["topology"] = n_input["topology"];

    const conduit::Node &n_input_values = n_input["values"];
    conduit::Node &n_output_values = n_output["values"];
    if(n_input_values.number_of_children() > 0)
    {
      for(conduit::index_t i = 0; i < n_input_values.number_of_children(); i++)
      {
        const conduit::Node &n_comp = n_input_values[i];
        conduit::Node &n_out_comp = n_output_values[n_comp.name()];
        blendSingleComponent(blend, n_comp, n_out_comp);
      }
    }
    else
    {
      blendSingleComponent(blend, n_input_values, n_output_values);
    }
  }

private:
  /**
   * \brief Blend data for a single field component.
   *
   * \param blend The BlendData that will be used to make the new field.
   * \param n_values The input values that we're blending.
   * \param n_output_values The output node that will contain the new field.
   */
  void blendSingleComponent(const BlendData &blend,
                            const conduit::Node &n_values,
                            conduit::Node &n_output_values) const
  {
    // We're allowing selectedIndicesView to be used to select specific blend
    // groups. If the user did not provide that, use all blend groups.
    const auto outputSize = SelectionPolicy::size(blend);

    // Allocate Conduit data through Axom.
    utilities::blueprint::ConduitAllocateThroughAxom<ExecSpace> c2a;
    n_output_values.set_allocator(c2a.getConduitAllocatorID());
    n_output_values.set(conduit::DataType(n_values.dtype().id(), outputSize));

    views::Node_to_ArrayView_same(
      n_values,
      n_output_values,
      [&](auto compView, auto outView) {
        using value_type = typename decltype(compView)::value_type;
        using accum_type =
          typename axom::mir::utilities::accumulation_traits<value_type>::type;

        IndexingPolicy deviceIndexing(m_indexing);
        const BlendData deviceBlend(blend);
        axom::for_all<ExecSpace>(
          outputSize,
          AXOM_LAMBDA(auto bgid) {
            // Get the index we want.
            const auto selectedIndex =
              SelectionPolicy::selectedIndex(deviceBlend, bgid);
            const auto start = deviceBlend.m_blendGroupStartView[selectedIndex];
            const auto end =
              start + deviceBlend.m_blendGroupSizesView[selectedIndex];

            accum_type blended = 0;
            for(IndexType i = start; i < end; i++)
            {
              const auto index = deviceBlend.m_blendIdsView[i];
              const auto weight = deviceBlend.m_blendCoeffView[i];
              const auto transformedIndex = deviceIndexing[index];
              blended +=
                static_cast<accum_type>(compView[transformedIndex]) * weight;
            }
            outView[bgid] = static_cast<value_type>(blended);
          });
      });
  }

private:
  IndexingPolicy m_indexing {};
};

}  // end namespace blueprint
}  // end namespace utilities
}  // end namespace mir
}  // end namespace axom

#endif

// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_FIELD_BLENDER_HPP_
#define AXOM_MIR_FIELD_BLENDER_HPP_

#include "axom/core.hpp"
#include "axom/mir/views/NodeArrayView.hpp"

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
  axom::ArrayView<uint64> m_blendNamesView;      // Contains unique hash ids for the sorted ids in each blend group.
  axom::ArrayView<int32>  m_uniqueIndicesView;   // Contains the index into the original blend group data for a unique blend group.

  axom::ArrayView<int32>  m_blendGroupSizesView; // The number of ids/weights in each blend group.
  axom::ArrayView<int32>  m_blendGroupStartView; // The starting offset for a blend group in the ids/weights.
  axom::ArrayView<int32>  m_blendIdsView;        // Contains ids that make up the blend groups
  axom::ArrayView<float>  m_blendCoeffView;      // Contains the weights that make up the blend groups.
};

/**
 * \accelerated
 * \class FieldBlender
 *
 * \brief This class uses BlendData to generate a new blended field from an input field.
 */
template <typename ExecSpace>
class FieldBlender
{ 
public:
  /**
   * \brief Create a new blended field from the \a n_input field and place it in \a n_output.
   *
   * \param blend The BlendData that will be used to make the new field.
   * \param n_input The input field that we're blending.
   * \param n_output The output node that will contain the new field.
   */
  void execute(const BlendData &blend, const conduit::Node &n_input, conduit::Node &n_output) const
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
  void blendSingleComponent(const BlendData &blend, const conduit::Node &n_values, conduit::Node &n_output_values) const
  {
    const auto allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    const auto outputSize = blend.m_uniqueIndices.size();
    n_output_values.set_allocator(allocatorID);
    n_output_values.set(conduit::DataType(n_values.dtype().id(), outputSize);

    views::Node_to_ArrayView_same(n_values, n_output_values, [&](auto compView, auto outView)
    {
      using value_type = typename decltype(compView)::value_type;
      using accum_type = axom::mir::utilities::accumulation_traits<value_type>::type;

      axom::for_all<ExecSpace>(nbg, AXOM_LAMBDA(auto bgid)
      {
        // Original blendGroup index.
        const auto origBGIdx =  blend.m_uniqueIndicesView[bgid];
        const auto start = blend.m_blendGroupStartView[origBGIdx];
        const auto end   = start + blend.m_blendGroupSizesView[origBGIdx];

        accum_type blended = 0;
        for(int32 i = start; i < end; i++)
        {
          const auto index = blend.m_blendIdsView[i];
          const auto weight = blend.m_blendCoeffView[i];
          blended += static_cast<accum_type>(compView[index]) * weight;
        }
        outView[bgid] = static_cast<value_type>(blended);
      });
    });
  }
};

} // end namespace blueprint
} // end namespace utilities
} // end namespace mir
} // end namespace axom

#endif

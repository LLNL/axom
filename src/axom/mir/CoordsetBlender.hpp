// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_COORDSET_BLENDER_HPP_
#define AXOM_MIR_COORDSET_BLENDER_HPP_

#include "axom/core.hpp"
#include "axom/mir/FieldBlender.hpp" // for BlendData
#include "axom/mir/blueprint_utilities.hpp" // for cpp2conduit
#include "axom/primal/Point.hpp"

#include <conduit/conduit.hpp>
#include <conduit/conduit_blueprint_mesh_utils.hpp>

namespace axom
{
namespace mir
{
namespace utilities
{
namespace blueprint
{

/**
 * \accelerated
 * \class FieldBlender
 *
 * \tparam ExecSpace The execution space for the algorithm.
 * \tparam CSVType The coordset view type.
 *
 * \brief This class uses BlendData to generate a new blended field from an input coordset.
 *        The output coordset will be explicit.
 */
template <typename ExecSpace, typename CSVType>
class CoordsetBlender
{ 
public:
  using CoordsetViewType = CSVType;

  /**
   * \brief Create a new blended field from the \a n_input field and place it in \a n_output.
   *
   * \param blend The BlendData that will be used to make the new coordset.
   * \param n_input The input coordset that we're blending.
   * \param n_output The output node that will contain the new coordset.
   *
   * \note The coordset view must agree with the coordset in n_input. We pass both
   *       a view and the coordset node since the view may not be able to contain
   *       come coordset metadata and remain trivially copyable.
   */
  void execute(const BlendData &blend, const CoordsetViewType &view, const conduit::Node &n_input, conduit::Node &n_output) const
  {
    using value_type = typename decltype(view)::value_type;
    using PointType = typename decltype(view)::PointType;

    const auto allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    const auto axes = conduit::blueprint::mesh::utils::coordset::axes(n_input);
    const auto nComponents = axes.size();
    assert(PointType::NDIMS == nComponents);

    n_output.reset();
    n_output["type"] = "explicit";
    conduit::Node &n_values = n_output["values"];

    // Make output nodes using axis names from the input coordset. Make array views too.
    const auto outputSize = blend.m_uniqueIndices.size();
    axom::StackArray<axom::ArrayView<value_type>, PointType::NDIMS> compViews;
    for(size_t i = 0; i < nComponents; i++)
    {
      // Use the view's value_type
      n_output[axisName[i]].set(conduit::DataType(axom::mir::utilities::blueprint::cpp2conduit<value_type>::id, outputSize));
      auto *comp_data = static_cast<value_type *>(n_output.fetch_existing(axisName[i]).data_ptr());
      compViews[i] = axom::ArrayView<value_type>(comp_data, outputSize);
    }

    // Iterate over each blend group.
    axom::for_all<ExecSpace>(outputSize, AXOM_LAMBDA(auto bgid)
    {
      // Original blendGroup index.
      const auto origBGIdx =  blend.m_uniqueIndicesView[bgid];
      const auto start = blend.m_blendGroupStartView[origBGIdx];
      const auto end   = start + blend.m_blendGroupSizesView[origBGIdx];

      // Blend points for this blend group.
      PointType blend{};
      for(int32 i = start; i < end; i++)
      {
        const auto index = blend.m_blendIdsView[i];
        const auto weight = blend.m_blendCoeffView[i];

        blended += (view[index] * weight);
      }

      // Store the point into the Conduit component arrays.
      for(int comp = 0; comp < PointType::NDIMS; comp++)
      {
        compViews[comp][bgid] = blended[comp];
      }
    });
  }
};

} // end namespace blueprint
} // end namespace utilities
} // end namespace mir
} // end namespace axom

#endif

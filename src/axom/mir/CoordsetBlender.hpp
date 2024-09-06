// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_COORDSET_BLENDER_HPP_
#define AXOM_MIR_COORDSET_BLENDER_HPP_

#include "axom/core.hpp"
#include "axom/mir/FieldBlender.hpp"         // for BlendData
#include "axom/mir/blueprint_utilities.hpp"  // for cpp2conduit
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/slic.hpp"

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
 * \tparam SelectionPolicy The selection policy to use.
 *
 * \brief This class uses BlendData to generate a new blended field from an input coordset.
 *        The output coordset will be explicit.
 */
template <typename ExecSpace, typename CSVType, typename SelectionPolicy>
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
  void execute(const BlendData &blend,
               const CoordsetViewType &view,
               const conduit::Node &n_input,
               conduit::Node &n_output) const
  {
    using value_type = typename CoordsetViewType::value_type;
    using PointType = typename CoordsetViewType::PointType;
    using VectorType = axom::primal::Vector<value_type, PointType::DIMENSION>;

    const auto axes = conduit::blueprint::mesh::utils::coordset::axes(n_input);
    const auto nComponents = axes.size();
    SLIC_ASSERT(PointType::DIMENSION == nComponents);

    // Get the ID of a Conduit allocator that will allocate through Axom with device allocator allocatorID.
    utilities::blueprint::ConduitAllocateThroughAxom<ExecSpace> c2a;

    n_output.reset();
    n_output["type"] = "explicit";
    conduit::Node &n_values = n_output["values"];

    // Determine output size.
    const auto origSize = blend.m_originalIdsView.size();
    const auto blendSize = SelectionPolicy::size(blend);
    const auto outputSize = origSize + blendSize;

    // Make output nodes using axis names from the input coordset. Make array views too.
    axom::StackArray<axom::ArrayView<value_type>, PointType::DIMENSION> compViews;
    for(size_t i = 0; i < nComponents; i++)
    {
      // Allocate data in the Conduit node and make a view.
      conduit::Node &comp = n_values[axes[i]];
      comp.set_allocator(c2a.getConduitAllocatorID());
      comp.set(conduit::DataType(
        axom::mir::utilities::blueprint::cpp2conduit<value_type>::id,
        outputSize));
      auto *comp_data = static_cast<value_type *>(comp.data_ptr());
      compViews[i] = axom::ArrayView<value_type>(comp_data, outputSize);
    }

    const CoordsetViewType deviceView(view);
    const BlendData deviceBlend(blend);

    // Copy over some original values to the start of the array.
    axom::for_all<ExecSpace>(
      origSize,
      AXOM_LAMBDA(auto index) {
        const auto srcIndex = deviceBlend.m_originalIdsView[index];
        const auto pt = deviceView[srcIndex];

        // Store the point into the Conduit component arrays.
        for(int comp = 0; comp < PointType::DIMENSION; comp++)
        {
          compViews[comp][index] = pt[comp];
        }
    });

    // Append blended values to the end of the array.
    axom::for_all<ExecSpace>(
      blendSize,
      AXOM_LAMBDA(auto bgid) {
        // Get the blend group index we want.
        const auto selectedIndex =
          SelectionPolicy::selectedIndex(deviceBlend, bgid);
        const auto start = deviceBlend.m_blendGroupStartView[selectedIndex];
        const auto nValues = deviceBlend.m_blendGroupSizesView[selectedIndex];
        const auto destIndex = bgid + origSize;

        VectorType blended {};
        if(nValues == 1)
        {
          const auto index = deviceBlend.m_blendIdsView[start];
          blended = VectorType(deviceView[index]);
        }
        else
        {
          const auto end = start + deviceBlend.m_blendGroupSizesView[selectedIndex];

          // Blend points for this blend group.
          for(IndexType i = start; i < end; i++)
          {
            const auto index = deviceBlend.m_blendIdsView[i];
            const auto weight = deviceBlend.m_blendCoeffView[i];
            blended += (VectorType(deviceView[index]) * static_cast<value_type>(weight));
          }
        }

        // Store the point into the Conduit component arrays.
        for(int comp = 0; comp < PointType::DIMENSION; comp++)
        {
          compViews[comp][destIndex] = blended[comp];
        }
      });
  }
};

}  // end namespace blueprint
}  // end namespace utilities
}  // end namespace mir
}  // end namespace axom

#endif

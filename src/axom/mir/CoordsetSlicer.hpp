// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_COORDSET_SLICER_HPP_
#define AXOM_MIR_COORDSET_SLICER_HPP_

#include "axom/core.hpp"
#include "axom/mir.hpp"

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
 * \accelerated
 * \class CoordsetSlicer
 *
 * \brief This class uses SliceData to generate a new sliced coordset (pulling out specific points from input coordset).
 *
 * \tparam ExecSpace The execution space where the algorithm will run.
 *
 */
template <typename ExecSpace, typename CoordsetView>
class CoordsetSlicer
{
public:
  /// Constructor
  CoordsetSlicer(const CoordsetView &coordsetView) : m_coordsetView(coordsetView) { }

  /**
   * \brief Execute the slice on the \a n_input coordset and store the new sliced coordset in \a n_output.
   *
   * \param slice    The slice data that indicates how the coordset will be sliced.
   * \param n_input  A Conduit node containing the coordset to be sliced.
   * \param n_output A node that will contain the sliced coordset.
   *
   * \note We assume for now that n_input != n_output.
   */
  void execute(const SliceData &slice,
               const conduit::Node &n_input,
               conduit::Node &n_output)
  {
    using value_type = typename CoordsetView::value_type;
    using PointType = typename CoordsetView::PointType;
    namespace bputils = axom::mir::utilities::blueprint;

    const auto axes = conduit::blueprint::mesh::utils::coordset::axes(n_input);
    const auto nComponents = axes.size();
    SLIC_ASSERT(PointType::DIMENSION == nComponents);

    // Get the ID of a Conduit allocator that will allocate through Axom with device allocator allocatorID.
    bputils::ConduitAllocateThroughAxom<ExecSpace> c2a;

    n_output.reset();
    n_output["type"] = "explicit";
    conduit::Node &n_values = n_output["values"];

    // Determine output size.
    const auto outputSize = slice.m_indicesView;

    // Make output nodes using axis names from the input coordset. Make array views too.
    axom::StackArray<axom::ArrayView<value_type>, PointType::DIMENSION> compViews;
    for(size_t i = 0; i < nComponents; i++)
    {
      // Allocate data in the Conduit node and make a view.
      conduit::Node &comp = n_values[axes[i]];
      comp.set_allocator(c2a.getConduitAllocatorID());
      comp.set(conduit::DataType(
        bputils::cpp2conduit<value_type>::id,
        outputSize));
      auto *comp_data = static_cast<value_type *>(comp.data_ptr());
      compViews[i] = axom::ArrayView<value_type>(comp_data, outputSize);
    }

    // Select the nodes we want in the output.
    const CoordsetView deviceView(m_coordsetView);
    const auto deviceIndicesView = slice.m_indicesView;
    axom::for_all<ExecSpace>(
      outputSize,
      AXOM_LAMBDA(auto index) {
        const auto srcIndex = deviceIndicesView[index];
        const auto pt = deviceView[srcIndex];

        // Store the point into the Conduit component arrays.
        for(int comp = 0; comp < PointType::DIMENSION; comp++)
        {
          compViews[comp][index] = pt[comp];
        }
      });
  }
private:
  CoordsetView m_coordsetView;
};

}  // end namespace blueprint
}  // end namespace utilities
}  // end namespace mir
}  // end namespace axom

#endif

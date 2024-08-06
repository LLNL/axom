// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_FIELD_SLICER_HPP_
#define AXOM_MIR_FIELD_SLICER_HPP_

#include "axom/core.hpp"
#include "axom/mir/views/NodeArrayView.hpp"
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
 * \brief Contains the indices to be sliced out of a Blueprint field.
 */
struct SliceData
{
  axom::ArrayView<IndexType> m_indicesView;
};

/**
 * \accelerated
 * \class FieldSlicer
 *
 * \brief This class uses SliceData to generate a new sliced field from an input field.
 *
 * \tparam ExecSpace The execution space where the algorithm will run.
 * \tparam IndexingPolicy A class that provides operator[] that can transform zone indices.
 *
 */
template <typename ExecSpace, typename IndexingPolicy = DirectIndexing>
class FieldSlicer
{
public:
  /// Constructor
  FieldSlicer() : m_indexing()
  {
  }

  /**
   * \brief Constructor
   * \param indexing An object used to transform node indices.
   */
  FieldSlicer(const IndexingPolicy &indexing) : m_indexing(indexing)
  {
  }

  /**
   * \brief Execute the slice on the \a n_input field and store the new sliced field in \a n_output.
   *
   * \param slice    The slice data that indicates how the field will be sliced.
   * \param n_input  A Conduit node containing the field to be sliced.
   * \param n_output A node that will contain the sliced field.
   *
   * \note We assume for now that n_input != n_output.
   */
  void execute(const SliceData &slice,
               const conduit::Node &n_input,
               conduit::Node &n_output)
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
        sliceSingleComponent(slice, n_input, n_comp, n_out_comp);
      }
    }
    else
    {
      sliceSingleComponent(slice, n_input, n_input_values, n_output_values);
    }
  }

private:
  /**
   * \brief Slice data for a single field component.
   *
   * \param slice The SliceData that will be used to make the new field.
   * \param n_input The node that contains the field.
   * \param n_values The input values that we're slicing.
   * \param n_output_values The output node that will contain the new field.
   */
  void sliceSingleComponent(const SliceData &slice,
                            const conduit::Node &n_input,
                            const conduit::Node &n_values,
                            conduit::Node &n_output_values) const
  {
    const auto outputSize = slice.m_indicesView.size();

    // Allocate Conduit data through Axom.
    utilities::blueprint::ConduitAllocateThroughAxom<ExecSpace> c2a;
    n_output_values.set_allocator(c2a.getConduitAllocatorID());
    n_output_values.set(conduit::DataType(n_values.dtype().id(), outputSize));

    views::Node_to_ArrayView_same(
      n_values,
      n_output_values,
      [&](auto valuesView, auto outputView) {

        // Let the indexing object update itself from the node strides/offsets.
        IndexingPolicy deviceIndexing(m_indexing);
        deviceIndexing.update(n_input);

        SliceData deviceSlice(slice);
        axom::for_all<ExecSpace>(
          outputSize,
          AXOM_LAMBDA(auto index) {
            const auto zoneIndex = deviceSlice.m_indicesView[index];
            const auto transformedIndex = deviceIndexing[zoneIndex];
            outputView[index] = valuesView[transformedIndex];
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

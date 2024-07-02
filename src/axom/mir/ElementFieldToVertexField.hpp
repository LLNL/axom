// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_ELEMENT_FIELD_TO_VERTEX_FIELD_HPP_
#define AXOM_MIR_ELEMENT_FIELD_TO_VERTEX_FIELD_HPP_

#include "axom/core.hpp"
#include "axom/mir/utilities.hpp"

#include <conduit/conduit.hpp>

#include <RAJA/RAJA.hpp>

namespace axom
{
namespace mir
{
namespace utilities
{

/**
 * \brief Convert an element field to a vertex field.
 */
template <typename ExecSpace>
class ElementFieldToVertexField
{
public:

  /**
   * \brief Convert an element field to a vertex field, storing the new field in a new Conduit node.
   *
   * \param field       The input field.
   * \param relation    The node that contains an o2m relation with nodes to zones.
   * \param outField[out] The node that will contain the new field.
   */
  void execute(const conduit::Node &field, const conduit::Node &relation, conduit::Node &outField);
};

template <typename ExecSpace>
void
ElementFieldToVertexField<ExecSpace>::execute(const conduit::Node &field, const conduit::Node &relation, conduit::Node &outField)
{
  using loop_policy = axom::execution_space<ExecSpace>::loop_policy;
  const std::string association = field.fetch_existing("association").as_string();
  const auto allocatorID = axom::execution_space<ExecSpace>::allocatorID();

  outField["association"] = "vertex";
  outField["topology"] = field["topology"];

  auto handleComponent = [](const conduit::Node &relation, const conduit::Node &n_comp, conduit::Node &n_out, int allocatorID)
  {
    n_out.set_allocator(allocatorID);
    n_out.set(n_comp.dtype().id(), n_comp.dtype().number_of_elements());

    const conduit::Node &n_zones = relation["zones"];
    const conduit::Node &n_sizes = relation["sizes"];
    const conduit::Node &n_offsets = relation["offsets"];
    views::IndexNode_to_ArrayView_same(n_zones, n_sizes, n_offsets, [&](auto zonesView, auto sizesView, auto offsetsView)
    {
      views::Node_to_ArrayView_same(n_comp, n_out, [&](auto compView, auto outView)
      {
        using Precision = typename compView::value_type;
        const auto nnodes = sizesView.size();
        axom::for_all<ExecSpace>(nnodes, AXOM_LAMBDA(int nodeIndex)
        {
          const auto n = sizesView[nodeIndex];
          const auto offset = offsetsView[nodeIndex];

          double sum = 0.;
          for(int i = 0; i < n; i++)
          {
            const auto zid = zonesView[offset + i];
            sum += compView[zid];
          }

          outView[nodeIndex] = static_cast<Precision>(sum / n);
        });
      });
    });
  };

  if(association == "element")
  {
    const conduit::Node &n_zones = relation["zones"];
    const conduit::Node &n_sizes = relation["sizes"];
    const conduit::Node &n_offsets = relation["offsets"];
    views::IndexNode_to_ArrayView_same(n_zones, n_sizes, n_offsets, [&](auto zonesView, auto sizesView, auto offsetsView)
    {
      const auto nnodes = sizesView.size();
      const conduit::Node &n_values = field["values"];
      if(n_values.number_of_children() > 0)
      {
        for(conduit::index_t c = 0; c < n_values.number_of_children(); c++)
        {
          const conduit::Node &n_comp = n_values[c];
          handleComponent(relation, n_comp, outField["values"][n_comp.name()], allocatorID);
        }
      }
      else
      {
        handleComponent(relation, n_values, outField["values"], allocatorID);
      }
    });
  }
  else if(association == "vertex")
  {
    outField["values"].set_external(field["values"]);
  }
}

} // end namespace utilities
} // end namespace mir
} // end namespace axom

#endif

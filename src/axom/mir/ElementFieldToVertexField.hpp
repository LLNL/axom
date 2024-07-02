// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_RECENTER_FIELD_HPP_
#define AXOM_MIR_RECENTER_FIELD_HPP_

#include "axom/core.hpp"
#include "axom/mir/views/NodeToArrayView.hpp"

#include <conduit/conduit.hpp>
#include <conduit/conduit_blueprint.hpp>

namespace axom
{
namespace mir
{

/**
 * \brief This struct contains the type that should be used to accumulate values of type T.
 */
template <typename T>
struct accumulate_traits
{
  using value_type = double;
};

template <>
struct accumulate_traits<float>
{
  using value_type = float;
};

/**
 * \brief Convert a field with one association type to a field of another association type using an o2mrelation.
 */
template <typename ExecSpace>
class RecenterField
{
public:

  /**
   * \brief Convert the input field to a different association type using the o2mrelation and store the new field in the output field.
   *
   * \param field       The input field.
   * \param relation    The node that contains an o2mrelation with nodes to zones.
   * \param outField[out] The node that will contain the new field.
   */
  static void execute(const conduit::Node &field, const conduit::Node &relation, conduit::Node &outField);
};

template <typename ExecSpace>
void
RecenterField<ExecSpace>::execute(const conduit::Node &field, const conduit::Node &relation, conduit::Node &outField)
{
  auto handleComponent = [](const conduit::Node &relation, const conduit::Node &n_comp, conduit::Node &n_out, int allocatorID)
  {
    // Get the data field for the o2m relation.
    const auto data_paths = conduit::blueprint::o2mrelation::data_paths(relation);

    // Use the o2mrelation to average data from n_comp to the n_out.
    const conduit::Node &n_relvalues = relation[data_paths[0]];
    const conduit::Node &n_sizes = relation["sizes"];
    const conduit::Node &n_offsets = relation["offsets"];
    views::IndexNode_to_ArrayView_same(n_relvalues, n_sizes, n_offsets, [&](auto relView, auto sizesView, auto offsetsView)
    {
      const auto relSize = sizesView.size();

      // Allocate data for n_out (same type as n_comp).
      n_out.set_allocator(allocatorID);
      n_out.set(n_comp.dtype().id(), relSize);

      views::Node_to_ArrayView_same(n_comp, n_out, [&](auto compView, auto outView)
      {
        using Precision = typename decltype(compView)::value_type;
        using AccumType = accumulate_traits<Precision>::value_type;
        axom::for_all<ExecSpace>(relSize, AXOM_LAMBDA(int relIndex)
        {
          const auto n = sizesView[relIndex];
          const auto offset = offsetsView[relIndex];

          AccumType sum = 0;
          for(int i = 0; i < n; i++)
          {
            const auto id = relView[offset + i];
            sum += static_cast<AccumType>(compView[id]);
          }

          outView[relIndex] = static_cast<Precision>(sum / n);
        });
      });
    });
  };

  const std::string association = field.fetch_existing("association").as_string();
  const auto allocatorID = axom::execution_space<ExecSpace>::allocatorID();

  // Assume that we're flipping the association.
  outField["association"] = (association == "element") ? "vertex" : "element";
  outField["topology"] = field["topology"];

  // Make output values.
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
}

} // end namespace mir
} // end namespace axom

#endif

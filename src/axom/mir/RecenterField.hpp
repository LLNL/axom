// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_RECENTER_FIELD_HPP_
#define AXOM_MIR_RECENTER_FIELD_HPP_

#include "axom/core.hpp"
#include "axom/mir/views/NodeArrayView.hpp"
#include "axom/mir/utilities.hpp"

#include <conduit/conduit.hpp>
#include <conduit/conduit_blueprint.hpp>

namespace axom
{
namespace mir
{
namespace utilities
{
namespace blueprint
{
/**
 * \brief Convert a field with one association type to a field of another association type using an o2mrelation.
 *
 * \tparam ExecSpace The execution space where the algorithm runs.
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
  void execute(const conduit::Node &field,
               const conduit::Node &relation,
               conduit::Node &outField) const;

private:
  /**
   * \brief Recenter a single field component.
   *
   * \param relation    The node that contains an o2mrelation with nodes to zones.
   * \param n_comp      The input component.
   * \param n_out[out] The node that will contain the new field.
   */
  void recenterSingleComponent(const conduit::Node &n_comp,
                               const conduit::Node &relation,
                               conduit::Node &n_out) const;
};

template <typename ExecSpace>
void RecenterField<ExecSpace>::execute(const conduit::Node &field,
                                       const conduit::Node &relation,
                                       conduit::Node &outField) const
{
  const std::string association = field.fetch_existing("association").as_string();

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
      recenterSingleComponent(n_comp, relation, outField["values"][n_comp.name()]);
    }
  }
  else
  {
    recenterSingleComponent(n_values, relation, outField["values"]);
  }
}

template <typename ExecSpace>
void RecenterField<ExecSpace>::recenterSingleComponent(
  const conduit::Node &n_comp,
  const conduit::Node &relation,
  conduit::Node &n_out) const
{
  // Get the data field for the o2m relation.
  const auto data_paths = conduit::blueprint::o2mrelation::data_paths(relation);

  // Use the o2mrelation to average data from n_comp to the n_out.
  const conduit::Node &n_relvalues = relation[data_paths[0]];
  const conduit::Node &n_sizes = relation["sizes"];
  const conduit::Node &n_offsets = relation["offsets"];
  views::IndexNode_to_ArrayView_same(
    n_relvalues,
    n_sizes,
    n_offsets,
    [&](auto relView, auto sizesView, auto offsetsView) {

      // Allocate Conduit data through Axom.
      const auto relSize = sizesView.size();
      utilities::blueprint::ConduitAllocateThroughAxom<ExecSpace> c2a;
      n_out.set_allocator(c2a.getConduitAllocatorID());
      n_out.set(conduit::DataType(n_comp.dtype().id(), relSize));

      views::Node_to_ArrayView_same(n_comp, n_out, [&](auto compView, auto outView) {
        using Precision = typename decltype(compView)::value_type;
        using AccumType =
          typename axom::mir::utilities::accumulation_traits<Precision>::value_type;
        axom::for_all<ExecSpace>(
          relSize,
          AXOM_LAMBDA(auto relIndex) {
            const auto n = static_cast<axom::IndexType>(sizesView[relIndex]);
            const auto offset = offsetsView[relIndex];

            AccumType sum {};
            for(axom::IndexType i = 0; i < n; i++)
            {
              const auto id = relView[offset + i];
              sum += static_cast<AccumType>(compView[id]);
            }

            outView[relIndex] = static_cast<Precision>(sum / n);
          });
      });
    });
}

}  // end namespace blueprint
}  // end namespace utilities
}  // end namespace mir
}  // end namespace axom

#endif

// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_DISPATCH_UNIFORM_TOPOLOGY_HPP_
#define AXOM_MIR_DISPATCH_UNIFORM_TOPOLOGY_HPP_

#include "axom/mir/views/StructuredTopologyView.hpp"
#include "axom/mir/views/dispatch_utilities.hpp"
#include <conduit/conduit.hpp>

namespace axom
{
namespace mir
{
namespace views
{
/**
 * \brief Creates a topology view compatible with uniform topologies and passes that view to the supplied function.
 *
 * \tparam FuncType            The function/lambda type to invoke on the view.
 * \tparam SelectedDimensions  An integer whose bits indicate which dimensions are set.
 *
 * \param coordset The coordset node that contains the topology dimensions.
 * \param func     The function to invoke using the view.
 */
template <int SelectedDimensions = select_dimensions(1, 2, 3), typename FuncType>
void dispatch_uniform_topology(const conduit::Node &AXOM_UNUSED_PARAM(topo),
                               const conduit::Node &coordset,
                               FuncType &&func)
{
  const conduit::Node &n_dims = coordset["dims"];
  switch(n_dims.dtype().number_of_elements())
  {
  default:
  case 3:
    if constexpr(dimension_selected(SelectedDimensions, 3))
    {
      axom::StackArray<axom::IndexType, 3> zoneDims;
      zoneDims[0] = n_dims.as_int_accessor()[0] - 1;
      zoneDims[1] = n_dims.as_int_accessor()[1] - 1;
      zoneDims[2] = n_dims.as_int_accessor()[2] - 1;
      views::StructuredTopologyView<views::StructuredIndexing<axom::IndexType, 3>>
        topoView(zoneDims);
      const std::string shape("hex");
      func(shape, topoView);
    }
    break;
  case 2:
    if constexpr(dimension_selected(SelectedDimensions, 2))
    {
      axom::StackArray<axom::IndexType, 2> zoneDims;
      zoneDims[0] = n_dims.as_int_accessor()[0] - 1;
      zoneDims[1] = n_dims.as_int_accessor()[1] - 1;
      views::StructuredTopologyView<views::StructuredIndexing<axom::IndexType, 2>>
        topoView(zoneDims);
      const std::string shape("quad");
      func(shape, topoView);
    }
    break;
  case 1:
    if constexpr(dimension_selected(SelectedDimensions, 3))
    {
      axom::StackArray<axom::IndexType, 1> zoneDims;
      zoneDims[0] = n_dims.as_int_accessor()[0] - 1;
      views::StructuredTopologyView<views::StructuredIndexing<axom::IndexType, 1>>
        topoView(zoneDims);
      const std::string shape("line");
      func(shape, topoView);
    }
    break;
  }
}

}  // end namespace views
}  // end namespace mir
}  // end namespace axom

#endif

// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_DISPATCH_RECTILINEAR_TOPOLOGY_HPP_
#define AXOM_MIR_DISPATCH_RECTILINEAR_TOPOLOGY_HPP_

#include "axom/mir/views/StructuredTopologyView.hpp"
#include "axom/mir/views/dispatch_utilities.hpp"
#include <conduit/conduit.hpp>
#include <conduit/conduit_blueprint_mesh_utils.hpp>

namespace axom
{
namespace mir
{
namespace views
{

/**
 * \brief Creates a topology view compatible with rectilinear topologies and passes that view to the supplied function.
 *
 * \tparam FuncType The function/lambda type to invoke on the view.
 * \tparam SelectedDimensions  An integer whose bits indicate which dimensions are set.
 *
 * \param topo     The node that contains the rectilinear topology.
 * \param coordset The coordset node that contains the topology dimensions.
 * \param func     The function to invoke using the view.
 */
template <int SelectedDimensions = select_dimensions(1,2,3), typename FuncType>
void dispatch_rectilinear_topology(const conduit::Node &topo, const conduit::Node &coordset, FuncType &&func)
{
  const auto axes = conduit::blueprint::mesh::utils::coordset::axes(coordset);
  switch(axes.size())
  {
  case 3:
    if constexpr (dimension_selected(SelectedDimensions, 3))
    {
      axom::StackArray<axom::IndexType, 3> dims;
      dims[0] = coordset.fetch_existing(axes[0]).dtype().number_of_elements();
      dims[1] = coordset.fetch_existing(axes[1]).dtype().number_of_elements();
      dims[2] = coordset.fetch_existing(axes[2]).dtype().number_of_elements();
      views::StructuredTopologyView<axom::IndexType, 3> topoView(dims);
      const std::string shape("hex");
      func(shape, topoView);
    }
    break;
  case 2:
    if constexpr (dimension_selected(SelectedDimensions, 2))
    {
      axom::StackArray<axom::IndexType, 2> dims;
      dims[0] = coordset.fetch_existing(axes[0]).dtype().number_of_elements();
      dims[1] = coordset.fetch_existing(axes[1]).dtype().number_of_elements();
      views::StructuredTopologyView<axom::IndexType, 2> topoView(dims);
      const std::string shape("quad");
      func(shape, topoView);
    }
    break;
  case 1:
    if constexpr (dimension_selected(SelectedDimensions, 1))
    {
      axom::StackArray<axom::IndexType, 1> dims;
      dims[0] = coordset.fetch_existing(axes[0]).dtype().number_of_elements();
      views::StructuredTopologyView<axom::IndexType, 1> topoView(dims);
      const std::string shape("line");
      func(shape, topoView);
    }
    break;
  default:
    break;
  }
}

} // end namespace views
} // end namespace mir
} // end namespace axom

#endif

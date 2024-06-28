// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_DISPATCH_RECTILINEAR_TOPOLOGY_HPP_
#define AXOM_MIR_DISPATCH_RECTILINEAR_TOPOLOGY_HPP_

#include "axom/mir/views/StructuredTopologyView.hpp"
#include <conduit/conduit.hpp>

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
 *
 * \param topo     The node that contains the rectilinear topology.
 * \param coordset The coordset node that contains the topology dimensions.
 * \param func     The function to invoke using the view.
 */
template <typename FuncType>
void dispatch_rectilinear_topology(const conduit::Node &topo, const conduit::Node &coordset, FuncType &&func)
{
  const auto axes = conduit::blueprint::mesh::utils::coordset::axes(coordset);
  if(axes.size() == 3)
  {
    axom::StackArray<axom::IndexType, 3> dims;
    dims[0] = coordset.fetch_existing(axes[0]).dtype().number_of_elements();
    dims[1] = coordset.fetch_existing(axes[1]).dtype().number_of_elements();
    dims[2] = coordset.fetch_existing(axes[2]).dtype().number_of_elements();
    views::StructuredTopologyView<axom::IndexType, 3> topoView(dims);
    const std::string shape("hex");
    func(shape, topoView);
  }
  else if(axes.size() == 2)
  {
    axom::StackArray<axom::IndexType, 2> dims;
    dims[0] = coordset.fetch_existing(axes[0]).dtype().number_of_elements();
    dims[1] = coordset.fetch_existing(axes[1]).dtype().number_of_elements();
    views::StructuredTopologyView<axom::IndexType, 2> topoView(dims);
    const std::string shape("quad");
    func(shape, topoView);
  }
  else if(axes.size() == 1)
  {
    axom::StackArray<axom::IndexType, 1> dims;
    dims[0] = coordset.fetch_existing(axes[0]).dtype().number_of_elements();
    views::StructuredTopologyView<axom::IndexType, 1> topoView(dims);
    const std::string shape("line");
    func(shape, topoView);
  }
}

} // end namespace views
} // end namespace mir
} // end namespace axom

#endif

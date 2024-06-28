// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_DISPATCH_UNIFORM_TOPOLOGY_HPP_
#define AXOM_MIR_DISPATCH_UNIFORM_TOPOLOGY_HPP_

#include "axom/mir/views/StructuredTopologyView.hpp"
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
 * \tparam FuncType The function/lambda type to invoke on the view.
 *
 * \param coordset The coordset node that contains the topology dimensions.
 * \param func     The function to invoke using the view.
 */
template <typename FuncType>
void dispatch_uniform_topology(const conduit::Node & AXOM_UNUSED_PARAM(topo), const conduit::Node &coordset, FuncType &&func)
{
  const conduit::Node &n_dims = coordset["dims"];
  const conduit::index_t ndims = n_dims.dtype().number_of_elements();
  if(ndims == 3)
  {
    axom::StackArray<axom::IndexType, 3> dims;
    dims[0] = n_dims.as_int_accessor()[0];
    dims[1] = n_dims.as_int_accessor()[1];
    dims[2] = n_dims.as_int_accessor()[2];
    views::StructuredTopologyView<axom::IndexType, 3> topoView(dims);
    const std::string shape("hex");
    func(shape, topoView);
  }
  else if(ndims == 2)
  {
    axom::StackArray<axom::IndexType, 2> dims;
    dims[0] = n_dims.as_int_accessor()[0];
    dims[1] = n_dims.as_int_accessor()[1];
    views::StructuredTopologyView<axom::IndexType, 2> topoView(dims);
    const std::string shape("quad");
    func(shape, topoView);
  }
  else if(ndims == 1)
  {
    axom::StackArray<axom::IndexType, 1> dims;
    dims[0] = n_dims.as_int_accessor()[0];
    views::StructuredTopologyView<axom::IndexType, 1> topoView(dims);
    const std::string shape("line");
    func(shape, topoView);
  }
}

} // end namespace views
} // end namespace mir
} // end namespace axom

#endif

// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_DISPATCH_TOPOLOGY_HPP_
#define AXOM_MIR_DISPATCH_TOPOLOGY_HPP_

#include "axom/mir/views/StructuredTopologyView.hpp"
#include "axom/mir/views/dispatch_uniform_topology.hpp"
#include "axom/mir/views/dispatch_rectilinear_topology.hpp"
#include "axom/mir/views/dispatch_structured_topology.hpp"
#include "axom/mir/views/dispatch_unstructured_topology.hpp"

#include <conduit/conduit.hpp>

namespace axom
{
namespace mir
{
namespace views
{

/**
 * \brief Creates a topology view and passes that view to the supplied function.
 *
 * \tparam FuncType The function/lambda type to invoke on the view.
 * \tparam SelectedDimensions  An integer whose bits indicate which dimensions are set. dimension
 *
 * \param topo     The node that contains the rectilinear topology.
 * \param coordset The coordset node that contains the topology dimensions.
 * \param func     The function to invoke using the view. It should accept a string with the shape name and an auto parameter for the view.
 */
template <int SelectedDimensions = select_dimensions(1,2,3), typename FuncType>
void dispatch_topology(const conduit::Node &topo, const conduit::Node &coordset, FuncType &&func)
{
  const auto type = topo.fetch_existing("type").as_string();

  if(type == "uniform")
    dispatch_uniform_topology<SelectedDimensions>(topo, coordset, func);
  else if(type == "rectilinear")
    dispatch_rectilinear_topology<SelectedDimensions>(topo, coordset, func);
  else if(type == "structured")
    dispatch_structured_topology<SelectedDimensions>(topo, func);
  else if(type == "unstructured")
    dispatch_unstructured_topology(topo, func);
}

} // end namespace views
} // end namespace mir
} // end namespace axom

#endif

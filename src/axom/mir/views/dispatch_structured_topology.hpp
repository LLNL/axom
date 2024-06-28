// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_DISPATCH_STRUCTURED_TOPOLOGY_HPP_
#define AXOM_MIR_DISPATCH_STRUCTURED_TOPOLOGY_HPP_

#include "axom/mir/views/StructuredTopologyView.hpp"
#include "axom/mir/views/dispatch_uniform_topology.hpp"
#include "axom/mir/views/dispatch_rectilinear_topology.hpp"

#include <conduit/conduit.hpp>

namespace axom
{
namespace mir
{
namespace views
{

/**
 * \brief Creates a topology view compatible with structured topologies and passes that view to the supplied function.
 *
 * \tparam FuncType The function/lambda type to invoke on the view.
 *
 * \param topo     The node that contains the rectilinear topology.
 * \param coordset The coordset node that contains the topology dimensions.
 * \param func     The function to invoke using the view. It should accept a string with the shape name and an auto parameter for the view.
 */
template <typename FuncType>
void dispatch_structured_topology(const conduit::Node &topo, FuncType &&func)
{
  // TODO: add strided structured variant
  //StructuredTopologyView<StructuredIndexing<axom::IndexType, 3>> topoView(StructuredIndexing<axom::IndexType, 3>(dims));
  //StructuredTopologyView<StridedStructuredIndexing<axom::IndexType, 3>> topoView(StridedStructuredIndexing<axom::IndexType, 3>(dims, offset, stride));


  if(topo.has_path("elements/dims/k"))
  {
    axom::StackArray<axom::IndexType, 3> dims;
    dims[0] = topo.fetch_existing("elements/dims/i").as_int();
    dims[1] = topo.fetch_existing("elements/dims/j").as_int();
    dims[2] = topo.fetch_existing("elements/dims/k").as_int();
    views::StructuredTopologyView<axom::IndexType, 3> topoView(dims);
    const std::string shape("hex");
    func(shape, topoView);
  }
  else if(topo.has_path("elements/dims/j"))
  {
    axom::StackArray<axom::IndexType, 2> dims;
    dims[0] = topo.fetch_existing("elements/dims/i").as_int();
    dims[1] = topo.fetch_existing("elements/dims/j").as_int();
    views::StructuredTopologyView<axom::IndexType, 2> topoView(dims);
    const std::string shape("quad");
    func(shape, topoView);
  }
  else
  {
    axom::StackArray<axom::IndexType, 1> dims;
    dims[0] = topo.fetch_existing("elements/dims/i").as_int();
    views::StructuredTopologyView<axom::IndexType, 1> topoView(dims);
    const std::string shape("line");
    func(shape, topoView);
  }
}

/**
 * \brief Creates a topology view compatible with various logically "structured" topologies (uniform, rectilinear, structured) and passes that view to the supplied function.
 *
 * \param topo     The node that contains the topology.
 * \param coordset The coordset node that contains the topology dimensions.
 * \param func     The function to invoke using the view. It should accept a string with the shape name and an auto parameter for the view.

 */
template <typename FuncType>
void dispatch_structured_topologies(const conduit::Node &topo, const conduit::Node &coordset, FuncType &&func)
{
  const std::string type = topo["type"].as_string();
  if(type == "uniform")
    dispatch_uniform_topology(topo, coordset, func);
  else if(type == "rectilinear")
    dispatch_rectilinear_topology(topo, coordset, func);
  else if(type == "structured")
    dispatch_structured_topology(topo, coordset, func);
}

} // end namespace views
} // end namespace mir
} // end namespace axom

#endif

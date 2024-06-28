// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_DISPATCH_STRUCTURED_TOPOLOGY_HPP_
#define AXOM_MIR_DISPATCH_STRUCTURED_TOPOLOGY_HPP_

#include "axom/mir/views/StructuredTopologyView.hpp"
#include "axom/mir/views/dispatch_utilities.hpp"
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
 * \tparam SelectedDimensions  An integer whose bits indicate which dimensions are set. dimension
 *
 * \param topo     The node that contains the rectilinear topology.
 * \param coordset The coordset node that contains the topology dimensions.
 * \param func     The function to invoke using the view. It should accept a string with the shape name and an auto parameter for the view.
 */
template <int SelectedDimensions = select_dimensions(1,2,3), typename FuncType>
void dispatch_structured_topology(const conduit::Node &topo, FuncType &&func)
{
  // TODO: add strided structured variant
  //StructuredTopologyView<StructuredIndexing<axom::IndexType, 3>> topoView(StructuredIndexing<axom::IndexType, 3>(dims));
  //StructuredTopologyView<StridedStructuredIndexing<axom::IndexType, 3>> topoView(StridedStructuredIndexing<axom::IndexType, 3>(dims, offset, stride));
  int ndims = 1;
  ndims += topo.has_path("elements/dims/j") ? 1 : 0;
  ndims += topo.has_path("elements/dims/k") ? 1 : 0;
  switch(ndims)
  {
  case 3:
    if constexpr (dimension_selected(SelectedDimensions, 3))
    {
      axom::StackArray<axom::IndexType, 3> dims;
      dims[0] = topo.fetch_existing("elements/dims/i").as_int();
      dims[1] = topo.fetch_existing("elements/dims/j").as_int();
      dims[2] = topo.fetch_existing("elements/dims/k").as_int();
      views::StructuredTopologyView<axom::IndexType, 3> topoView(dims);
      const std::string shape("hex");
      func(shape, topoView);
    }
    break;
  case 2:
    if constexpr (dimension_selected(SelectedDimensions, 2))
    {
      axom::StackArray<axom::IndexType, 2> dims;
      dims[0] = topo.fetch_existing("elements/dims/i").as_int();
      dims[1] = topo.fetch_existing("elements/dims/j").as_int();
      views::StructuredTopologyView<axom::IndexType, 2> topoView(dims);
      const std::string shape("quad");
      func(shape, topoView);
    }
    break;
  case 1:
    if constexpr (dimension_selected(SelectedDimensions, 1))
    {
      axom::StackArray<axom::IndexType, 1> dims;
      dims[0] = topo.fetch_existing("elements/dims/i").as_int();
      views::StructuredTopologyView<axom::IndexType, 1> topoView(dims);
      const std::string shape("line");
      func(shape, topoView);
    }
  }
}

/**
 * \brief Creates a topology view compatible with various logically "structured" topologies (uniform, rectilinear, structured) and passes that view to the supplied function.
 *
 * \tparam FuncType The function/lambda type to invoke on the view.
 * \tparam SelectedDimensions  An integer whose bits indicate which dimensions are set. dimension
 *
 * \param topo     The node that contains the topology.
 * \param coordset The coordset node that contains the topology dimensions.
 * \param func     The function to invoke using the view. It should accept a string with the shape name and an auto parameter for the view.

 */
template <int SelectedDimensions = select_dimensions(1,2,3), typename FuncType>
void dispatch_structured_topologies(const conduit::Node &topo, const conduit::Node &coordset, FuncType &&func)
{
  const std::string type = topo["type"].as_string();
  if(type == "uniform")
    dispatch_uniform_topology<SelectedDimensions>(topo, coordset, func);
  else if(type == "rectilinear")
    dispatch_rectilinear_topology<SelectedDimensions>(topo, coordset, func);
  else if(type == "structured")
    dispatch_structured_topology<SelectedDimensions>(topo, func);
}

} // end namespace views
} // end namespace mir
} // end namespace axom

#endif

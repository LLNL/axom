// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_DISPATCH_UNIFORM_TOPOLOGY_HPP_
#define AXOM_MIR_DISPATCH_UNIFORM_TOPOLOGY_HPP_

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
 * \brief Base template for uniform topology creation
 */
template <int NDIMS>
struct make_uniform
{ };

/**
 * \brief Create a 3D structured topology view with normal structured indexing.
 */
template <>
struct make_uniform<3>
{
  using Indexing = views::StructuredIndexing<axom::IndexType, 3>;
  using LogicalIndex = typename Indexing::LogicalIndex;
  using TopoView = views::StructuredTopologyView<Indexing>;

  /**
   * \brief Create the indexing and initialize it from the topology.
   * \param topo The node containing the topology.
   * \return The indexing.
   */
  static Indexing indexing(const conduit::Node &topo)
  {
    const conduit::Node *coordset =
      conduit::blueprint::mesh::utils::find_reference_node(topo, "coordset");
    SLIC_ASSERT(coordset != nullptr);
    const conduit::Node &n_dims = coordset->fetch_existing("dims");
    LogicalIndex zoneDims;
    zoneDims[0] = n_dims[0].to_index_t() - 1;
    zoneDims[1] = n_dims[1].to_index_t() - 1;
    zoneDims[2] = n_dims[2].to_index_t() - 1;
    return Indexing(zoneDims);
  }

  /**
   * \brief Create the topology view and initialize it from the topology.
   * \param topo The node containing the topology.
   * \return The topology view.
   */
  static TopoView view(const conduit::Node &topo)
  {
    return TopoView(indexing(topo));
  }
};

/**
 * \brief Create a 2D structured topology view with normal structured indexing.
 */
template <>
struct make_uniform<2>
{
  using Indexing = views::StructuredIndexing<axom::IndexType, 2>;
  using LogicalIndex = typename Indexing::LogicalIndex;
  using TopoView = views::StructuredTopologyView<Indexing>;

  /**
   * \brief Create the indexing and initialize it from the topology.
   * \param topo The node containing the topology.
   * \return The indexing.
   */
  static Indexing indexing(const conduit::Node &topo)
  {
    const conduit::Node *coordset =
      conduit::blueprint::mesh::utils::find_reference_node(topo, "coordset");
    SLIC_ASSERT(coordset != nullptr);
    const conduit::Node &n_dims = coordset->fetch_existing("dims");
    LogicalIndex zoneDims;
    zoneDims[0] = n_dims[0].to_index_t() - 1;
    zoneDims[1] = n_dims[1].to_index_t() - 1;
    return Indexing(zoneDims);
  }

  /**
   * \brief Create the topology view and initialize it from the topology.
   * \param topo The node containing the topology.
   * \return The topology view.
   */
  static TopoView view(const conduit::Node &topo)
  {
    return TopoView(indexing(topo));
  }
};

/**
 * \brief Create a 1D structured topology view with normal structured indexing.
 */
template <>
struct make_uniform<1>
{
  using Indexing = views::StructuredIndexing<axom::IndexType, 1>;
  using LogicalIndex = typename Indexing::LogicalIndex;
  using TopoView = views::StructuredTopologyView<Indexing>;

  /**
   * \brief Create the indexing and initialize it from the topology.
   * \param topo The node containing the topology.
   * \return The indexing.
   */
  static Indexing indexing(const conduit::Node &topo)
  {
    const conduit::Node *coordset =
      conduit::blueprint::mesh::utils::find_reference_node(topo, "coordset");
    SLIC_ASSERT(coordset != nullptr);
    const conduit::Node &n_dims = coordset->fetch_existing("dims");
    LogicalIndex zoneDims;
    zoneDims[0] = n_dims[0].to_index_t() - 1;
    return Indexing(zoneDims);
  }

  /**
   * \brief Create the topology view and initialize it from the topology.
   * \param topo The node containing the topology.
   * \return The topology view.
   */
  static TopoView view(const conduit::Node &topo)
  {
    return TopoView(indexing(topo));
  }
};

/**
 * \brief Creates a topology view compatible with uniform topologies and passes that view to the supplied function.
 *
 * \tparam FuncType            The function/lambda type to invoke on the view.
 * \tparam SelectedDimensions  An integer whose bits indicate which dimensions are set.
 *
 * \param topo The node that contains the topology.
 * \param func The function to invoke using the view.
 */
template <int SelectedDimensions = select_dimensions(1, 2, 3), typename FuncType>
void dispatch_uniform_topology(const conduit::Node &topo, FuncType &&func)
{
  const conduit::Node *coordset =
    conduit::blueprint::mesh::utils::find_reference_node(topo, "coordset");
  SLIC_ASSERT(coordset != nullptr);
  const conduit::Node &n_dims = coordset->fetch_existing("dims");
  switch(n_dims.dtype().number_of_elements())
  {
  default:
  case 3:
    if constexpr(dimension_selected(SelectedDimensions, 3))
    {
      auto topoView = make_uniform<3>::view(topo);
      const std::string shape("hex");
      func(shape, topoView);
    }
    break;
  case 2:
    if constexpr(dimension_selected(SelectedDimensions, 2))
    {
      auto topoView = make_uniform<2>::view(topo);
      const std::string shape("quad");
      func(shape, topoView);
    }
    break;
  case 1:
    if constexpr(dimension_selected(SelectedDimensions, 3))
    {
      auto topoView = make_uniform<1>::view(topo);
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
